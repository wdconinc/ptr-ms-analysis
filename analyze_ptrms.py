#!/usr/bin/env python3
"""
analyze_ptrms.py — PTR-MS data analysis for Ionikon HDF5 files.

Usage:
    python analyze_ptrms.py <file.h5> [options]

Examples:
    python analyze_ptrms.py Data_*.h5
    python analyze_ptrms.py Data_*.h5 --min-conc 0.5
    python analyze_ptrms.py Data_*.h5 --compounds "CH4O" "C2H6O" --plot
    python analyze_ptrms.py Data_*.h5 --avg-spectrum --plot
"""

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# HDF5 reading
# ---------------------------------------------------------------------------

LABVIEW_EPOCH_OFFSET = datetime(1904, 1, 1, tzinfo=timezone.utc).timestamp()


def decode(arr):
    """Decode a byte-string numpy array to a list of Python strings."""
    return [x.decode("latin-1").strip() for x in arr]


def read_file(path: Path) -> dict:
    """Read all relevant data from an Ionikon PTR-MS HDF5 file."""
    data = {}
    with h5py.File(path, "r") as f:
        # ── Timestamps ──────────────────────────────────────────────────────
        times_raw = f["SPECdata/Times"][:]          # shape (N, 4)
        # Column 2 is LabVIEW absolute seconds since 1904-01-01
        lv_ts = times_raw[:, 2]
        unix_ts = lv_ts + LABVIEW_EPOCH_OFFSET
        data["timestamps"] = pd.to_datetime(unix_ts, unit="s", utc=True)
        data["elapsed_s"] = unix_ts - unix_ts[0]

        # ── Trace data ──────────────────────────────────────────────────────
        info = f["TRACEdata/TraceInfo"][:]           # shape (8, N_masses)
        data["trace_names"]  = np.array(decode(info[1]))
        data["trace_masses"] = np.array([float(x) for x in decode(info[2])])
        data["trace_mass_lo"] = np.array([float(x) for x in decode(info[3])])
        data["trace_mass_hi"] = np.array([float(x) for x in decode(info[4])])
        data["conc"]  = f["TRACEdata/TraceConcentration"][:]   # ppbv
        data["raw"]   = f["TRACEdata/TraceRaw"][:]             # cps
        data["corr"]  = f["TRACEdata/TraceCorrected"][:]       # corrected cps

        # ── Instrument / reaction conditions ────────────────────────────────
        rxn_info = f["AddTraces/PTR-Reaction/Info"][:]
        rxn_data = f["AddTraces/PTR-Reaction/Data"][:]
        rxn_names = decode(rxn_info[0])
        rxn_units = decode(rxn_info[1])
        data["rxn"] = pd.DataFrame(
            rxn_data,
            columns=[f"{n} [{u}]" for n, u in zip(rxn_names, rxn_units)],
            index=data["timestamps"],
        )

        # ── Primary ion diagnostics (CalcTraces) ────────────────────────────
        ci = f["CalcTraces/Info"][:]
        cd = f["CalcTraces/Data"][:]
        calc_names = decode(ci[0])
        calc_units = decode(ci[1])
        data["calc"] = pd.DataFrame(
            cd,
            columns=[f"{n} [{u}]" for n, u in zip(calc_names, calc_units)],
            index=data["timestamps"],
        )

        # ── Average spectrum ─────────────────────────────────────────────────
        data["avg_spectrum"] = f["SPECdata/AverageSpec"][:]
        data["cal_mapping"]  = f["CALdata/Mapping"][:]   # mass calibration (a, b)

    return data


# ---------------------------------------------------------------------------
# Mass calibration
# ---------------------------------------------------------------------------

def calibrated_mz(n_bins: int, mapping: np.ndarray) -> np.ndarray:
    """Convert spectrum bin index to m/z using the stored calibration points.

    Ionikon stores calibration as (mass, bin) reference pairs in CALdata/Mapping.
    The TOF sqrt-law gives:  sqrt(mz) = a * bin + b
    We fit a and b by least-squares from the reference points, then evaluate
    mz(bin) = (a * bin + b)^2  for every bin in the spectrum.
    """
    ref_mass = mapping[:, 0]
    ref_bin  = mapping[:, 1]
    # Fit sqrt(mass) = a * bin + b
    a, b = np.polyfit(ref_bin, np.sqrt(ref_mass), 1)
    bins = np.arange(n_bins)
    return (a * bins + b) ** 2


# ---------------------------------------------------------------------------
# Peak list
# ---------------------------------------------------------------------------

KNOWN_PRIMARY_IONS = {
    "H3O 18+", "O[18O]+", "O2+", "(H2O)H+", "cluster h2oH+",
    "H5O[18O]+", "NO+", "(H2O)+", "(H3N)+", "NO+ Isotope",
    "(N2)+", "(N2)H+", "H3O+",
}


def build_peak_list(data: dict, min_conc: float = 0.0,
                    exclude_primary: bool = True) -> pd.DataFrame:
    """Build a peak list DataFrame sorted by mean concentration."""
    conc = data["conc"]
    raw  = data["raw"]
    names  = data["trace_names"]
    masses = data["trace_masses"]

    df = pd.DataFrame({
        "name":          names,
        "mz":            masses,
        "mz_lo":         data["trace_mass_lo"],
        "mz_hi":         data["trace_mass_hi"],
        "mean_ppbv":     conc.mean(axis=0),
        "std_ppbv":      conc.std(axis=0),
        "min_ppbv":      conc.min(axis=0),
        "max_ppbv":      conc.max(axis=0),
        "mean_raw_cps":  raw.mean(axis=0),
    }).sort_values("mean_ppbv", ascending=False)

    if exclude_primary:
        df = df[~df["name"].isin(KNOWN_PRIMARY_IONS)]
    if min_conc > 0:
        df = df[df["mean_ppbv"] >= min_conc]

    return df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Instrument summary
# ---------------------------------------------------------------------------

def instrument_summary(data: dict) -> dict:
    """Return a dict of key instrument parameters."""
    rxn = data["rxn"]
    calc = data["calc"]
    t0 = data["timestamps"][0]
    t1 = data["timestamps"][-1]
    duration_min = (t1 - t0).total_seconds() / 60

    summary = {
        "start_time":     t0.isoformat(),
        "end_time":       t1.isoformat(),
        "duration_min":   round(duration_min, 2),
        "n_spectra":      len(data["timestamps"]),
        "n_masses":       len(data["trace_names"]),
    }

    for col in rxn.columns:
        summary[col] = round(float(rxn[col].mean()), 3)
    for col in calc.columns:
        summary[col] = round(float(calc[col].mean()), 3)

    return summary


# ---------------------------------------------------------------------------
# Time series helpers
# ---------------------------------------------------------------------------

def get_traces(data: dict, patterns: list[str]) -> pd.DataFrame:
    """Return a DataFrame of concentration time series matching name patterns."""
    names = data["trace_names"]
    ts    = data["timestamps"]
    conc  = data["conc"]

    cols = {}
    for pat in patterns:
        pat_lower = pat.lower()
        matches = [i for i, n in enumerate(names) if pat_lower in n.lower()]
        if not matches:
            print(f"  [warn] no trace matching '{pat}'", file=sys.stderr)
        for i in matches:
            cols[f"{names[i]} (m/z {data['trace_masses'][i]:.3f})"] = conc[:, i]

    return pd.DataFrame(cols, index=ts)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_peak_list(peaks: pd.DataFrame, title: str, top_n: int = 30,
                   output: Path | None = None):
    import matplotlib.pyplot as plt

    top = peaks.head(top_n)
    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(top))
    ax.bar(x, top["mean_ppbv"], yerr=top["std_ppbv"], capsize=3,
           color="steelblue", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{r['name']}\n{r['mz']:.2f}" for _, r in top.iterrows()],
        rotation=45, ha="right", fontsize=7,
    )
    ax.set_ylabel("Mean concentration (ppbv)")
    ax.set_yscale("log")
    ax.set_title(title)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    _save_or_show(fig, output)


def plot_time_series(traces_df: pd.DataFrame, title: str,
                     output: Path | None = None):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(12, 5))
    for col in traces_df.columns:
        elapsed = (traces_df.index - traces_df.index[0]).total_seconds() / 60
        ax.plot(elapsed, traces_df[col], label=col, linewidth=1.2)
    ax.set_xlabel("Elapsed time (min)")
    ax.set_ylabel("Concentration (ppbv)")
    ax.set_title(title)
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    _save_or_show(fig, output)


def plot_avg_spectrum(data: dict, title: str, mz_range: tuple | None = None,
                      output: Path | None = None):
    import matplotlib.pyplot as plt

    spec = data["avg_spectrum"]
    mz   = calibrated_mz(len(spec), data["cal_mapping"])

    fig, ax = plt.subplots(figsize=(14, 5))
    if mz_range:
        mask = (mz >= mz_range[0]) & (mz <= mz_range[1])
        ax.plot(mz[mask], spec[mask], linewidth=0.6, color="steelblue")
        ax.set_xlim(mz_range)
    else:
        ax.plot(mz, spec, linewidth=0.4, color="steelblue")
    ax.set_xlabel("m/z")
    ax.set_ylabel("Signal (counts)")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    _save_or_show(fig, output)


def _save_or_show(fig, output: Path | None):
    import matplotlib.pyplot as plt
    if output:
        fig.savefig(output, dpi=150)
        print(f"  Saved: {output}")
        plt.close(fig)
    else:
        plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Analyze Ionikon PTR-MS HDF5 files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("files", nargs="+", type=Path, help="HDF5 file(s) to analyze")
    p.add_argument("--min-conc", type=float, default=0.1, metavar="PPB",
                   help="Minimum mean concentration to include in peak list (default: 0.1 ppbv)")
    p.add_argument("--top-n", type=int, default=30,
                   help="Number of top peaks to show/plot (default: 30)")
    p.add_argument("--include-primary", action="store_true",
                   help="Include known primary ions in peak list")
    p.add_argument("--compounds", nargs="+", metavar="NAME",
                   help="Plot time series for compounds matching these name substrings")
    p.add_argument("--avg-spectrum", action="store_true",
                   help="Plot the average mass spectrum")
    p.add_argument("--mz-range", nargs=2, type=float, metavar=("MZ_LO", "MZ_HI"),
                   help="m/z range for average spectrum plot")
    p.add_argument("--plot", action="store_true",
                   help="Show interactive plots (default: save PNG files)")
    p.add_argument("--no-csv", action="store_true",
                   help="Skip writing peak list CSV")
    return p.parse_args()


def analyze_file(path: Path, args):
    print(f"\n{'='*60}")
    print(f"File: {path.name}")
    print(f"{'='*60}")

    data = read_file(path)
    stem = path.stem

    # ── Instrument summary ───────────────────────────────────────────────
    summ = instrument_summary(data)
    print(f"  Start:    {summ['start_time']}")
    print(f"  Duration: {summ['duration_min']} min  ({summ['n_spectra']} spectra)")
    print(f"  E/N:      {summ.get('E/N_Act [Td]', 'n/a')} Td")
    print(f"  H3O+ purity: {summ.get('H3O+ [%]', 'n/a')} %")
    print(f"  O2+:         {summ.get('O2+ [%]', 'n/a')} %")
    print(f"  PI total:    {summ.get('PI (total) [x1E6]', 'n/a')} ×10⁶ cps")

    # ── Peak list ────────────────────────────────────────────────────────
    peaks = build_peak_list(
        data,
        min_conc=args.min_conc,
        exclude_primary=not args.include_primary,
    )
    print(f"\n  Peak list (> {args.min_conc} ppbv): {len(peaks)} compounds")
    print(
        peaks[["name", "mz", "mean_ppbv", "std_ppbv", "max_ppbv"]]
        .head(args.top_n)
        .to_string(index=False, float_format="{:.3f}".format)
    )

    if not args.no_csv:
        csv_path = path.parent / f"{stem}_peak_list.csv"
        peaks.to_csv(csv_path, index=False)
        print(f"\n  Peak list saved: {csv_path}")

    # ── Peak list bar chart ──────────────────────────────────────────────
    out = None if args.plot else path.parent / f"{stem}_peaks.png"
    plot_peak_list(peaks, title=f"Peak list — {stem}", top_n=args.top_n, output=out)

    # ── Optional time series ─────────────────────────────────────────────
    if args.compounds:
        traces = get_traces(data, args.compounds)
        if not traces.empty:
            out = None if args.plot else path.parent / f"{stem}_timeseries.png"
            plot_time_series(traces, title=f"Time series — {stem}", output=out)

    # ── Optional average spectrum ────────────────────────────────────────
    if args.avg_spectrum:
        mz_range = tuple(args.mz_range) if args.mz_range else None
        out = None if args.plot else path.parent / f"{stem}_spectrum.png"
        plot_avg_spectrum(data, title=f"Avg spectrum — {stem}",
                          mz_range=mz_range, output=out)

    return peaks, summ


def main():
    args = parse_args()
    all_peaks = []

    for path in args.files:
        if not path.exists():
            print(f"[error] file not found: {path}", file=sys.stderr)
            continue
        peaks, summ = analyze_file(path, args)
        peaks.insert(0, "file", path.name)
        all_peaks.append(peaks)

    # Multi-file combined peak list
    if len(all_peaks) > 1:
        combined = pd.concat(all_peaks, ignore_index=True)
        out = args.files[0].parent / "combined_peak_list.csv"
        combined.to_csv(out, index=False)
        print(f"\nCombined peak list saved: {out}")


if __name__ == "__main__":
    main()
