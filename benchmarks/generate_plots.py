#!/usr/bin/env python3
"""Generate PNG visualization plots from benchmark and convergence CSV results.

Reads CSV files from benchmarks/results/ and produces:
  1. speedup_comparison.png  — Bar chart: R vs Rust timing across scales and configs
  2. scaling_curves.png      — Log-log timing vs dataset size with scaling exponents
  3. convergence_bias.png    — Bias vs temperature with 95% CI bands
  4. ci_coverage.png         — Bar chart of CI coverage (%) per approach and config
  5. memory_scaling.png      — Estimated memory usage vs scale for R vs Rust crossbasis

All output saved to benchmarks/results/report/.
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")
REPORT_DIR = os.path.join(RESULTS_DIR, "report")

SPEEDUP_CSV = os.path.join(RESULTS_DIR, "speedup_comparison.csv")
OPTIMIZED_CSV = os.path.join(RESULTS_DIR, "optimized_timing_results.csv")
TIMING_CSV = os.path.join(RESULTS_DIR, "timing_results.csv")
CONVERGENCE_METRICS_CSV = os.path.join(RESULTS_DIR, "convergence_metrics.csv")
CONVERGENCE_DETAILS_CSV = os.path.join(RESULTS_DIR, "convergence_details.csv")

# ---------------------------------------------------------------------------
# Style helpers
# ---------------------------------------------------------------------------
PALETTE = {
    "R": "#4C72B0",       # steel blue
    "Rust": "#DD8452",    # warm orange
    "single_stage": "#4C72B0",
    "two_stage": "#C44E52",
}

SCALE_ORDER = ["10mb", "100mb", "1gb"]
SCALE_LABELS = {"10mb": "10 MB", "100mb": "100 MB", "1gb": "1 GB"}
SCALE_ROWS = {"10mb": 122736, "100mb": 1201790, "1gb": 12017900, "10gb": 120179000}

STAGE_ORDER = ["crossbasis", "glm", "crosspred", "crossreduce"]
STAGE_LABELS = {
    "crossbasis": "Cross-basis",
    "glm": "GLM",
    "crosspred": "Prediction",
    "crossreduce": "Reduction",
}


def _apply_academic_style(ax):
    """Apply a clean academic style: remove top/right spines, light grid."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=4)
    ax.grid(axis="y", linewidth=0.3, alpha=0.5, color="#cccccc")


def _ensure_report_dir():
    os.makedirs(REPORT_DIR, exist_ok=True)


# ===================================================================
# Plot 1: Speedup Comparison (Bar chart)
# ===================================================================
def plot_speedup_comparison():
    """Bar chart: R vs Rust timing across scales (10MB, 100MB, 1GB) and configs,
    grouped by stage."""
    df = pd.read_csv(SPEEDUP_CSV)

    # Determine available configs and order them
    configs = sorted(df["config"].unique())
    stages = [s for s in STAGE_ORDER if s in df["stage"].unique()]

    fig, axes = plt.subplots(1, len(stages), figsize=(4.5 * len(stages), 5),
                             sharey=False)
    if len(stages) == 1:
        axes = [axes]

    for ax, stage in zip(axes, stages):
        sub = df[df["stage"] == stage].copy()
        # Build combined label: scale + config
        sub["label"] = sub.apply(
            lambda r: f"{SCALE_LABELS.get(r['scale'], r['scale'])}\n{r['config']}",
            axis=1,
        )
        # Sort by scale order then config
        sub["scale_idx"] = sub["scale"].map(
            {s: i for i, s in enumerate(SCALE_ORDER)}
        )
        sub = sub.sort_values(["scale_idx", "config"]).reset_index(drop=True)

        x = np.arange(len(sub))
        w = 0.35

        ax.bar(x - w / 2, sub["r_time_sec"], w, label="R (baseline)",
               color=PALETTE["R"], edgecolor="white", linewidth=0.5)
        ax.bar(x + w / 2, sub["rust_time_sec"], w, label="Rust (optimized)",
               color=PALETTE["Rust"], edgecolor="white", linewidth=0.5)

        # Annotate speedup factor
        for i, row in sub.iterrows():
            idx = sub.index.get_loc(i)
            speedup = row["speedup_factor"]
            ymax = max(row["r_time_sec"], row["rust_time_sec"])
            if speedup > 1:
                ax.text(idx, ymax * 1.05, f"{speedup:.1f}x",
                        ha="center", va="bottom", fontsize=7, fontweight="bold",
                        color="#2a7f2a")
            else:
                ax.text(idx, ymax * 1.05, f"{speedup:.2f}x",
                        ha="center", va="bottom", fontsize=7, color="#888888")

        ax.set_xticks(x)
        ax.set_xticklabels(sub["label"], fontsize=8)
        ax.set_ylabel("Time (seconds)" if ax == axes[0] else "")
        ax.set_title(STAGE_LABELS.get(stage, stage), fontsize=11, fontweight="bold")
        _apply_academic_style(ax)

    # Single legend for the figure
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, 1.02))

    fig.suptitle("R vs Rust Timing Comparison by Stage", fontsize=13,
                 fontweight="bold", y=1.08)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    path = os.path.join(REPORT_DIR, "speedup_comparison.png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {path}")


# ===================================================================
# Plot 2: Scaling Curves (Log-log)
# ===================================================================
def plot_scaling_curves():
    """Log-log plot: timing vs dataset size for each stage, R and Rust lines
    with scaling exponents annotated."""
    df_r = pd.read_csv(TIMING_CSV)
    df_r["backend"] = "R"
    df_rust = pd.read_csv(OPTIMIZED_CSV)
    df_rust["backend"] = "Rust"
    df = pd.concat([df_r, df_rust], ignore_index=True)

    # Focus on crossbasis stage for the clearest scaling story,
    # but also show glm if available
    stages_to_plot = [s for s in ["crossbasis", "glm"] if s in df["stage"].unique()]
    configs_to_plot = sorted(
        set(df[df["scale"].isin(SCALE_ORDER)]["config"].unique()) &
        set(["C1", "C2", "C3"])
    )

    fig, axes = plt.subplots(1, len(stages_to_plot),
                             figsize=(6 * len(stages_to_plot), 5))
    if len(stages_to_plot) == 1:
        axes = [axes]

    markers = {"C1": "o", "C2": "s", "C3": "^", "C4": "D", "C5": "v"}

    for ax, stage in zip(axes, stages_to_plot):
        for cfg in configs_to_plot:
            for backend, color, ls in [("R", PALETTE["R"], "-"),
                                       ("Rust", PALETTE["Rust"], "--")]:
                sub = df[(df["stage"] == stage) &
                         (df["config"] == cfg) &
                         (df["backend"] == backend) &
                         (df["scale"].isin(SCALE_ORDER))].copy()
                if sub.empty:
                    continue
                sub["n"] = sub["scale"].map(SCALE_ROWS)
                sub = sub.sort_values("n")
                ax.plot(sub["n"], sub["median_time_sec"],
                        marker=markers.get(cfg, "o"), markersize=6,
                        linestyle=ls, color=color, linewidth=1.5,
                        label=f"{backend} {cfg}")

                # Compute and annotate scaling exponent if 2+ points
                if len(sub) >= 2:
                    log_n = np.log10(sub["n"].values.astype(float))
                    log_t = np.log10(sub["median_time_sec"].values.astype(float))
                    # Use least-squares for exponent
                    if len(log_n) >= 2:
                        slope, _ = np.polyfit(log_n, log_t, 1)
                        # Annotate near the last point
                        last_n = sub["n"].iloc[-1]
                        last_t = sub["median_time_sec"].iloc[-1]
                        ax.annotate(
                            f"  a={slope:.2f}",
                            xy=(last_n, last_t),
                            fontsize=7, color=color, alpha=0.8,
                        )

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Number of Rows", fontsize=10)
        ax.set_ylabel("Time (seconds)" if ax == axes[0] else "")
        ax.set_title(f"{STAGE_LABELS.get(stage, stage)} Scaling",
                     fontsize=11, fontweight="bold")
        _apply_academic_style(ax)

        # Format x ticks nicely
        ax.xaxis.set_major_formatter(FuncFormatter(
            lambda x, _: f"{x/1e6:.0f}M" if x >= 1e6 else f"{x/1e3:.0f}K"
        ))

    # Legend below
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False,
               fontsize=8, bbox_to_anchor=(0.5, -0.05))

    fig.suptitle("Timing Scaling Curves (log-log)", fontsize=13,
                 fontweight="bold", y=1.02)
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])

    path = os.path.join(REPORT_DIR, "scaling_curves.png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {path}")


# ===================================================================
# Plot 3: Convergence Bias
# ===================================================================
def plot_convergence_bias():
    """Line plot: estimated RR vs temperature for single-stage and two-stage,
    overlaid with true RR. 95% CI bands shown."""
    df = pd.read_csv(CONVERGENCE_DETAILS_CSV)

    # Get unique scale values
    scales = sorted(df["scale"].unique())

    fig, axes = plt.subplots(1, len(scales), figsize=(7 * len(scales), 5),
                             sharey=True)
    if len(scales) == 1:
        axes = [axes]

    approach_labels = {
        "single_stage": "Single-Stage Optimized DLNM",
        "two_stage": "Two-Stage DLNM + Mixmeta",
    }

    for ax, scale in zip(axes, scales):
        sub = df[df["scale"] == scale]

        # Plot true RR (take from any approach — it's the same)
        for approach in ["single_stage", "two_stage"]:
            asub = sub[sub["approach"] == approach].sort_values("temperature")
            if asub.empty:
                continue

            color = PALETTE[approach]
            label = approach_labels.get(approach, approach)
            temp = asub["temperature"].values
            rr = asub["estimated_rr"].values
            ci_lo = asub["ci_lower"].values
            ci_hi = asub["ci_upper"].values

            ax.plot(temp, rr, color=color, linewidth=1.8, label=label)
            ax.fill_between(temp, ci_lo, ci_hi, color=color, alpha=0.15)

        # True RR — use data from the first available approach
        true_sub = sub[sub["approach"] == sub["approach"].iloc[0]].sort_values(
            "temperature"
        )
        ax.plot(
            true_sub["temperature"], true_sub["true_rr"],
            color="#333333", linewidth=1.5, linestyle=":", label="True RR",
        )

        ax.axhline(1.0, color="#999999", linewidth=0.5, linestyle="-")
        ax.set_xlabel("Temperature (°C)", fontsize=10)
        if ax == axes[0]:
            ax.set_ylabel("Relative Risk (RR)", fontsize=10)
        ax.set_title(f"Scale: {SCALE_LABELS.get(scale, scale)}",
                     fontsize=11, fontweight="bold")
        _apply_academic_style(ax)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, -0.06))

    fig.suptitle("Convergence: Estimated vs True Relative Risk",
                 fontsize=13, fontweight="bold", y=1.02)
    fig.tight_layout(rect=[0, 0.06, 1, 0.95])

    path = os.path.join(REPORT_DIR, "convergence_bias.png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {path}")


# ===================================================================
# Plot 4: CI Coverage
# ===================================================================
def plot_ci_coverage():
    """Bar chart: CI coverage (%) for each approach and config/scale."""
    df = pd.read_csv(CONVERGENCE_METRICS_CSV)

    fig, ax = plt.subplots(figsize=(8, 5))

    # Create labels combining config + scale
    df["label"] = df.apply(
        lambda r: f"{r['config']}\n{SCALE_LABELS.get(r['scale'], r['scale'])}",
        axis=1,
    )

    approaches = sorted(df["approach"].unique())
    n_approaches = len(approaches)
    x = np.arange(len(df["label"].unique()))
    labels_unique = df.drop_duplicates("label")["label"].values
    w = 0.35

    approach_labels = {
        "single_stage": "Single-Stage DLNM",
        "two_stage": "Two-Stage Mixmeta",
    }

    for i, approach in enumerate(approaches):
        sub = df[df["approach"] == approach]
        # Align bars with labels
        coverage_vals = []
        for lbl in labels_unique:
            match = sub[sub["label"] == lbl]
            coverage_vals.append(
                match["ci_coverage"].values[0] * 100 if len(match) > 0 else 0
            )
        offset = (i - (n_approaches - 1) / 2) * w
        bars = ax.bar(x + offset, coverage_vals, w,
                      label=approach_labels.get(approach, approach),
                      color=PALETTE.get(approach, f"C{i}"),
                      edgecolor="white", linewidth=0.5)
        # Annotate bar values
        for bar_rect, val in zip(bars, coverage_vals):
            if val > 0:
                ax.text(bar_rect.get_x() + bar_rect.get_width() / 2,
                        bar_rect.get_height() + 1.5,
                        f"{val:.0f}%", ha="center", va="bottom",
                        fontsize=8, fontweight="bold")

    # Reference line at 95%
    ax.axhline(95, color="#999999", linewidth=0.8, linestyle="--", label="Nominal 95%")

    ax.set_xticks(x)
    ax.set_xticklabels(labels_unique, fontsize=9)
    ax.set_ylabel("CI Coverage (%)", fontsize=10)
    ax.set_ylim(0, 105)
    ax.set_title("95% CI Coverage by Approach and Scale",
                 fontsize=13, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    _apply_academic_style(ax)

    fig.tight_layout()
    path = os.path.join(REPORT_DIR, "ci_coverage.png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {path}")


# ===================================================================
# Plot 5: Memory Scaling (Estimated)
# ===================================================================
def plot_memory_scaling():
    """Line plot: estimated memory usage vs scale for R vs Rust crossbasis.

    Memory estimates:
    - R crossbasis: materializes full lag matrix per basis column
        mem_R = n_rows * cb_cols * lag_range * 8 bytes  (full lag matrix)
    - Rust crossbasis: fused kernel, no lag matrix materialization
        mem_Rust = n_rows * cb_cols * 8 bytes  (output only)

    We use representative values: cb_cols from C2 config (20 cols, lag 0-21 => 22 lags).
    """
    # Read optimized timing for n_rows / cb_cols info
    df = pd.read_csv(OPTIMIZED_CSV)
    cb = df[df["stage"] == "crossbasis"].copy()

    # Use C2 config (representative, has data across scales)
    c2 = cb[cb["config"] == "C2"]
    if c2.empty:
        c2 = cb.drop_duplicates("scale")  # fallback to whatever is available

    lag_range = 22  # lag 0-21 for C2

    scales = []
    mem_r = []
    mem_rust = []

    for _, row in c2.iterrows():
        n = row["n_rows"]
        cols = row["cb_cols"]
        # R: materializes n x lag_range intermediate per basis col,
        # total ~ n_rows * lag_range * n_var_basis_cols * 8 bytes
        # For C2: n_var_basis = 5 (ns(5)), so intermediate = n * 22 * 5 * 8
        # Plus output matrix: n * 20 * 8
        n_var = max(1, int(cols / (lag_range - 1 + 1))) if lag_range > 0 else 1
        r_mem = (n * lag_range * n_var * 8) + (n * cols * 8)
        # Rust: only output matrix (fused kernel, no lag matrix)
        rust_mem = n * cols * 8
        scales.append(row["scale"])
        mem_r.append(r_mem / 1e9)       # GB
        mem_rust.append(rust_mem / 1e9)  # GB

    # Also add 10GB estimate if not present
    if "10gb" not in scales:
        n_10gb = 120179000
        cols_10gb = 20
        n_var_10gb = 5
        r_mem_10gb = (n_10gb * lag_range * n_var_10gb * 8) + (n_10gb * cols_10gb * 8)
        rust_mem_10gb = n_10gb * cols_10gb * 8
        scales.append("10gb")
        mem_r.append(r_mem_10gb / 1e9)
        mem_rust.append(rust_mem_10gb / 1e9)

    # Sort by scale
    scale_sort = {"10mb": 0, "100mb": 1, "1gb": 2, "10gb": 3}
    order = sorted(range(len(scales)), key=lambda i: scale_sort.get(scales[i], 99))
    scales = [scales[i] for i in order]
    mem_r = [mem_r[i] for i in order]
    mem_rust = [mem_rust[i] for i in order]

    n_rows_list = [SCALE_ROWS.get(s, 0) for s in scales]

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(n_rows_list, mem_r, marker="o", markersize=7, linewidth=2,
            color=PALETTE["R"], label="R (with lag matrix)")
    ax.plot(n_rows_list, mem_rust, marker="s", markersize=7, linewidth=2,
            color=PALETTE["Rust"], label="Rust (fused, output only)")

    # Annotate values
    for i, s in enumerate(scales):
        ax.annotate(f"  {mem_r[i]:.2f} GB",
                    xy=(n_rows_list[i], mem_r[i]),
                    fontsize=7, color=PALETTE["R"])
        ax.annotate(f"  {mem_rust[i]:.2f} GB",
                    xy=(n_rows_list[i], mem_rust[i]),
                    fontsize=7, color=PALETTE["Rust"])

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.xaxis.set_major_formatter(FuncFormatter(
        lambda x, _: f"{x/1e6:.0f}M" if x >= 1e6 else f"{x/1e3:.0f}K"
    ))
    ax.set_xlabel("Number of Rows", fontsize=10)
    ax.set_ylabel("Estimated Memory (GB)", fontsize=10)
    ax.set_title("Crossbasis Memory Usage: R vs Rust (C2 Config, lag 0-21)",
                 fontsize=12, fontweight="bold")
    ax.legend(frameon=False, fontsize=10)
    _apply_academic_style(ax)

    fig.tight_layout()
    path = os.path.join(REPORT_DIR, "memory_scaling.png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {path}")


# ===================================================================
# Main
# ===================================================================
def main():
    _ensure_report_dir()
    print("Generating plots...")

    print("\n[1/5] Speedup comparison")
    plot_speedup_comparison()

    print("[2/5] Scaling curves")
    plot_scaling_curves()

    print("[3/5] Convergence bias")
    plot_convergence_bias()

    print("[4/5] CI coverage")
    plot_ci_coverage()

    print("[5/5] Memory scaling")
    plot_memory_scaling()

    print("\nAll plots saved to", REPORT_DIR)


if __name__ == "__main__":
    main()
