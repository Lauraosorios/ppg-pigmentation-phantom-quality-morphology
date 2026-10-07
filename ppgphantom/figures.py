"""
figures.py -- manuscript figures: pulse templates (Fig 9, S1 Fig), SQI
boxplots (Fig 10), effect-size summary (Fig 11) and morphology boxplots (S2 Fig).
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches   # noqa: E402
import matplotlib.pyplot as plt         # noqa: E402
import numpy as np                      # noqa: E402
import pandas as pd                     # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (CHANNELS, SKIN_TONES, SKIN_LABELS,      # noqa: E402
                    HEART_RATES, TARGET_FLOWS)
from ppgphantom.features import MORPHOLOGY_FAMILIES         # noqa: E402

SUMMARY_FAMILY_ORDER = ["SQI", "Amplitude/area", "Composite index", "Width/duration", "Timing"]
EFFECT_CLASS_ORDER = ["negligible/flat", "small", "moderate", "large"]
EFFECT_CLASS_COLOR = {"negligible/flat": "#d9d9d9", "small": "#fdbb84",
                      "moderate": "#e34a33", "large": "#b30000"}
WAVELENGTH_COLOR = {"RED": "#e03428", "IR": "#4a90d9", "GREEN": "#27a348"}
MODE_WAVELENGTHS = {"Reflectance": ["RED", "IR", "GREEN"], "Transmittance": ["RED", "IR"]}
FAMILY_TITLE_COLOR = {"Timing": "#1f4fd1", "Width/duration": "#1fa02a",
                      "Amplitude/area": "#e0662a", "Composite index": "#c020c0"}

FS_PANEL_TITLE, FS_AXIS_LABEL, FS_TICK = 20, 15, 13
FS_FAMILY_LABEL_A, FS_FAMILY_LABEL_B = 14, 15
FS_FEATURE_MAIN = 13
FS_LEGEND, FS_CBAR, FS_N_ANNOT = 14, 13, 13


def _family_bars(jt_df: pd.DataFrame, ax, legend_ax) -> None:
    """Panel A: proportion of feature x channel tests in each effect-size class, per family."""
    props = (jt_df.groupby(["feature_family", "effect_size_class"]).size()
             .unstack(fill_value=0)
             .reindex(index=SUMMARY_FAMILY_ORDER, columns=EFFECT_CLASS_ORDER, fill_value=0))
    props = props.div(props.sum(axis=1), axis=0)

    y = np.arange(len(SUMMARY_FAMILY_ORDER))
    left = np.zeros(len(SUMMARY_FAMILY_ORDER))
    handles = []
    for cls in EFFECT_CLASS_ORDER:
        vals = props[cls].to_numpy()
        handles.append(ax.barh(y, vals, left=left, color=EFFECT_CLASS_COLOR[cls], edgecolor="white",
                               linewidth=0.8, height=0.62, label=cls))
        left += vals

    ax.set_yticks(y)
    ax.set_yticklabels(SUMMARY_FAMILY_ORDER, fontsize=FS_FAMILY_LABEL_A)
    ax.set_xlim(0, 1.38)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.tick_params(axis="x", labelsize=FS_TICK)
    ax.set_xlabel("Proportion of feature × channel tests", fontsize=FS_AXIS_LABEL, labelpad=8)
    ax.invert_yaxis()
    ax.set_title("A", loc="left", fontsize=FS_PANEL_TITLE, fontweight="bold")
    n_per_family = jt_df.groupby("feature_family").size().reindex(SUMMARY_FAMILY_ORDER)
    for yi, fam in zip(y, SUMMARY_FAMILY_ORDER):
        ax.text(1.05, yi, f"N={n_per_family[fam]}", va="center", fontsize=FS_N_ANNOT, color="#333333")

    legend_ax.legend(handles, EFFECT_CLASS_ORDER, loc="center", ncol=2, fontsize=FS_LEGEND,
                     frameon=False, handlelength=1.6, columnspacing=1.4)
    legend_ax.axis("off")


def _effect_heatmap(jt_df: pd.DataFrame, features: list, ax, label_ax, feature_fontsize: float):
    """Panel B: JT r per feature (rows, grouped by family) x channel (columns)."""
    feat_family = jt_df.drop_duplicates("feature").set_index("feature")["feature_family"]

    ordered, bounds, idx = [], [], 0
    for fam in SUMMARY_FAMILY_ORDER:
        feats = sorted(f for f in features if feat_family.get(f) == fam)
        if not feats:
            continue
        ordered.extend(feats)
        bounds.append((idx, idx + len(feats), fam))
        idx += len(feats)

    piv = jt_df.pivot(index="feature", columns="channel", values="r").reindex(index=ordered, columns=CHANNELS)
    im = ax.imshow(piv.to_numpy(dtype=float), aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(CHANNELS)))
    ax.set_xticklabels([c.replace("Reflectance", "Refl.").replace("Transmittance", "Trans.") for c in CHANNELS],
                       fontsize=FS_TICK, rotation=30, ha="right")
    ax.set_yticks(range(len(ordered)))
    ax.set_yticklabels(ordered, fontsize=feature_fontsize)
    ax.set_title("B", loc="left", fontsize=FS_PANEL_TITLE, fontweight="bold")
    ax.set_xlim(-0.5, len(CHANNELS) - 0.5)
    for start, _, _ in bounds:
        if start > 0:
            ax.axhline(start - 0.5, color="black", linewidth=1.2)

    label_ax.set_ylim(ax.get_ylim())
    label_ax.set_xlim(0, 1)
    label_ax.axis("off")
    for start, end, fam in bounds:
        if start > 0:
            label_ax.axhline(start - 0.5, color="black", linewidth=1.2, xmin=0, xmax=0.15)
        label_ax.text(0.1, (start + end - 1) / 2, fam, va="center", ha="left",
                      fontsize=FS_FAMILY_LABEL_B, fontweight="bold")
    return im


def _top_colorbar(fig, im, cax) -> None:
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
    cax.xaxis.set_ticks_position("top")
    cax.xaxis.set_label_position("top")
    cbar.set_label("JT effect size  r", fontsize=FS_CBAR, labelpad=6)
    cbar.ax.tick_params(labelsize=FS_TICK)
    cbar.ax.annotate("decreasing", xy=(0, 1), xycoords="axes fraction", xytext=(0, 42),
                     textcoords="offset points", ha="left", va="bottom", fontsize=FS_TICK - 1, color="#1a3a7a")
    cbar.ax.annotate("increasing", xy=(1, 1), xycoords="axes fraction", xytext=(0, 42),
                     textcoords="offset points", ha="right", va="bottom", fontsize=FS_TICK - 1, color="#7a1a1a")


def plot_effect_size_summary(jt_df: pd.DataFrame, out_path: str) -> None:
    """Fig 11: (A) effect-size classes per family, (B) heatmap of the features
    reaching |r| >= 0.10 in at least one channel."""
    features = (jt_df.groupby("feature")["r"].apply(lambda s: s.abs().max())
                .loc[lambda s: s >= 0.10].index.tolist())
    n_rows = len(features)
    fig = plt.figure(figsize=(16, 1.3 + 0.42 * n_rows))
    outer = fig.add_gridspec(2, 2, height_ratios=[0.5, n_rows], width_ratios=[3.0, 6.4],
                             hspace=0.14, wspace=0.55)
    ax_legend = fig.add_subplot(outer[0, 0])
    ax_a = fig.add_subplot(outer[1, 0])
    inner = outer[:, 1].subgridspec(2, 2, height_ratios=[0.5, n_rows], width_ratios=[5, 1.3], wspace=0.1)
    ax_cbar = fig.add_subplot(inner[0, 0])
    ax_b = fig.add_subplot(inner[1, 0])
    ax_fam = fig.add_subplot(inner[1, 1])
    fig.add_subplot(inner[0, 1]).axis("off")

    _family_bars(jt_df, ax_a, ax_legend)
    im = _effect_heatmap(jt_df, features, ax_b, ax_fam, FS_FEATURE_MAIN)
    _top_colorbar(fig, im, ax_cbar)

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {os.path.basename(out_path)}  ({n_rows} features with |r| >= 0.10)")


SQI_DISPLAY_ORDER = ["TMCC", "Skewness", "Kurtosis", "Shannon Entropy", "Sample Entropy",
                     "AC/DC", "SNR", "PSD ratio", "Spec. Flatness"]
SQI_YLIMS = {"AC/DC": (0, None), "PSD ratio": (0.15, 0.50),
             "Shannon Entropy": (3.8, 5.0), "Sample Entropy": (0.07, 0.19)}


def plot_sqi_by_skin(sqi_df: pd.DataFrame, mode: str, out_path: str) -> None:
    """Fig 10: 3 x 3 grid of SQIs by skin tone, one box per wavelength. Each box
    is the 27 recording-level medians (3 replicates x 9 HR x flow states)."""
    from ppgphantom.sqi import SQI_LABELS, SQI_UNITS

    wavelengths = MODE_WAVELENGTHS[mode]
    step = len(wavelengths) + 1
    centers = [g * step + (len(wavelengths) - 1) / 2.0 for g in range(len(SKIN_TONES))]

    fig, axs = plt.subplots(3, 3, figsize=(13, 7.5))
    fig.suptitle(f"SQI — {mode}", fontsize=18, fontweight="bold")

    for ax, label in zip(axs.flat, SQI_DISPLAY_ORDER):
        for g_idx, skin in enumerate(SKIN_TONES):
            for w_idx, wl in enumerate(wavelengths):
                color = WAVELENGTH_COLOR[wl]
                vals = sqi_df.loc[(sqi_df["skin"] == skin) & (sqi_df["channel"] == f"{mode} {wl}"),
                                  label].dropna().to_numpy()
                bp = ax.boxplot([vals] if len(vals) else [[]], positions=[g_idx * step + w_idx], widths=0.72,
                                patch_artist=True, medianprops=dict(color="black", linewidth=1.5),
                                flierprops=dict(marker="o", markersize=3, alpha=0.4,
                                                markerfacecolor=color, markeredgecolor=color),
                                manage_ticks=False)
                for patch in bp["boxes"]:
                    patch.set_facecolor(color)
                    patch.set_alpha(0.75)

        ax.set_xticks(centers)
        ax.set_xticklabels([SKIN_LABELS[s] for s in SKIN_TONES], fontsize=13)
        ax.set_xlim(-0.8, len(SKIN_TONES) * step - 0.5)
        if label in SQI_YLIMS:
            lo, hi = SQI_YLIMS[label]
            cur = ax.get_ylim()
            ax.set_ylim(lo if lo is not None else cur[0], hi if hi is not None else cur[1])
        ax.set_title(f"{label}  {SQI_UNITS[SQI_LABELS.index(label)]}", fontsize=15, fontweight="bold")
        ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=3, min_n_ticks=3))
        ax.tick_params(axis="y", labelsize=12)
        ax.grid(True, axis="y", alpha=0.35, linewidth=0.5)

    fig.tight_layout(rect=[0, 0.08, 1, 0.93])
    fig.legend(handles=[mpatches.Patch(facecolor=WAVELENGTH_COLOR[wl], alpha=0.75, label=wl) for wl in wavelengths],
               loc="lower center", ncol=len(wavelengths), fontsize=14, frameon=True, bbox_to_anchor=(0.5, 0.0))
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"  Saved {os.path.basename(out_path)}")


def plot_templates(templates: pd.DataFrame, skin: str, out_path: str) -> None:
    """Fig 9 / S1 Fig: median pulse template with 5th-95th percentile envelope
    for one skin tone, every HR x target flow state; (a) reflectance, (b) transmittance."""
    t = templates[templates["skin"] == skin]
    n_hr, n_fl = len(HEART_RATES), len(TARGET_FLOWS)
    modes = ("Reflectance", "Transmittance")

    # Each mode is a block of (wavelengths x flows) columns, thin spacers between
    # wavelengths and a wider gap between the two modes.
    wl_spacer, mode_gap = 0.18, 0.45
    width_ratios, mode_start, col = [], {}, 0
    for mi, mode in enumerate(modes):
        mode_start[mode] = col
        n_wl = len(MODE_WAVELENGTHS[mode])
        for wi in range(n_wl):
            width_ratios.extend([1.0] * n_fl)
            col += n_fl
            if wi < n_wl - 1:
                width_ratios.append(wl_spacer)
                col += 1
        if mi == 0:
            width_ratios.append(mode_gap)
            col += 1

    fig = plt.figure(figsize=(0.66 * sum(width_ratios), 1.75 * n_hr))
    gs = fig.add_gridspec(n_hr, col, width_ratios=width_ratios, wspace=0.22, hspace=0.4,
                          top=0.78, bottom=0.1, left=0.06, right=0.99)

    for mode, panel in zip(modes, ("a)", "b)")):
        mode_axes = []
        for wi, wl in enumerate(MODE_WAVELENGTHS[mode]):
            color = WAVELENGTH_COLOR[wl]
            channel = f"{mode} {wl}"
            c0 = mode_start[mode] + wi * (n_fl + 1)
            top_axes = []
            for hi, hr in enumerate(HEART_RATES):
                for fj, flow in enumerate(TARGET_FLOWS):
                    ax = fig.add_subplot(gs[hi, c0 + fj])
                    ax.set_xticks([]); ax.set_yticks([])
                    mode_axes.append(ax)
                    if hi == 0:
                        top_axes.append(ax)
                    s = t[(t["hr"] == hr) & (t["co_lmin"] == flow) & (t["channel"] == channel)].sort_values("point")
                    n = int(s["n_pulses"].iloc[0]) if len(s) else 0
                    if len(s):
                        x = np.linspace(0, 1, len(s))
                        ax.fill_between(x, s["p05"], s["p95"], alpha=0.2, color=color)
                        ax.plot(x, s["median"], color=color, linewidth=1.1)
                    else:
                        ax.text(0.5, 0.5, "—", ha="center", va="center", fontsize=9, color="grey",
                                transform=ax.transAxes)
                    ax.set_title(f"{flow}L (n={n})", fontsize=5.8)
                    if mode == "Reflectance" and wi == 0 and fj == 0:
                        b = ax.get_position()
                        fig.text(0.02, (b.y0 + b.y1) / 2, f"{hr}\nbpm", ha="center", va="center",
                                 fontsize=9, fontweight="bold", rotation=90)
            x0, x1 = top_axes[0].get_position().x0, top_axes[-1].get_position().x1
            fig.text((x0 + x1) / 2, 0.865, wl, ha="center", fontsize=10.5, fontweight="bold", color=color)
        px0, px1 = mode_axes[0].get_position().x0, mode_axes[-1].get_position().x1
        fig.text((px0 + px1) / 2, 0.955, panel, ha="center", fontsize=15, fontweight="bold")

    fig.savefig(out_path, dpi=600)
    plt.close(fig)
    print(f"  Saved {os.path.basename(out_path)}")


def plot_morphology_boxplots(morph_df: pd.DataFrame, mode: str, features: dict, out_path: str) -> None:
    """S2 Fig: recording-level distribution of each morphology feature by skin
    tone, one box per wavelength. `features` maps feature -> (family, description)."""
    wavelengths = MODE_WAVELENGTHS[mode]
    df = morph_df[morph_df["channel"].isin([f"{mode} {wl}" for wl in wavelengths])]

    ordered, fam_of = [], {}
    for fam in MORPHOLOGY_FAMILIES:
        feats = sorted(f for f, (fm, _) in features.items() if fm == fam)
        ordered.extend(feats)
        fam_of.update({f: fam for f in feats})

    n_feat = len(ordered)
    n_cols = min(7, max(3, int(np.ceil(np.sqrt(n_feat)))))
    n_rows = int(np.ceil(n_feat / n_cols))
    step = len(wavelengths) + 1
    centers = [g * step + (len(wavelengths) - 1) / 2.0 for g in range(len(SKIN_TONES))]

    fig, axs = plt.subplots(n_rows, n_cols, figsize=(3.72 * n_cols, 2.5 * n_rows), squeeze=False)
    fig.suptitle(f"Morphology feature distributions — {mode}", fontsize=19, fontweight="bold", y=1.005)

    for i, ax in enumerate(axs.flat):
        if i >= n_feat:
            ax.set_visible(False)
            continue
        feat = ordered[i]
        for g_idx, skin in enumerate(SKIN_TONES):
            for w_idx, wl in enumerate(wavelengths):
                vals = df.loc[(df["skin"] == skin) & (df["channel"] == f"{mode} {wl}"), feat].dropna().to_numpy()
                bp = ax.boxplot([vals] if len(vals) else [[]], positions=[g_idx * step + w_idx], widths=0.8,
                                patch_artist=True, showfliers=True,
                                medianprops=dict(color="black", linewidth=1.0),
                                flierprops=dict(marker="o", markersize=1.5, alpha=0.4,
                                                markerfacecolor=WAVELENGTH_COLOR[wl],
                                                markeredgecolor=WAVELENGTH_COLOR[wl]),
                                manage_ticks=False)
                for patch in bp["boxes"]:
                    patch.set_facecolor(WAVELENGTH_COLOR[wl])
                    patch.set_alpha(0.78)
        ax.set_xticks(centers)
        ax.set_xticklabels(SKIN_TONES, fontsize=11.5)
        ax.set_xlim(-0.8, len(SKIN_TONES) * step - 0.5)
        ax.set_title(feat, fontsize=13, fontweight="bold", color=FAMILY_TITLE_COLOR[fam_of[feat]])
        ax.tick_params(axis="y", labelsize=10)
        ax.grid(True, axis="y", alpha=0.3, linewidth=0.5)

    handles = ([plt.Rectangle((0, 0), 1, 1, fc=WAVELENGTH_COLOR[wl], alpha=0.78) for wl in wavelengths] +
               [plt.Line2D([0], [0], color=FAMILY_TITLE_COLOR[c], lw=3) for c in MORPHOLOGY_FAMILIES])
    fig.legend(handles, wavelengths + MORPHOLOGY_FAMILIES, loc="lower center",
               ncol=len(wavelengths) + len(MORPHOLOGY_FAMILIES), fontsize=13, frameon=True,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=[0, 0.04, 1, 0.97])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {os.path.basename(out_path)}  ({n_feat} features)")
