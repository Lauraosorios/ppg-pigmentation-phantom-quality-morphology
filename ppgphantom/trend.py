"""
trend.py -- ordered pigmentation trend (pale -> medium -> dark) and its
consistency across the nine heart rate x target flow states.

Per feature x channel:
  1. each recording is already reduced to its median over pulses
     (recording-level tables);
  2. the three replicates of each skin x HR x flow condition are reduced to
     their median -> 27 observations (9 states x 3 skin tones);
  3. Jonckheere-Terpstra test with the order P < M < D; effect size r = Z / sqrt(N),
     classified negligible (<0.10), small (<0.30), moderate (<0.50), large (>=0.50);
  4. per state: sign of the pale -> dark change compared with the pooled
     direction, and the standardised mean difference
     SMD = (dark - pale) / SD of the 27 aggregated values.
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import CHANNELS, HEART_RATES, TARGET_FLOWS, SKIN_TONES   # noqa: E402

EFFECT_BINS = [(0.10, "negligible/flat"), (0.30, "small"), (0.50, "moderate")]


def jonckheere_terpstra(groups: list) -> tuple:
    """Returns (JT, E_JT, Var_JT, Z, two-sided p). Positive Z = increasing
    across the group order. Normal approximation, ties counted as 0.5,
    no tie correction to the variance."""
    groups = [np.asarray(g, dtype=float) for g in groups]
    groups = [g[np.isfinite(g)] for g in groups]
    groups = [g for g in groups if len(g) > 0]
    k = len(groups)
    if k < 2:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    n = np.array([len(g) for g in groups], dtype=float)
    N = int(n.sum())
    if N < 3:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    JT = 0.0
    for i in range(k - 1):
        for j in range(i + 1, k):
            diff = groups[j][:, None] - groups[i][None, :]
            JT += float(np.sum(diff > 0)) + 0.5 * float(np.sum(diff == 0))

    E_JT = (N ** 2 - float(np.sum(n ** 2))) / 4.0
    Var_JT = (N ** 2 * (2 * N + 3) - float(np.sum(n ** 2 * (2 * n + 3)))) / 72.0
    if Var_JT <= 0:
        return JT, E_JT, Var_JT, np.nan, np.nan

    z = (JT - E_JT) / float(np.sqrt(Var_JT))
    p = 2.0 * float(norm.sf(abs(z)))
    return JT, E_JT, Var_JT, float(z), float(p)


def effect_size_class(r: float) -> str:
    if not np.isfinite(r):
        return "n/a"
    for cutoff, label in EFFECT_BINS:
        if abs(r) < cutoff:
            return label
    return "large"


def _direction(r: float) -> str:
    if not np.isfinite(r) or r == 0:
        return "flat"
    return "increasing" if r > 0 else "decreasing"


def pivot_states(df: pd.DataFrame, feature: str, channel: str) -> pd.DataFrame:
    """Replicate medians: index (hr, co_lmin) -> columns P / M / D."""
    sub = df[df["channel"] == channel]
    med = sub.groupby(["skin", "hr", "co_lmin"])[feature].median().reset_index()
    piv = med.pivot(index=["hr", "co_lmin"], columns="skin", values=feature)
    for s in SKIN_TONES:
        if s not in piv.columns:
            piv[s] = np.nan
    return piv[SKIN_TONES].sort_index()


def run_trend_analysis(features: pd.DataFrame, morph_df: pd.DataFrame,
                       sqi_df: pd.DataFrame) -> tuple:
    """Returns (jt_df, state_df): one row per feature x channel, and one row
    per feature x channel x haemodynamic state."""
    jt_rows, state_rows = [], []

    for frow in features.itertuples(index=False):
        feat, fam = frow.feature, frow.feature_family
        df_src = morph_df if frow.source_table == "morphology" else sqi_df
        if feat not in df_src.columns:
            print(f"  [skip] {feat} not found in {frow.source_table} table")
            continue

        for channel in CHANNELS:
            if channel not in df_src["channel"].unique():
                continue
            piv = pivot_states(df_src, feat, channel)
            groups = [piv[s].dropna().to_numpy() for s in SKIN_TONES]
            JT, _, _, Z, p = jonckheere_terpstra(groups)
            N = sum(len(g) for g in groups)
            # r = Z / sqrt(N) is not strictly bounded by 1 (max ~1.04 for 3 x 9);
            # values are clipped to [-1, 1] as reported in the manuscript.
            r = float(np.clip(Z / np.sqrt(N), -1, 1)) if np.isfinite(Z) and N > 0 else np.nan
            direction = _direction(r)
            mode, wavelength = channel.split(" ", 1)

            jt_rows.append({
                "feature": feat, "feature_family": fam, "mode": mode, "wavelength": wavelength,
                "channel": channel, "N": N, "JT_statistic": JT, "Z": Z, "p_value": p, "r": r,
                "direction": direction, "effect_size_class": effect_size_class(r),
            })

            all_vals = np.concatenate(groups)
            sd27 = float(np.std(all_vals, ddof=1)) if len(all_vals) > 1 else np.nan

            for (hr, flow), row in piv.iterrows():
                pale, medium, dark = row["P"], row["M"], row["D"]
                if pd.isna(pale) or pd.isna(medium) or pd.isna(dark):
                    shape = "incomplete"
                elif pale < medium < dark:
                    shape = "increasing"
                elif pale > medium > dark:
                    shape = "decreasing"
                else:
                    shape = "non-monotonic"

                # Headline consistency measure: sign of the pale -> dark change.
                if pd.isna(pale) or pd.isna(dark):
                    pd_dir = "incomplete"
                elif dark > pale:
                    pd_dir = "increasing"
                elif dark < pale:
                    pd_dir = "decreasing"
                else:
                    pd_dir = "flat"

                smd = ((dark - pale) / sd27
                       if np.isfinite(sd27) and sd27 > 0 and not pd.isna(dark) and not pd.isna(pale)
                       else np.nan)

                state_rows.append({
                    "feature": feat, "feature_family": fam, "mode": mode, "wavelength": wavelength,
                    "channel": channel, "hr": hr, "co_lmin": flow,
                    "pale": pale, "medium": medium, "dark": dark,
                    "pale_dark_direction": pd_dir,
                    "agrees_pale_dark_with_overall": pd_dir in ("increasing", "decreasing") and pd_dir == direction,
                    "pale_to_dark_smd": smd,
                    "monotonic_shape": shape,
                    "agrees_monotonic_with_overall": shape in ("increasing", "decreasing") and shape == direction,
                    "overall_r": r, "overall_direction": direction,
                })

    return pd.DataFrame(jt_rows), pd.DataFrame(state_rows)


def state_consistency_summary(state_df: pd.DataFrame) -> pd.DataFrame:
    """One row per feature x channel: proportion of the 9 states agreeing with
    the pooled direction, and mean / min / max pale -> dark SMD."""
    rows = []
    for (feat, channel), g in state_df.groupby(["feature", "channel"], sort=False):
        eff = g["pale_to_dark_smd"].dropna()
        rows.append({
            "feature": feat, "feature_family": g["feature_family"].iloc[0],
            "mode": g["mode"].iloc[0], "wavelength": g["wavelength"].iloc[0], "channel": channel,
            "overall_r": g["overall_r"].iloc[0], "overall_direction": g["overall_direction"].iloc[0],
            "n_states": len(g),
            "n_pale_dark_increasing": int((g["pale_dark_direction"] == "increasing").sum()),
            "n_pale_dark_decreasing": int((g["pale_dark_direction"] == "decreasing").sum()),
            "proportion_agree_pale_dark": int(g["agrees_pale_dark_with_overall"].sum()) / 9,
            "mean_pale_to_dark_smd": eff.mean() if len(eff) else np.nan,
            "min_pale_to_dark_smd": eff.min() if len(eff) else np.nan,
            "max_pale_to_dark_smd": eff.max() if len(eff) else np.nan,
            "n_increasing_monotonic": int((g["monotonic_shape"] == "increasing").sum()),
            "n_decreasing_monotonic": int((g["monotonic_shape"] == "decreasing").sum()),
            "n_nonmonotonic": int((g["monotonic_shape"] == "non-monotonic").sum()),
            "n_incomplete": int((g["monotonic_shape"] == "incomplete").sum()),
            "proportion_agree_monotonic": int(g["agrees_monotonic_with_overall"].sum()) / 9,
        })
    return pd.DataFrame(rows)


def moderate_large_table(jt_df: pd.DataFrame, summary_df: pd.DataFrame) -> pd.DataFrame:
    """Manuscript Table 3: feature x channel combinations with |r| >= 0.30,
    sorted by |r|, with state agreement and mean (range) SMD."""
    t = jt_df[jt_df["effect_size_class"].isin(["moderate", "large"])].copy()
    t = t.assign(_abs_r=t["r"].abs()).sort_values("_abs_r", ascending=False).drop(columns="_abs_r")
    key = ["feature", "mode", "wavelength"]
    t = t.merge(summary_df[key + ["proportion_agree_pale_dark", "mean_pale_to_dark_smd",
                                  "min_pale_to_dark_smd", "max_pale_to_dark_smd"]],
                on=key, how="left").reset_index(drop=True)
    t["state_agreement"] = (t["proportion_agree_pale_dark"] * 9).round().astype(int).astype(str) + "/9"
    t["smd_display"] = t.apply(
        lambda row: f"{row.mean_pale_to_dark_smd:.2f} ({row.min_pale_to_dark_smd:.2f} to "
                    f"{row.max_pale_to_dark_smd:.2f})", axis=1)
    return t[["feature", "feature_family", "mode", "wavelength", "p_value", "r", "effect_size_class",
              "state_agreement", "mean_pale_to_dark_smd", "min_pale_to_dark_smd", "max_pale_to_dark_smd",
              "smd_display"]]


def smd_by_state_table(table3: pd.DataFrame, state_df: pd.DataFrame) -> pd.DataFrame:
    """Per-state pale -> dark SMD for the Table 3 combinations (S3 Table)."""
    key = ["feature", "mode", "wavelength"]
    state = state_df.merge(table3[key], on=key, how="inner").copy()
    state["state_label"] = state["hr"].astype(str) + "bpm_" + state["co_lmin"].astype(str) + "Lmin"
    piv = state.pivot(index=key, columns="state_label", values="pale_to_dark_smd")
    piv = piv[[f"{hr}bpm_{f}Lmin" for hr in HEART_RATES for f in TARGET_FLOWS]].reset_index()
    return table3[key + ["feature_family", "r", "p_value", "effect_size_class"]].merge(piv, on=key, how="left")
