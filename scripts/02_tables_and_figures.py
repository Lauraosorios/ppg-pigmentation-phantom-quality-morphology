"""
02_tables_and_figures.py -- every table and figure of the article from the
recording-level tables produced by step 01.

Tables (results/tables_manuscript/):
    Table3_moderate_large_effects.csv
    S1_Table_feature_classification.csv
    S2_Table_JT_trend_results.csv
    S3_Table_pale_to_dark_SMD_by_state.csv
Figures (results/figures/):
    Fig9_templates_dark.png, S1a_Fig_templates_pale.png, S1b_Fig_templates_medium.png
    Fig10_SQI_reflectance.png, Fig10_SQI_transmittance.png
    Fig11_effect_size_summary.png
    S2a_Fig_morphology_reflectance.png, S2b_Fig_morphology_transmittance.png

Input: results/tables/recording_{sqi,morphology}_table.csv, pulse_templates.csv
Does not need pyPPG or the raw data; runs in about a minute.

Usage:
    python scripts/02_tables_and_figures.py
"""

import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (RECORDING_SQI_CSV, RECORDING_MORPHOLOGY_CSV, TEMPLATES_CSV,   # noqa: E402
                    MANUSCRIPT_TABLES_DIR, FIGURES_DIR)
from ppgphantom.features import (MORPHOLOGY_FEATURES, MORPHOLOGY_FAMILIES,      # noqa: E402
                                 feature_table)
from ppgphantom import trend, figures                                          # noqa: E402

CLASS_ORDER = ["negligible/flat", "small", "moderate", "large"]


def print_summary(jt_df: pd.DataFrame, table3: pd.DataFrame) -> None:
    n = len(jt_df)
    print(f"\n{n} feature x channel combinations")
    for cls, c in jt_df["effect_size_class"].value_counts().reindex(CLASS_ORDER, fill_value=0).items():
        print(f"  {cls:16s} {c:4d}  ({100 * c / n:.1f}%)")
    print(f"Moderate or large: {len(table3)} combinations, {table3['feature'].nunique()} features\n")


def main():
    os.makedirs(MANUSCRIPT_TABLES_DIR, exist_ok=True)
    os.makedirs(FIGURES_DIR, exist_ok=True)

    features = feature_table()
    morph_df = pd.read_csv(RECORDING_MORPHOLOGY_CSV)
    sqi_df = pd.read_csv(RECORDING_SQI_CSV)
    templates = pd.read_csv(TEMPLATES_CSV)

    # ---- Tables ----------------------------------------------------------
    order = MORPHOLOGY_FAMILIES + ["SQI"]
    s1 = (features.assign(_o=features["feature_family"].map(order.index))
          .sort_values(["_o", "feature"])[["feature", "feature_family", "description"]])

    jt_df, state_df = trend.run_trend_analysis(features, morph_df, sqi_df)
    summary_df = trend.state_consistency_summary(state_df)
    table3 = trend.moderate_large_table(jt_df, summary_df)
    s3 = trend.smd_by_state_table(table3, state_df)

    s1.to_csv(os.path.join(MANUSCRIPT_TABLES_DIR, "S1_Table_feature_classification.csv"), index=False)
    jt_df[["feature", "feature_family", "mode", "wavelength", "N", "JT_statistic", "Z", "p_value",
           "r", "direction", "effect_size_class"]].round(6).to_csv(
        os.path.join(MANUSCRIPT_TABLES_DIR, "S2_Table_JT_trend_results.csv"), index=False)
    table3.round(6).to_csv(os.path.join(MANUSCRIPT_TABLES_DIR, "Table3_moderate_large_effects.csv"), index=False)
    s3.round(3).to_csv(os.path.join(MANUSCRIPT_TABLES_DIR, "S3_Table_pale_to_dark_SMD_by_state.csv"), index=False)
    print("Tables: Table 3, S1, S2, S3 saved")
    print_summary(jt_df, table3)

    # ---- Figures ---------------------------------------------------------
    print("Figures:")
    for skin, name in (("D", "Fig9_templates_dark"), ("P", "S1a_Fig_templates_pale"),
                       ("M", "S1b_Fig_templates_medium")):
        figures.plot_templates(templates, skin, os.path.join(FIGURES_DIR, f"{name}.png"))
    for mode in ("Reflectance", "Transmittance"):
        figures.plot_sqi_by_skin(sqi_df, mode, os.path.join(FIGURES_DIR, f"Fig10_SQI_{mode.lower()}.png"))
    figures.plot_effect_size_summary(jt_df, os.path.join(FIGURES_DIR, "Fig11_effect_size_summary.png"))
    for mode, tag in (("Reflectance", "S2a"), ("Transmittance", "S2b")):
        figures.plot_morphology_boxplots(morph_df, mode, MORPHOLOGY_FEATURES,
                                         os.path.join(FIGURES_DIR, f"{tag}_Fig_morphology_{mode.lower()}.png"))


if __name__ == "__main__":
    main()
