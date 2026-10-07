# PPG signal quality and morphology across skin pigmentation levels: finger phantom analysis code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23210570.svg)](https://doi.org/10.5281/zenodo.23210570)

Code for the article *"Photoplethysmography (PPG) signal quality and morphology features across skin pigmentation levels using a multilayer vascular finger phantom"*, by L. Osorio-Sanchez, J. M. May and P. Kyriacou (City St George's, University of London).

The code takes the 81 phantom PPG recordings (3 skin tones × 3 heart rates × 3 target flows × 3 replicates) through the full pipeline:

1. AC / DC / noise separation and pyPPG fiducial detection.
2. Nine pulse-level signal quality indices (SQIs) and 38 pyPPG morphology features.
3. The Jonckheere–Terpstra pigmentation trend test (pale → medium → dark), plus its consistency across the nine haemodynamic states.

It reproduces the article's Figs 9-11, Table 3 and the supplementary figures and tables (S1-S2 Figs, S1-S3 Tables).

## Data

The recordings are on Zenodo: [https://doi.org/10.5281/zenodo.23208189](https://doi.org/10.5281/zenodo.23208189)

Download the dataset and either:

- unzip it into `data/` in this repository, so that you have `data/PPG/` and `data/Pressure/`; or
- set the environment variable `PPG_PHANTOM_DATA` to the folder that contains them.

The dataset `README.txt` describes the file format, units and acquisition details.

## Installation

pyPPG 1.0.73 requires **Python 3.10**.

```bash
python3.10 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Reproducing the results

| Step | Command | Needs | Produces |
|---|---|---|---|
| 1 | `python scripts/01_extract_features.py` | dataset, pyPPG; long run (81 recordings × 5 channels) | `results/tables/`: recording-level SQI and morphology tables, pulse templates, pulse counts |
| 2 | `python scripts/02_tables_and_figures.py` | step 1 output only; about a minute | every table and figure listed below |
| 3 | `python scripts/03_haemodynamics.py` | dataset (pressure) | MAP per recording and per state |

The step 1 outputs are included in `results/tables/`, so **step 2 can be run straight after cloning, without the dataset or a long pyPPG run.** The figures and tables themselves are not stored in the repository; running the scripts creates them in the locations below.

### Where each output appears in the article

| Article | File |
|---|---|
| Fig 9 (pulse templates, dark) | `results/figures/Fig9_templates_dark.png` |
| S1a / S1b Fig (pulse templates, pale / medium) | `results/figures/S1a_Fig_templates_pale.png`, `S1b_Fig_templates_medium.png` |
| Fig 10 (SQIs by skin tone) | `results/figures/Fig10_SQI_{reflectance,transmittance}.png` |
| Fig 11 (effect sizes by family and channel) | `results/figures/Fig11_effect_size_summary.png` |
| S2 Fig (morphology distributions) | `results/figures/S2{a,b}_Fig_morphology_*.png` |
| Table 3 (moderate / large effects) | `results/tables_manuscript/Table3_moderate_large_effects.csv` |
| S1 Table (feature classification) | `results/tables_manuscript/S1_Table_feature_classification.csv` |
| S2 Table (all 235 JT tests) | `results/tables_manuscript/S2_Table_JT_trend_results.csv` |
| S3 Table (per-state pale → dark SMD) | `results/tables_manuscript/S3_Table_pale_to_dark_SMD_by_state.csv` |
| Retained pulse counts | `results/tables/pulse_counts.csv` |
| MAP range of the protocol | `results/haemodynamics/state_haemodynamics.csv` |

## Repository layout

```
config.py                    paths, protocol and filter parameters
ppgphantom/
    filtering.py             AC (0.5-8 Hz), DC (<0.5 Hz) and noise (>8 Hz) components
    pyppg_features.py        pyPPG fiducials, pulse segmentation, morphology biomarkers
    sqi.py                   the nine SQIs (article Table 1)
    segmentation.py          pulse container and pooled template
    templates.py             per-state pulse templates (Fig 9, S1 Fig)
    features.py              analysed features, families and descriptions
    trend.py                 Jonckheere-Terpstra test and state consistency
    figures.py               Fig 9-11, S1 Fig, S2 Fig
scripts/                     01-03 entry points (above)
tools/
    build_zenodo_dataset.py  assembles the Zenodo deposit from the lab's working folders
results/tables/              step 1 output, read by step 2 (included)
results/figures/, tables_manuscript/, haemodynamics/   created by steps 2-3 (not included)
```

## Method notes

- The AC component is inverted in the raw recordings (front-end sign convention). `filtering.ac_component` negates it once.
- pyPPG runs on the recombined, correctly oriented AC + DC signal, using its own internal 0.5–8 Hz filter. The first detected pulse of each recording is discarded, because its onset cannot be verified.
- SQIs are computed on exactly the pulses pyPPG produced biomarkers for. No additional quality gating is applied.
- Transmittance green is in the dataset but excluded from the analysis: it contains no pulsatile signal.
- The effect size is `r = Z / sqrt(N)` with N = 27, clipped to [−1, 1].

## Licence

Code: MIT (see `LICENSE`). Data: CC BY 4.0 (see the Zenodo record).

## Citation

Please cite the article and, as appropriate:

- **Code:** https://doi.org/10.5281/zenodo.23210570 (all versions; v1.0.0 is 10.5281/zenodo.23210571)
- **Data:** https://doi.org/10.5281/zenodo.23208189

See also `CITATION.cff`.
