"""
features.py -- the analysed feature set (38 morphology features + 9 SQIs),
their families (manuscript S1 Table) and one-line descriptions.

Morphology features exported by pyppg_features.py but NOT analysed are listed
in NUMERICALLY_UNSTABLE: scale-dependent (retain a raw amplitude term) or with
>5% of recordings outside 1.5 x IQR.
"""

NUMERICALLY_UNSTABLE = {
    "RIp1", "RIp2", "SC", "AGImod", "Ta/Tpi", "Ac/Aa", "Ad/Aa", "Tp1-dp/Tpi",
    "Tdw50/Tsw50", "Tpw50/Tsp", "(Tv-Tb)/Tpi", "Tb-d/Tpi", "Tpw75/Tpi",
    "Tdw66/Tsw66", "Tdw75/Tsw75", "Tdw90/Tsw90", "Asp/Aoff", "Tpw75/Tsp", "Tsp/Asp",
    "Asp/(Tpi-Tsp)", "Asp/deltaT",
}

# feature -> (family, description). Order fixes the row order of every output table.
MORPHOLOGY_FEATURES = {
    "CT": ("Timing", "Time from pulse onset to systolic peak, as a fraction of pulse interval"),
    "APG b/a": ("Amplitude/area", "Ratio of APG turning-point amplitudes b and a"),
    "Ts/Td": ("Timing", "Ratio of systolic to diastolic phase duration"),
    "RI": ("Amplitude/area", "Reflection index: ratio of diastolic to systolic peak amplitude"),
    "AUC": ("Amplitude/area", "Total area under the pulse, normalised to amplitude x duration"),
    "FWHM": ("Width/duration", "Pulse width at half maximum amplitude"),
    "IPA": ("Amplitude/area", "Ratio of diastolic to systolic pulse area"),
    "AI": ("Composite index", "Augmentation index: relative amplitude of the late (p2) vs early (p1) systolic wave"),
    "AGI": ("Composite index", "Ageing index: combination of APG amplitudes b, c, d, e relative to a"),
    "(Ac-Ab)/Aa": ("Composite index", "Composite ratio combining APG amplitudes a, b, c"),
    "(Ad-Ab)/Aa": ("Composite index", "Composite ratio combining APG amplitudes a, b, d"),
    "(Tu-Ta)/Tpi": ("Timing", "Interval between VPG point u and APG point a"),
    "AGIinf": ("Composite index", "Modified ageing index variant"),
    "AUCdia": ("Amplitude/area", "Area under the diastolic portion of the pulse"),
    "AUCsys": ("Amplitude/area", "Area under the systolic portion of the pulse"),
    "Ae/Aa": ("Amplitude/area", "Ratio of APG turning-point amplitudes e to a"),
    "Af/Aa": ("Amplitude/area", "Ratio of APG turning-point amplitudes f to a"),
    "Ap2/Ap1": ("Amplitude/area", "Ratio of JPG turning-point amplitudes p2 to p1"),
    "Au/Asp": ("Amplitude/area", "Ratio of VPG turning-point amplitude u to systolic peak amplitude"),
    "Av/Au": ("Amplitude/area", "Ratio of VPG turning-point amplitudes v to u"),
    "Aw/Au": ("Amplitude/area", "Ratio of VPG turning-point amplitudes w to u"),
    "IPAD": ("Composite index", "Combines the systolic/diastolic area ratio (IPA) with the normalised d-point amplitude"),
    "Tb-c/Tpi": ("Timing", "Interval between APG turning points b and c"),
    "Tb/Tpi": ("Timing", "Timing of APG (2nd derivative) turning point b"),
    "Tc/Tpi": ("Timing", "Timing of APG turning point c"),
    "Td/Tpi": ("Timing", "Timing of APG turning point d"),
    "Tdw10/Tsw10": ("Width/duration", "Ratio of diastolic- to systolic-side pulse width at 10% amplitude"),
    "Tdw25/Tsw25": ("Width/duration", "Ratio of diastolic- to systolic-side pulse width at 25% amplitude"),
    "Tdw33/Tsw33": ("Width/duration", "Ratio of diastolic- to systolic-side pulse width at 33% amplitude"),
    "Te/Tpi": ("Timing", "Timing of APG turning point e"),
    "Tf/Tpi": ("Timing", "Timing of APG turning point f"),
    "Tp2-dp/Tpi": ("Timing", "Interval between JPG (3rd derivative) point p2 and the diastolic peak"),
    "Tpw25/Tpi": ("Width/duration", "Pulse width at 25% peak amplitude, relative to pulse interval"),
    "Tpw25/Tsp": ("Width/duration", "Pulse width at 25% peak amplitude, relative to time-to-systolic-peak"),
    "Tsp/Tpi": ("Timing", "Time to systolic peak, as a fraction of pulse interval"),
    "Tu/Tpi": ("Timing", "Timing of VPG (1st derivative) turning point u"),
    "Tv/Tpi": ("Timing", "Timing of VPG turning point v"),
    "Tw/Tpi": ("Timing", "Timing of VPG turning point w"),
}

SQI_FEATURES = {
    "TMCC": "Template matching correlation coefficient: waveform-shape similarity to the pulse template",
    "AC/DC": "Ratio of pulsatile (AC) to non-pulsatile (DC) signal amplitude",
    "SNR": "Signal-to-noise ratio of the pulsatile signal",
    "Skewness": "Asymmetry of the pulse waveform amplitude distribution",
    "PSD ratio": "Power spectral density ratio: concentration of spectral energy in the pulsatile band",
    "Shannon Entropy": "Waveform complexity / information content of the pulsatile signal",
    "Kurtosis": "Peakedness of the pulse waveform amplitude distribution",
    "Spec. Flatness": "Spectral flatness (Wiener entropy): how noise-like versus tonal the pulsatile spectral content is",
    "Sample Entropy": "Waveform regularity / complexity of the pulsatile signal",
}

MORPHOLOGY_FAMILIES = ["Timing", "Width/duration", "Amplitude/area", "Composite index"]


def feature_table():
    """One row per analysed feature: feature, feature_family, source_table, description."""
    import pandas as pd
    rows = [{"feature": f, "feature_family": fam, "source_table": "morphology", "description": desc}
            for f, (fam, desc) in MORPHOLOGY_FEATURES.items()]
    rows += [{"feature": f, "feature_family": "SQI", "source_table": "sqi", "description": desc}
             for f, desc in SQI_FEATURES.items()]
    return pd.DataFrame(rows)
