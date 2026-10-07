"""
templates.py -- pooled pulse templates per skin tone x heart rate x target
flow x channel (replicates pooled), stored as a small long-format table so the
template figures can be drawn without re-running pyPPG.
"""

import numpy as np
import pandas as pd

from ppgphantom.segmentation import pooled_template, resample_pulse


class TemplateCollector:
    """Accumulates resampled normalised pulses per (skin, hr, co_lmin, channel)."""

    def __init__(self):
        self._pulses = {}

    def add(self, skin: str, hr: int, flow: int, channel: str, ac_norms) -> None:
        stack = self._pulses.setdefault((skin, hr, flow, channel), [])
        stack.extend(resample_pulse(a) for a in ac_norms)

    def to_table(self) -> pd.DataFrame:
        frames = []
        for (skin, hr, flow, channel), pulses in self._pulses.items():
            med, p05, p95 = pooled_template(np.vstack(pulses))
            frames.append(pd.DataFrame({
                "skin": skin, "hr": hr, "co_lmin": flow, "channel": channel,
                "n_pulses": len(pulses), "point": np.arange(len(med)),
                "median": med, "p05": p05, "p95": p95,
            }))
        return pd.concat(frames, ignore_index=True)
