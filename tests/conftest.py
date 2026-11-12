"""Synthetic channels in the input-contract format."""

import numpy as np
import pytest
import yaml


def th1(values, variances=None, edges=None, sumw2=True):
    import hist
    values = np.atleast_1d(np.asarray(values, dtype=float))
    edges = np.arange(len(values) + 1, dtype=float) if edges is None else np.asarray(edges, float)
    if not sumw2:
        return values, edges
    h = hist.Hist(hist.axis.Variable(edges), storage=hist.storage.Weight())
    h.view()["value"] = values
    h.view()["variance"] = values * 0.01 if variances is None else np.atleast_1d(variances)
    return h


def write_channel(directory, tag, hists, abcd=None, raw_keys=(), theory="acceptance"):
    """``hists``: ``{"era/region/name": values | (values, variances) | th1}`` -> channel.yaml path."""
    import uproot
    directory.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(directory / "shapes.root") as f:
        for key, content in hists.items():
            if isinstance(content, tuple):
                content = th1(*content)
            elif not hasattr(content, "view") and not isinstance(content, tuple):
                content = th1(content, np.atleast_1d(content) if key.endswith("data_obs") else None)
            f[key] = content
        for key, content in raw_keys:
            f[key] = content
    cfg = {"channel": tag}
    if theory:
        cfg["theory"] = theory
    if abcd:
        cfg["abcd"] = abcd
    path = directory / "channel.yaml"
    path.write_text(yaml.safe_dump(cfg))
    return path


def abcd_hists(eras, sig=1.0, bkg=(2.0, 20.0, 30.0, 300.0), data=(0.0, 21.0, 29.0, 310.0), nuis=()):
    """A counting ABCD channel: signal + data-driven ``bkg`` + data in regions A, B, C, D."""
    out = {}
    for era in eras:
        for i, r in enumerate("ABCD"):
            s = sig * (1.0 if r == "A" else 0.2)
            out[f"{era}/{r}/Signal_C2V_1p5_C3_1p0"] = s
            out[f"{era}/{r}/bkg"] = bkg[i]
            out[f"{era}/{r}/data_obs"] = data[i]
            for n, k in nuis:
                n = n.format(era=era)
                out[f"{era}/{r}/Signal_C2V_1p5_C3_1p0__{n}Up"] = s * k
                out[f"{era}/{r}/Signal_C2V_1p5_C3_1p0__{n}Down"] = s / k
    return out


ABCD = {"background": "bkg", "A": "A", "B": "B", "C": "C", "D": "D"}


@pytest.fixture
def analysis():
    from vbsvvh_stat.config import Analysis
    return Analysis()
