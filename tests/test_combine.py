"""The impact figures: the ABCD terms apart from the other nuisances."""

from vbsvvh_stat.combine import split_impacts


def _params(*names):
    return [{"name": n} for n in names]


def test_abcd_terms_get_their_own_figure():
    params = _params("CMS_HIG26013_MET1FJ_Run2_abcd_B", "CMS_HIG26013_MET1FJ_Run2_abcd_closure", "ps_fsr", "lumi_13TeV_1516_l")
    figs = {k: [p["name"] for p in v] for k, v in split_impacts(params).items()}
    assert figs[""] == [p["name"] for p in params]
    assert figs["_abcd"] == ["CMS_HIG26013_MET1FJ_Run2_abcd_B", "CMS_HIG26013_MET1FJ_Run2_abcd_closure"]
    assert figs["_noabcd"] == ["ps_fsr", "lumi_13TeV_1516_l"]


def test_without_abcd_terms_one_figure():
    assert list(split_impacts(_params("ps_fsr", "ps_isr"))) == [""]
