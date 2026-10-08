"""Tests of the result classes and their methods."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from focus_cpt import (
    Detector,
    DetectorStatistics,
    DetectorSummary,
    OfflineResult,
    OfflineSummary,
    ProjectionIndexes,
    focus_offline,
    generate_projection_indexes,
)


def make_data(seed=1):
    rng = np.random.default_rng(seed)
    return np.concatenate([rng.normal(size=100), rng.normal(1, size=50)])


def make_detector(type="univariate", Y=None, **kwargs):
    det = Detector(type=type, **kwargs)
    for y in make_data() if Y is None else Y:
        det.update(y)
    return det


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def test_statistics():
    det = make_detector()
    r = det.get_statistics(family="gaussian")
    assert isinstance(r, DetectorStatistics)
    assert r.family == "gaussian"
    assert r.stat == r["stat"]
    assert repr(r).startswith("focus statistics (family: gaussian)")


def test_detector_summary():
    det = make_detector()
    s = det.summary()
    assert isinstance(s, DetectorSummary)
    assert s["n"] == 150
    assert s["n_candidates"] == det.get_n_candidates()
    np.testing.assert_array_equal(s["candidates"], np.unique(det.get_candidates()["tau"]))
    assert s["statistics"] is None
    assert repr(s).startswith("focus detector: summary")

    s2 = det.summary(family="gaussian")
    assert s2["statistics"] == det.get_statistics(family="gaussian")
    assert "focus statistics" in repr(s2)

    det_arp = make_detector("arp", rho=[0.5])
    s_arp = det_arp.summary(family="arp")
    assert s_arp["ar_order"] == 1
    assert s_arp["n_candidates"] == det_arp.get_n_candidates() > 0
    np.testing.assert_array_equal(s_arp["candidates"], np.unique(det_arp.get_candidates()["tau"]))
    assert "AR order" in repr(s_arp)

    s_side = make_detector("univariate_one_sided", side="left").summary()
    assert 'side = "left"' in repr(s_side)


def test_offline_result():
    Y = make_data(123)
    res = focus_offline(Y, threshold=20, type="univariate", family="gaussian")
    assert isinstance(res, OfflineResult)
    assert repr(res).startswith("focus offline detection")
    s = res.summary()
    assert isinstance(s, OfflineSummary)
    assert s["statistics"][0]["max"] == pytest.approx(np.max(res.stat))

    res.plot()
    assert len(res.plot(data=Y)) == 2

    res_np = focus_offline(Y, threshold=[10, 5], type="npfocus", family="npfocus",
                           quantiles=[-0.67, 0, 0.67])
    assert len(res_np.plot(data=Y)) == 2

    # Multivariate data: one panel per dimension
    rng = np.random.default_rng(1)
    Y_mv = np.vstack([rng.normal(size=(100, 3)), rng.normal(1, size=(50, 3))])
    res_mv = focus_offline(Y_mv, threshold=30, type="multivariate", family="gaussian")
    assert len(res_mv.plot(data=Y_mv)) == 4


def test_projection_indexes():
    proj = generate_projection_indexes(6, 2)
    assert isinstance(proj, ProjectionIndexes)
    assert (proj.d, proj.p) == (6, 2)
    assert len(proj) == 6
    np.testing.assert_array_equal(proj[5], [5, 0])
    assert "6 of size 2" in repr(proj)

    first = proj[:2]
    assert isinstance(first, ProjectionIndexes)
    assert len(first) == 2
    assert proj.to_array().shape == (6, 2)
    assert first.to_array().shape == (2, 2)

    rng = np.random.default_rng(4)
    Y = rng.normal(size=(100, 6))
    det_cls = make_detector("multivariate", Y=Y, dim_indexes=proj)
    det_lst = make_detector("multivariate", Y=Y, dim_indexes=[list(idx) for idx in proj])
    assert det_cls.get_statistics("gaussian") == det_lst.get_statistics("gaussian")
