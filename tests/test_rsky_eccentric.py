"""Accuracy of the eccentric-orbit sky separation (_rsky), against an extended-precision reference."""
import numpy as np
import pytest

from batman import _rsky

LD = np.longdouble


def reference(t, tc, per, a, inc, e, w):
    """Separation of the centres in units of the stellar radius, in long double precision."""
    t, tc, per, a, inc, e, w = (LD(x) if not isinstance(x, np.ndarray) else x.astype(LD)
                                for x in (t, tc, per, a, inc, e, w))
    pi = LD("3.14159265358979323846264338327950288")
    f0 = pi / 2 - w
    E0 = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(f0 / 2))
    tp = tc - per * (E0 - e * np.sin(E0)) / (2 * pi)
    x = (t - tp) / per
    M = 2 * pi * (x - np.round(x))
    sg = np.where(M < 0, -1, 1).astype(LD)
    Ma = np.abs(M)
    E = Ma.copy()
    for _ in range(80):                                   # Newton from a safe start
        E = E - (E - e * np.sin(E) - Ma) / (1 - e * np.cos(E))
        E = np.clip(E, 0, pi)
    E = sg * E
    f = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))
    s = np.sin(f + w)
    return a * (1 - e * e) / (1 + e * np.cos(f)) * np.sqrt(1 - s * s * np.sin(inc) ** 2), s * np.sin(inc)


@pytest.mark.parametrize("e", [1e-4, 0.05, 0.3, 0.6, 0.9, 0.99])
def test_separation_matches_reference(e):
    rng = np.random.default_rng(7)
    per, a, inc, w, tc = 3.7, 11.0, np.radians(87.5), np.radians(63.0), 0.21
    t = np.concatenate([rng.uniform(tc - 0.1 * per, tc + 0.1 * per, 200),
                        rng.uniform(-4 * per, 4 * per, 200), tc + per * np.arange(-3, 4)])
    d = _rsky._rsky(t, tc, per, a, inc, e, w, 1, 1)
    ref, z = reference(t, tc, per, a, inc, e, w)
    on = np.asarray(z > 1e-9)                             # primary-transit half of the orbit
    assert np.max(np.abs(d[on] - ref[on].astype(float))) < 1e-11


def test_block_independence_and_nan():
    rng = np.random.default_rng(3)
    t = rng.uniform(-20, 20, 101)
    args = (0.1, 3.3, 9.0, np.radians(87.0), 0.45, np.radians(120.0), 1)
    full = _rsky._rsky(t, *args, 1)
    for n in (1, 2, 7, 8, 9, 16, 17, 100):                # lengths around the block size
        assert np.array_equal(_rsky._rsky(t[:n], *args, 1), full[:n])
    perm = rng.permutation(t.size)
    assert np.array_equal(_rsky._rsky(t[perm], *args, 1), full[perm])
    tt = np.array([0.0, np.nan, 1.0, np.inf, 2.0])
    r = _rsky._rsky(tt, *args, 1)
    assert np.isnan(r[1]) and np.isnan(r[3]) and np.isfinite(r[[0, 2, 4]]).all()
