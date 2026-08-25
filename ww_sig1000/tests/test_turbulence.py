"""Data-free unit tests for ww_sig1000.turbulence pure functions."""
import numpy as np

from ww_sig1000.turbulence import (fit_kolmogorov, epsilon_from_A, angular_demean,
                                    _raw_spectrum, _k_grid, hr_bins, C_K,
                                    _structure_function_profile)


def test_sf_white_noise_is_all_noise_floor():
    # white noise -> D(r) = 2 sigma^2 flat -> all variance in the N intercept, A ~ 0
    rng = np.random.default_rng(0)
    cs, npings, ncols, sig = 0.06, 600, 64, 0.01
    z = np.tile(cs * np.arange(ncols), (npings, 1))
    w = rng.normal(0, sig, (npings, ncols))
    eps, N, A = _structure_function_profile(w, z, np.array([2.0]), cellsize=cs, dep_res=1.0)
    assert abs(A[0]) < 5e-5                       # negligible turbulence slope
    assert np.isclose(N[0], 2 * sig ** 2, rtol=0.3)   # intercept = 2 sigma^2


def test_sf_recovers_turbulent_slope():
    # a k^-5/3 field (Hurst 1/3) has D(r) ~ r^2/3 -> positive A, finite eps
    rng = np.random.default_rng(1)
    cs, npings, ncols, sig = 0.06, 800, 64, 0.01
    z = np.tile(cs * np.arange(ncols), (npings, 1))
    k = np.fft.rfftfreq(ncols); k[0] = k[1]
    ph = rng.uniform(0, 2 * np.pi, (npings, len(k)))
    f = np.fft.irfft((k ** (-5.0 / 6.0))[None, :] * np.exp(1j * ph), n=ncols, axis=1)
    f *= sig / f.std()
    eps, N, A = _structure_function_profile(f, z, np.array([2.0]), cellsize=cs, dep_res=1.0)
    assert A[0] > 0 and np.isfinite(eps[0])


def test_fit_kolmogorov_recovers_N_A():
    k = np.linspace(1.0, 20.0, 60)
    N, A = 1.0e-8, 1.0e-6
    S = N + A * k ** (-5.0 / 3.0)
    fN, fA = fit_kolmogorov(k, S)
    np.testing.assert_allclose(fN, N, rtol=1e-6)
    np.testing.assert_allclose(fA, A, rtol=1e-6)
    np.testing.assert_allclose(epsilon_from_A(fA), (A / C_K) ** 1.5, rtol=1e-6)


def test_epsilon_from_A_sign():
    assert np.isnan(epsilon_from_A(-1.0))
    assert epsilon_from_A(0.53) > 0
    out = epsilon_from_A(np.array([0.53, -0.1, 4 * 0.53]))
    assert np.isnan(out[1]) and out[0] > 0 and out[2] > out[0]


def test_angular_demean_removes_large_offset():
    # A constant velocity larger than the ambiguity wraps; demean -> ~0.
    v_a = 0.18
    vel = np.full((1, 30), 0.5)              # 0.5 m/s > v_a -> aliased
    out = angular_demean(vel, v_a, skip_cells=15)
    np.testing.assert_allclose(out, 0.0, atol=1e-9)


def test_angular_demean_preserves_fluctuation():
    # offset (aliased) + small fluctuation -> demean recovers fluctuation minus its
    # interior mean (the demean is over cells skip:-1).
    v_a = 0.18
    n = 30
    fluct = 0.02 * np.sin(2 * np.pi * np.arange(n) / 10.0)
    vel = (0.5 + fluct)[None, :]
    out = angular_demean(vel, v_a, skip_cells=15)[0]
    expected = fluct - fluct[15:-1].mean()
    np.testing.assert_allclose(out, expected, atol=2e-3)


def test_raw_spectrum_peaks_at_input_wavenumber():
    cellsize = 0.06
    n = 128
    k0 = 2.0                                  # cyc/m
    x = cellsize * np.arange(n)
    v = np.sin(2 * np.pi * k0 * x)
    f, ft = _raw_spectrum(v, cellsize)
    assert abs(f[np.nanargmax(ft)] - k0) < 0.2


def test_k_grid_edges_bracket_ks():
    # band edges must bracket their own ks node and be monotonic
    ks, lo, hi = _k_grid(cellsize=0.06, ks_len=50)
    assert ks.size == 25
    assert np.all(hi > lo)
    assert np.all(np.diff(ks) > 0)


def test_hr_bins_matches_nortek_grid():
    # dep_res 3 m, cellsize 0.06 -> centres 0.75, 2.25, 3.75, ... (67 bins to 100 m)
    centres, win = hr_bins(cellsize=0.06, dep_res=3.0, max_dep=100.0)
    assert win == 3.0
    np.testing.assert_allclose(centres[:3], [0.75, 2.25, 3.75], atol=1e-9)
    assert centres.size == 67
    np.testing.assert_allclose(centres[1] - centres[0], 1.5)   # dep_res/2 spacing
