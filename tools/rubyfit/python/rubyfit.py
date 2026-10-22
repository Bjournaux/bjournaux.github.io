#!/usr/bin/env python3
"""Fast batch fitting of ruby fluorescence (R1/R2) spectra and pressure calculation.

Same physical model as the "Ruby Pressure Calibration GUI" notebook (two Pearson VII
peaks on a cubic background, R1 shift of a sample ruby relative to a reference ruby
at the same temperature), reworked for speed and robustness:

  * peaks are parametrised by height / center / HWHM / shape (no gamma functions);
  * R2 is tied to R1 through a bounded R1-R2 splitting, so the peaks cannot swap;
  * the background polynomial uses a coordinate scaled to [-1, 1] (well conditioned);
  * analytic Jacobian + scipy.optimize.least_squares, initial guesses from the data;
  * the fit is restricted to a window around the doublet (configurable);
  * files can be processed in parallel.

File names are parsed for the experimental conditions, e.g.
    ice_VI_H2O_inside_297K_ruby3 044.csv
      side        -> "inside"  = sample (high pressure), "outside" = reference
      temperature -> 297 K     (sample and reference are matched by temperature)
      ruby        -> 3         (repeat spectra of the same ruby are averaged)

Command line:
    python rubyfit.py rubyfiles/                       # results in ./rubyfit_results
    python rubyfit.py rubyfiles/ -o out --plots all --jobs 8
    python rubyfit.py "rubyfiles/*297K*.csv" --ref-wavelength 694.22   # single-ruby mode
    python rubyfit.py rubyfiles/ --profile split_pearson7              # asymmetric peaks

Python / Jupyter:
    import rubyfit
    fits = rubyfit.fit_files(rubyfit.find_files(["rubyfiles"]))
    fits, per_ruby, summary = rubyfit.pressure_tables(fits)
"""
from __future__ import annotations

import argparse
import glob
import math
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from functools import partial

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.signal import medfilt

LN2 = math.log(2.0)

DEFAULT_WINDOW = (-4.5, 4.0)      # fit window relative to the R1 peak (nm)
DEFAULT_MIN_R2 = 0.98             # r-squared (inside the fit window) below which a fit is flagged
DEFAULT_REF_KEYS = ("outside", "backside", "reference", "ref")
DEFAULT_SAMPLE_KEYS = ("inside", "chamber", "sample", "hp")

PROFILES = {
    "pearson7": "Pearson VII",
    "split_pearson7": "Split Pearson VII (asymmetric)",
}
DEFAULT_PROFILE = "pearson7"
# half widths of each peak: one for a symmetric peak, left/right (below/above the maximum) for a split one
PEAK_WIDTHS = {"pearson7": ("hwhm",), "split_pearson7": ("hwhm_left", "hwhm_right")}


def param_names(profile=DEFAULT_PROFILE):
    """Names of the fit parameters, in order, for a line-shape profile."""
    def peak(n):
        return tuple(f"{n}_{w}" for w in PEAK_WIDTHS[profile]) + (f"{n}_shape",)
    return (("r1_height", "r1_center") + peak("r1") + ("r2_height", "r1_r2_split") + peak("r2")
            + ("bg0", "bg1", "bg2", "bg3"))


PARAM_NAMES = param_names()

# Flags that make a fit unusable (excluded from the pressure averages).
# Shape parameters at their bounds or R2 > R1 are reported but not disqualifying.
SERIOUS_FLAGS = {"fit_not_converged", "low_r2", "r1_center_at_bound", "r1_r2_split_at_bound",
                 "r1_hwhm_at_bound", "r2_hwhm_at_bound", "no_uncertainty"}

# Columns produced by RubyFit.row() (filled with NaN when a file cannot be fitted)
FIT_COLUMNS = ("r1_nm", "r1_err_nm", "r2_nm", "r2_err_nm", "r1_fwhm_nm", "r2_fwhm_nm",
               "r1_shape", "r2_shape", "r1_height", "r2_height", "r1_r2_split_nm",
               "r_squared", "redchi", "nfev", "r1_asymmetry", "profile")


# ----------------------------------------------------------------------------------------
# Pressure scales
# ----------------------------------------------------------------------------------------

SCALES = {
    "mao1986": "Mao et al. (1986), quasi-hydrostatic: P = A/B [(l/l0)^B - 1], A = 1904 GPa, B = 7.665",
    "shen2020": "Shen et al. (2020), ruby2020: P = A (dl/l0) [1 + B (dl/l0)], A = 1870 GPa, B = 5.63",
}


def pressure_gpa(lam, lam0, scale="mao1986"):
    """Pressure (GPa) from the R1 wavelength `lam` and the reference R1 wavelength `lam0` (nm).

    Returns (P, dP/dlam, dP/dlam0) so that uncertainties can be propagated.
    """
    lam = np.asarray(lam, dtype=float)
    lam0 = np.asarray(lam0, dtype=float)
    ratio = lam / lam0
    if scale == "mao1986":
        A, B = 1904.0, 7.665
        P = A / B * (ratio ** B - 1.0)
        dP_dlam = A * ratio ** (B - 1.0) / lam0
        dP_dlam0 = -A * ratio ** B / lam0
    elif scale == "shen2020":
        A, B = 1870.0, 5.63
        r = ratio - 1.0
        P = A * r * (1.0 + B * r)
        dP_dr = A * (1.0 + 2.0 * B * r)
        dP_dlam = dP_dr / lam0
        dP_dlam0 = -dP_dr * lam / lam0 ** 2
    else:
        raise ValueError(f"unknown pressure scale {scale!r}; choose from {sorted(SCALES)}")
    return P, dP_dlam, dP_dlam0


# ----------------------------------------------------------------------------------------
# Input
# ----------------------------------------------------------------------------------------

def find_files(inputs, pattern="*.csv"):
    """Expand folders (using `pattern`), glob patterns and plain paths into a sorted file list."""
    files = []
    for item in inputs:
        if os.path.isdir(item):
            files.extend(glob.glob(os.path.join(item, pattern)))
        elif any(ch in item for ch in "*?["):
            files.extend(glob.glob(item))
        else:
            files.append(item)
    return sorted(dict.fromkeys(os.path.abspath(f) for f in files))


_NUMBER_RE = re.compile(r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?")
_DATA_LINE_RE = re.compile(r"[+-]?\.?\d")
_COMMA_NUMBER_RE = re.compile(r"[+-]?\d+(?:,\d+)?(?:[eE][+-]?\d+)?")
_X_NAME_RE = re.compile(r"wave|lambda|nm\b", re.IGNORECASE)
_Y_NAME_RE = re.compile(r"intens|count|signal|cts", re.IGNORECASE)


def _decimal_comma_columns(line):
    """True for a space-separated line whose numbers use a decimal comma ("694,12 1234")."""
    f = line.split()
    return len(f) >= 2 and "," in f[0] + f[1] and all(_COMMA_NUMBER_RE.fullmatch(v) for v in f[:2])


def parse_spectrum_text(text):
    """Wavelength and intensity columns from the text of a CSV or TXT spectrum, sorted by wavelength.

    Columns may be separated by commas, semicolons, tabs or spaces; with semicolons, tabs or
    spaces a decimal comma is accepted. Header, comment and metadata lines are skipped. If a header
    names the columns ("Wavelength", "Intensity", "Counts", ...) those are used, otherwise the
    first two columns. (The notebook's pd.read_csv treated the first data row as a header.)
    """
    lines = [line.strip() for line in text.lstrip("﻿").splitlines()]
    first = next((i for i, line in enumerate(lines) if _DATA_LINE_RE.match(line)), None)
    if first is None:
        raise ValueError("no numeric data found: expected two columns (wavelength in nm, intensity)")
    data = [line for line in lines[first:] if _DATA_LINE_RE.match(line)]

    def share(test):
        return sum(1 for line in data if test(line)) / len(data)

    if share(lambda s: ";" in s) >= 0.5:
        sep, decimal_comma = r"\s*;\s*", True
    elif share(lambda s: "\t" in s) >= 0.5:
        sep, decimal_comma = r"\s*\t\s*", True
    elif share(lambda s: "," in s) >= 0.5 and share(_decimal_comma_columns) < 0.5:
        sep, decimal_comma = r"\s*,\s*", False
    else:
        sep, decimal_comma = r"\s+", True

    def fields(line):
        out = [f.strip().strip("\"'") for f in re.split(sep, line)]
        return [f.replace(",", ".") for f in out] if decimal_comma else out

    ix, iy = 0, 1                                        # a header naming the columns overrides this
    header = next((lines[i] for i in range(first - 1, -1, -1) if lines[i]), None)
    if header:
        names = [f.strip().strip("\"'") for f in re.split(sep, header)]
        if len(names) == len(fields(data[0])):
            hx = next((i for i, n in enumerate(names) if _X_NAME_RE.search(n)), None)
            hy = next((i for i, n in enumerate(names) if _Y_NAME_RE.search(n) and i != hx), None)
            if hx is not None and hy is not None:
                ix, iy = hx, hy
    xs, ys = [], []
    for line in data:
        f = fields(line)
        if len(f) > max(ix, iy) and _NUMBER_RE.fullmatch(f[ix]) and _NUMBER_RE.fullmatch(f[iy]):
            xs.append(float(f[ix]))
            ys.append(float(f[iy]))
    x, y = np.array(xs, dtype=float), np.array(ys, dtype=float)
    order = np.argsort(x, kind="stable")
    return x[order], y[order]


def read_spectrum(path):
    """Read a CSV or TXT spectrum file: wavelength (nm) and intensity, see parse_spectrum_text."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        return parse_spectrum_text(fh.read())


# A temperature is a number directly followed by a capital K ("297K", "297.5K", "297 K", "297_K"),
# not glued to other letters or digits (a leading "T" is allowed: "T297K"); "297,5K" means 297.5 K.
_TEMP_RE = re.compile(r"(?:^|[^A-Za-z0-9.,]|[Tt])(\d+(?:[.,]\d+)?)[\s_]?K(?![A-Za-z])")
_RUBY_RE = re.compile(r"ruby[\s_-]?(\d+)", re.IGNORECASE)


def parse_name(path, ref_keys=DEFAULT_REF_KEYS, sample_keys=DEFAULT_SAMPLE_KEYS):
    """Side (reference/sample), temperature and ruby number from a file name."""
    stem = os.path.splitext(os.path.basename(path))[0]
    words = re.sub(r"([a-z])([A-Z])", r"\1_\2", stem)    # "insideDAC" -> "inside_DAC"
    tokens = set(re.split(r"[^a-z0-9]+", words.lower()))
    if tokens & {k.lower() for k in ref_keys}:
        side = "reference"
    elif tokens & {k.lower() for k in sample_keys}:
        side = "sample"
    else:
        side = ""
    m = _TEMP_RE.search(stem)
    temperature = float(m.group(1).replace(",", ".")) if m else math.nan
    m = _RUBY_RE.search(stem)
    ruby = str(int(m.group(1))) if m else ""             # "ruby01" and "ruby1" are the same ruby
    return {"file": os.path.basename(path), "side": side, "temperature_K": temperature, "ruby": ruby}


# ----------------------------------------------------------------------------------------
# Model
# ----------------------------------------------------------------------------------------

def pearson7(x, height, center, hwhm, shape):
    """Pearson VII peak: height * [1 + (2^(1/shape) - 1) ((x - center)/hwhm)^2]^(-shape).

    shape = 1 is a Lorentzian, shape -> infinity tends to a Gaussian.
    """
    k = 2.0 ** (1.0 / shape) - 1.0
    s = (x - center) / hwhm
    return height * (1.0 + k * s * s) ** (-shape)


def peak_profile(x, height, center, widths, shape):
    """Pearson VII with one half width, or split Pearson VII with (left, right) half widths."""
    x = np.asarray(x, dtype=float)
    hwhm = widths[0] if len(widths) == 1 else np.where(x < center, widths[0], widths[1])
    return pearson7(x, height, center, hwhm, shape)


def _peak_derivs(x, height, center, widths, shape):
    """Peak value and derivatives w.r.t. height, center, each half width and shape."""
    left = x < center
    hwhm = widths[0] if len(widths) == 1 else np.where(left, widths[0], widths[1])
    k1 = 2.0 ** (1.0 / shape)
    k = k1 - 1.0
    s = (x - center) / hwhm
    s2 = s * s
    u = 1.0 + k * s2
    g = u ** (-shape)
    f = height * g
    q = 2.0 * shape * k * f / (u * hwhm)
    d_width = q * s2
    d_widths = [d_width] if len(widths) == 1 else [np.where(left, d_width, 0.0), np.where(left, 0.0, d_width)]
    d_shape = f * (k1 * LN2 * s2 / (shape * u) - np.log(u))
    return f, g, q * s, d_widths, d_shape


def _background(t, b0, b1, b2, b3):
    return b0 + t * (b1 + t * (b2 + t * b3))


class _Problem:
    """Residuals and Jacobian of the R1 + R2 + cubic background model (shared evaluation)."""

    def __init__(self, x, y, t, n_widths=1):
        self.x, self.y, self.t, self.nw = x, y, t, n_widths
        self.npk = 3 + n_widths                          # height, center (or splitting), widths, shape
        self._p = None
        self._jac = np.empty((x.size, 2 * self.npk + 4))
        bg = 2 * self.npk
        self._jac[:, bg] = 1.0
        self._jac[:, bg + 1] = t
        self._jac[:, bg + 2] = t * t
        self._jac[:, bg + 3] = t ** 3

    def residuals(self, p):
        nw, npk = self.nw, self.npk
        h1, c1, w1, m1 = p[0], p[1], p[2:2 + nw], p[2 + nw]
        h2, d, w2, m2 = p[npk], p[npk + 1], p[npk + 2:npk + 2 + nw], p[npk + 2 + nw]
        f1, g1, dc1, dw1, dm1 = _peak_derivs(self.x, h1, c1, w1, m1)
        f2, g2, dc2, dw2, dm2 = _peak_derivs(self.x, h2, c1 - d, w2, m2)
        J = self._jac
        J[:, 0], J[:, 1], J[:, 2 + nw] = g1, dc1 + dc2, dm1
        J[:, npk], J[:, npk + 1], J[:, npk + 2 + nw] = g2, -dc2, dm2
        for i in range(nw):
            J[:, 2 + i], J[:, npk + 2 + i] = dw1[i], dw2[i]
        self._p = p.copy()
        return f1 + f2 + _background(self.t, *p[2 * npk:]) - self.y

    def jac(self, p):
        if self._p is None or not np.array_equal(p, self._p):
            self.residuals(p)
        return self._jac.copy()


@dataclass
class RubyFit:
    """Result of fitting one spectrum."""
    x: np.ndarray            # wavelengths used in the fit (nm)
    y: np.ndarray            # intensities used in the fit
    params: dict             # best-fit values, see PARAM_NAMES (+ derived r2_center)
    errors: dict             # 1-sigma standard errors (NaN when not available)
    r_squared: float
    redchi: float
    success: bool
    nfev: int
    flags: list
    t_center: float          # background coordinate t = (x - t_center) / t_halfwidth
    t_halfwidth: float
    profile: str = DEFAULT_PROFILE

    @property
    def ok(self):
        return not (set(self.flags) & SERIOUS_FLAGS)

    def widths(self, peak):
        """Half widths of 'r1' or 'r2': (hwhm,) or (left, right) for a split profile."""
        return [self.params[f"{peak}_{w}"] for w in PEAK_WIDTHS[self.profile]]

    def components(self, x=None):
        """R1, R2 and background curves evaluated at `x` (defaults to the fitted points)."""
        x = self.x if x is None else np.asarray(x, dtype=float)
        p = self.params
        t = (x - self.t_center) / self.t_halfwidth
        return {
            "R1": peak_profile(x, p["r1_height"], p["r1_center"], self.widths("r1"), p["r1_shape"]),
            "R2": peak_profile(x, p["r2_height"], p["r2_center"], self.widths("r2"), p["r2_shape"]),
            "background": _background(t, p["bg0"], p["bg1"], p["bg2"], p["bg3"]),
        }

    def model(self, x=None):
        return sum(self.components(x).values())

    def row(self):
        p, e = self.params, self.errors
        w1, w2 = self.widths("r1"), self.widths("r2")
        return {
            "r1_nm": p["r1_center"], "r1_err_nm": e["r1_center"],
            "r2_nm": p["r2_center"], "r2_err_nm": e["r2_center"],
            "r1_fwhm_nm": w1[0] + w1[-1], "r2_fwhm_nm": w2[0] + w2[-1],
            "r1_shape": p["r1_shape"], "r2_shape": p["r2_shape"],
            "r1_height": p["r1_height"], "r2_height": p["r2_height"],
            "r1_r2_split_nm": p["r1_r2_split"],
            "r_squared": self.r_squared, "redchi": self.redchi, "nfev": self.nfev,
            "r1_asymmetry": w1[-1] / w1[0], "profile": self.profile,      # right / left half width
            "ok": self.ok, "flags": ";".join(self.flags),
        }


def _initial_guess(x, ys, c1):
    """Starting values from the (median-filtered) data; c1 is the R1 position."""
    bg = float(np.percentile(ys, 5))
    i1 = int(np.argmin(np.abs(x - c1)))
    h1 = max(ys[i1] - bg, 1.0)
    # HWHM from the high-wavelength flank of R1 (the low side overlaps with R2)
    below = np.flatnonzero((x > c1) & (ys < bg + h1 / 2))
    w1 = float(np.clip(x[below[0]] - c1, 0.05, 2.0)) if below.size else 0.3
    # R2: highest point 0.9-2.2 nm below R1
    sel = np.flatnonzero((x > c1 - 2.2) & (x < c1 - 0.9))
    if sel.size:
        j = sel[np.argmax(ys[sel])]
        split, h2 = c1 - x[j], ys[j] - bg
    else:
        split, h2 = 1.42, 0.6 * h1
    h2 = max(h2, 0.05 * h1)
    return np.array([h1, c1, w1, 1.5, h2, split, w1, 1.5, bg, 0.0, 0.0, 0.0])


def _standard_errors(jac, redchi):
    """Parameter standard errors from the Jacobian (same approach as scipy's curve_fit)."""
    _, s, vt = np.linalg.svd(jac, full_matrices=False)
    keep = s > np.finfo(float).eps * max(jac.shape) * s[0]
    if not keep.all():
        return None
    vt = vt / s[:, None]
    return vt.T @ vt * redchi


def fit_spectrum(x, y, window=DEFAULT_WINDOW, search=None, min_r2=DEFAULT_MIN_R2, profile=DEFAULT_PROFILE):
    """Fit R1 + R2 + cubic background to one spectrum.

    window  : (lo, hi) in nm relative to the R1 maximum, or None to fit the full range.
    search  : optional (lo, hi) absolute wavelength range (nm) in which to look for R1.
    profile : peak line shape, see PROFILES ("pearson7" or "split_pearson7"). For a split
              profile the center is the peak maximum.
    """
    if profile not in PROFILES:
        raise ValueError(f"unknown profile {profile!r}; choose from {sorted(PROFILES)}")
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 20:
        raise ValueError(f"only {x.size} data points")
    order = np.argsort(x, kind="stable")
    x, y = x[order], y[order]
    step = float(np.median(np.diff(x)))
    if step >= 0.9:
        raise ValueError(f"wavelength step of {step:.3g}: the first column should be the wavelength in nm "
                         "(this looks like pixel numbers)")
    ys = medfilt(y, 5)                                   # rejects single-pixel spikes

    cand = np.ones(x.size, bool) if search is None else (x >= search[0]) & (x <= search[1])
    if not cand.any():
        raise ValueError(f"no data in the R1 search range {search}")
    idx = np.flatnonzero(cand)
    c1 = float(x[idx[np.argmax(ys[idx])]])
    if not 600 <= c1 <= 900:
        raise ValueError(_units_hint(c1))

    if window is not None:
        sel = (x >= c1 + window[0]) & (x <= c1 + window[1])
        x, y, ys = x[sel], y[sel], ys[sel]
    if x.size < 20:
        raise ValueError(f"only {x.size} points in the fit window")

    t_center, t_halfwidth = 0.5 * (x[0] + x[-1]), 0.5 * (x[-1] - x[0])
    nw = len(PEAK_WIDTHS[profile])
    prob = _Problem(x, y, (x - t_center) / t_halfwidth, nw)
    inf = np.inf
    lb = np.array([0, c1 - 1.5, *[0.01] * nw, 0.6, 0, 0.8, *[0.01] * nw, 0.6, -inf, -inf, -inf, -inf])
    ub = np.array([inf, c1 + 1.5, *[3.0] * nw, 30, inf, 2.5, *[3.0] * nw, 30, inf, inf, inf, inf])
    margin = np.where(np.isfinite(ub - lb), 1e-6 * (ub - lb), 0.0)
    g = _initial_guess(x, ys, c1)                        # symmetric layout; widths repeated for a split peak
    p0 = np.array([g[0], g[1], *[g[2]] * nw, g[3], g[4], g[5], *[g[6]] * nw, g[7], *g[8:]])
    i_shape1, i_shape2 = 2 + nw, 5 + 2 * nw
    ss_tot = float(np.sum((y - y.mean()) ** 2))

    best = None
    for shape0 in (1.5, 1.0, 3.0):                       # retry other line shapes if needed
        p0[i_shape1] = p0[i_shape2] = shape0
        res = least_squares(prob.residuals, np.clip(p0, lb + margin, ub - margin), jac=prob.jac,
                            bounds=(lb, ub), method="trf", x_scale="jac",
                            ftol=1e-10, xtol=1e-10, gtol=1e-10, max_nfev=400)
        fit = _build_fit(res, x, y, lb, ub, ss_tot, min_r2, t_center, t_halfwidth, profile)
        if best is None or fit.redchi < best.redchi:
            best = fit
        if best.ok:
            break
    return best


def _units_hint(c1):
    """Error message when the R1 maximum is not where a wavelength in nm would put it."""
    if 6000 <= c1 <= 9000:
        return f"R1 found at {c1:.6g}: the wavelengths look like ångström; the tool expects nm"
    if 10000 <= c1 <= 17000:
        return f"R1 found at {c1:.6g}: the axis looks like wavenumbers (cm-1); the tool expects wavelengths in nm"
    if 0 < c1 < 6000:
        return (f"R1 found at {c1:.6g}: the axis looks like a Raman shift (cm-1); "
                "the tool expects wavelengths in nm")
    return f"R1 found at {c1:.6g}, outside 600-900 nm: the first column should be the wavelength in nm"


def _build_fit(res, x, y, lb, ub, ss_tot, min_r2, t_center, t_halfwidth, profile=DEFAULT_PROFILE):
    names = param_names(profile)
    p = res.x
    n, k = x.size, p.size
    ssr = 2.0 * res.cost
    redchi = ssr / max(n - k, 1)
    cov = _standard_errors(res.jac, redchi)
    params = dict(zip(names, p))
    params["r2_center"] = params["r1_center"] - params["r1_r2_split"]
    i_split = names.index("r1_r2_split")
    if cov is None:
        errors = {name: math.nan for name in params}
    else:
        sd = np.sqrt(np.clip(np.diag(cov), 0, None))
        errors = dict(zip(names, sd))
        errors["r2_center"] = math.sqrt(max(cov[1, 1] + cov[i_split, i_split] - 2 * cov[1, i_split], 0.0))
    r_squared = 1.0 - ssr / ss_tot if ss_tot > 0 else math.nan

    flags = []
    if not res.success:
        flags.append("fit_not_converged")
    if not r_squared >= min_r2:
        flags.append("low_r2")
    if cov is None:
        flags.append("no_uncertainty")
    tol = 1e-4 * np.where(np.isfinite(ub - lb), ub - lb, 1.0)
    for i, name in enumerate(names):                     # positions, widths, shapes, splitting
        if name.endswith("_height") or name.startswith("bg"):
            continue
        if p[i] - lb[i] < tol[i] or ub[i] - p[i] < tol[i]:
            flag = f"{name[:2]}_hwhm_at_bound" if "_hwhm" in name else f"{name}_at_bound"
            if flag not in flags:
                flags.append(flag)
    if params["r2_height"] > params["r1_height"]:
        flags.append("r2_stronger_than_r1")
    return RubyFit(x, y, params, errors, r_squared, redchi, bool(res.success), int(res.nfev),
                   flags, t_center, t_halfwidth, profile)


# ----------------------------------------------------------------------------------------
# Batch processing
# ----------------------------------------------------------------------------------------

def fit_file(path, plot_dir=None, plot_mode="flagged", ref_keys=DEFAULT_REF_KEYS,
             sample_keys=DEFAULT_SAMPLE_KEYS, **fit_kw):
    """Read, fit and (optionally) plot one file; returns one row of the results table."""
    row = parse_name(path, ref_keys, sample_keys)
    fit = None
    t0 = time.perf_counter()
    try:
        x, y = read_spectrum(path)
        t0 = time.perf_counter()
        fit = fit_spectrum(x, y, **fit_kw)
        row.update(fit.row())
    except Exception as exc:                              # keep going with the rest of the batch
        row.update(dict.fromkeys(FIT_COLUMNS, math.nan), profile=fit_kw.get("profile", DEFAULT_PROFILE),
                   ok=False, flags=f"error: {exc}")
    row["fit_ms"] = 1000 * (time.perf_counter() - t0)
    if plot_dir and fit is not None and (plot_mode == "all" or (plot_mode == "flagged" and not fit.ok)):
        os.makedirs(plot_dir, exist_ok=True)
        stem = os.path.splitext(os.path.basename(path))[0]
        plot_fit(fit, os.path.join(plot_dir, f"{stem}_fit.png"), title=os.path.basename(path))
    row["path"] = path
    return row


def fit_files(paths, jobs=1, **kw):
    """Fit many files (in parallel with jobs > 1). Keyword arguments go to fit_file()."""
    work = partial(fit_file, **kw)
    if jobs > 1 and len(paths) > 1:
        with ProcessPoolExecutor(max_workers=jobs) as pool:
            rows = list(pool.map(work, paths, chunksize=max(1, len(paths) // (4 * jobs))))
    else:
        rows = [work(p) for p in paths]
    return pd.DataFrame(rows)


def _mean_and_error(values, fit_errors):
    """Mean of independent measurements and its standard error.

    With a single measurement, fall back to its own (fit) uncertainty.
    """
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return math.nan, math.nan
    if values.size == 1:
        return float(values[0]), float(np.asarray(fit_errors, dtype=float)[0])
    return float(values.mean()), float(values.std(ddof=1) / math.sqrt(values.size))


def pressure_tables(fits, scale="mao1986", ref_wavelength=None, ref_err=0.02):
    """Pressures from fitted R1 positions.

    Dual-ruby mode (default): the reference wavelength at each temperature is the mean R1
    of the reference rubies at that temperature. Single-ruby mode: pass `ref_wavelength`
    (and its uncertainty `ref_err`), as the notebook's "Single Ruby, Ambient" mode with
    694.22 +/- 0.02 nm.

    Repeat spectra of the same ruby are averaged first, so every ruby counts once.
    Fits that are not ok are excluded. Returns (fits, per_ruby, summary):
      fits     : input table + lambda0_nm, lambda0_err_nm, P_GPa, P_err_GPa per spectrum
      per_ruby : one row per ruby (side, temperature, ruby number), with an outlier test
                 against the other rubies of the same side and temperature
      summary  : one row per temperature: reference wavelength, sample mean pressure,
                 its uncertainty, the spread between sample rubies and any outlier rubies
    """
    fits = fits.copy()
    good = fits[fits["ok"].astype(bool)].copy()
    good["chip"] = good["ruby"].where(good["ruby"] != "", good["file"])

    per_ruby = (good.groupby(["temperature_K", "side", "chip"], dropna=False, sort=True)
                .agg(n_spectra=("r1_nm", "size"), r1_mean_nm=("r1_nm", "mean"),
                     r1_std_nm=("r1_nm", "std"), r1_fit_err_nm=("r1_err_nm", "mean"),
                     files=("file", ";".join))
                .reset_index().rename(columns={"chip": "ruby"}))
    n = per_ruby["n_spectra"]
    per_ruby["r1_err_nm"] = np.where(n > 1, per_ruby["r1_std_nm"] / np.sqrt(n), per_ruby["r1_fit_err_nm"])

    # Reference wavelength for every temperature
    lam0 = {}
    for T, grp in per_ruby.groupby("temperature_K", dropna=False, sort=True):
        if ref_wavelength is not None:
            lam0[T] = (float(ref_wavelength), float(ref_err), 0)
        else:
            ref = grp[grp["side"] == "reference"]
            lam0[T] = (*_mean_and_error(ref["r1_mean_nm"], ref["r1_err_nm"]), len(ref))

    def attach(df, lam_col, err_col):
        key = df["temperature_K"]
        get = lambda i: key.map(lambda T: _lookup(lam0, T)[i])
        df["lambda0_nm"], df["lambda0_err_nm"] = get(0), get(1)
        P, dP, dP0 = pressure_gpa(df[lam_col], df["lambda0_nm"], scale)
        df["P_GPa"] = P
        df["P_err_GPa"] = np.hypot(dP * df[err_col], dP0 * df["lambda0_err_nm"])

    attach(fits, "r1_nm", "r1_err_nm")
    fits.loc[~fits["ok"].astype(bool), ["P_GPa", "P_err_GPa"]] = math.nan
    attach(per_ruby, "r1_mean_nm", "r1_err_nm")

    is_sample = per_ruby["side"] == "sample"
    if ref_wavelength is not None:                       # single-ruby mode: every non-reference file
        is_sample = per_ruby["side"] != "reference"
    per_ruby = _flag_outliers(per_ruby, is_sample)
    rows = []
    for T, (l0, l0_err, n_ref) in lam0.items():
        smp = per_ruby[is_sample & _same_T(per_ruby["temperature_K"], T)]
        lam, lam_err = _mean_and_error(smp["r1_mean_nm"], smp["r1_err_nm"])
        P, dP, dP0 = pressure_gpa(lam, l0, scale)
        spectra = good[_same_T(good["temperature_K"], T)]
        ref_spec = spectra[spectra["side"] == "reference"]["r1_nm"].to_numpy()
        smp_spec = spectra[spectra["side"] == "sample"]["r1_nm"].to_numpy()
        pairwise = (float(np.mean(pressure_gpa(smp_spec[:, None], ref_spec[None, :], scale)[0]))
                    if ref_spec.size and smp_spec.size else math.nan)
        flagged = per_ruby[_same_T(per_ruby["temperature_K"], T) & per_ruby["outlier"]]
        rows.append({
            "temperature_K": T,
            "n_reference_rubies": n_ref, "lambda0_nm": l0, "lambda0_err_nm": l0_err,
            "n_sample_rubies": len(smp), "n_sample_spectra": int(smp["n_spectra"].sum()),
            "r1_sample_nm": lam, "r1_sample_err_nm": lam_err,
            "P_GPa": float(P), "P_err_GPa": float(np.hypot(dP * lam_err, dP0 * l0_err)),
            "P_ruby_std_GPa": float(smp["P_GPa"].std(ddof=1)) if len(smp) > 1 else math.nan,
            "P_ruby_min_GPa": float(smp["P_GPa"].min()) if len(smp) else math.nan,
            "P_ruby_max_GPa": float(smp["P_GPa"].max()) if len(smp) else math.nan,
            "P_notebook_pairwise_GPa": pairwise,
            "n_bad_fits_excluded": int((~fits["ok"].astype(bool) & _same_T(fits["temperature_K"], T)).sum()),
            "outlier_rubies": "; ".join(f"{r.side or 'sample'} {r.ruby}" for r in flagged.itertuples()),
        })
    return fits, per_ruby, pd.DataFrame(rows)


def _flag_outliers(per_ruby, is_sample, level=0.995):
    """Leave-one-out check of every ruby against the other rubies of the same role and temperature.

    A ruby is flagged when its R1 lies outside the 99% prediction interval of the others:
        |x - mean_others| > t(0.995, n - 2) * sd_others * sqrt(1 + 1/(n - 1))
    This needs at least 4 rubies on that side; with fewer, nothing is tested.
    """
    from scipy.stats import t as student_t

    per_ruby = per_ruby.assign(r1_dev_from_others_nm=math.nan, outlier_t=math.nan,
                               outlier_t_limit=math.nan, outlier=False)
    role = np.where(per_ruby["side"] == "reference", "reference", np.where(is_sample, "sample", ""))
    for T in pd.unique(per_ruby["temperature_K"]):
        for r in ("reference", "sample"):
            idx = per_ruby.index[_same_T(per_ruby["temperature_K"], T) & (role == r)]
            n = len(idx)
            if n < 4:
                continue
            x = per_ruby.loc[idx, "r1_mean_nm"].to_numpy()
            limit = float(student_t.ppf(level, n - 2))
            for k, i in enumerate(idx):
                others = np.delete(x, k)
                dev = x[k] - others.mean()
                sd = others.std(ddof=1)
                t = dev / (sd * math.sqrt(1 + 1 / (n - 1))) if sd > 0 else (0.0 if dev == 0 else math.inf)
                per_ruby.loc[i, ["r1_dev_from_others_nm", "outlier_t", "outlier_t_limit"]] = [dev, t, limit]
                per_ruby.loc[i, "outlier"] = abs(t) > limit
    return per_ruby


def cross_matrix(per_ruby, scale="mao1986"):
    """Pressure of every reference x sample ruby pairing at each temperature (long format).

    The mean over all pairings at a temperature equals the summary pressure; the pattern of
    rows (reference rubies) and columns (sample rubies) shows which rubies drive the scatter.
    """
    rows = []
    for T, g in per_ruby.groupby("temperature_K", dropna=False, sort=True):
        ref, smp = g[g["side"] == "reference"], g[g["side"] == "sample"]
        for r in ref.itertuples():
            for s in smp.itertuples():
                P, dP, dP0 = pressure_gpa(s.r1_mean_nm, r.r1_mean_nm, scale)
                rows.append({
                    "temperature_K": T, "reference_ruby": r.ruby, "sample_ruby": s.ruby,
                    "lambda0_nm": r.r1_mean_nm, "lambda0_err_nm": r.r1_err_nm,
                    "r1_sample_nm": s.r1_mean_nm, "r1_sample_err_nm": s.r1_err_nm,
                    "P_GPa": float(P), "P_err_GPa": float(np.hypot(dP * s.r1_err_nm, dP0 * r.r1_err_nm)),
                })
    return pd.DataFrame(rows, columns=["temperature_K", "reference_ruby", "sample_ruby", "lambda0_nm",
                                       "lambda0_err_nm", "r1_sample_nm", "r1_sample_err_nm", "P_GPa", "P_err_GPa"])


def error_budget(fits, per_ruby, summary, cross, scale="mao1986", ref_wavelength=None):
    """Uncertainty components and R1 line widths per temperature, converted to pressure (GPa).

    *_std_GPa columns are 1-sigma spreads of single measurements (one spectrum, one ruby, one
    reference + sample pairing); *_part_GPa columns are the contributions of the reference and
    sample rubies to the uncertainty of the mean pressure (P_err_GPa, the reported +/-).
    """
    good = fits[fits["ok"].astype(bool)]
    rows = []
    for s in summary.itertuples():
        T = s.temperature_K
        _, dP, dP0 = pressure_gpa(s.r1_sample_nm, s.lambda0_nm, scale)
        dP, dP0 = abs(float(dP)), abs(float(dP0))
        pr = per_ruby[_same_T(per_ruby["temperature_K"], T)]
        sp = good[_same_T(good["temperature_K"], T)]
        if ref_wavelength is None:
            ref, smp = pr[pr["side"] == "reference"], pr[pr["side"] == "sample"]
            ref_sp, smp_sp = sp[sp["side"] == "reference"], sp[sp["side"] == "sample"]
        else:                                            # fixed reference: no reference rubies
            ref, smp = pr.iloc[:0], pr[pr["side"] != "reference"]
            ref_sp, smp_sp = sp.iloc[:0], sp[sp["side"] != "reference"]
        multi = pr[pr["n_spectra"] > 1]
        c = cross[_same_T(cross["temperature_K"], T)]
        pair_std = float(c["P_GPa"].std(ddof=1)) if len(c) > 1 else math.nan
        rows.append({
            "temperature_K": T, "P_GPa": s.P_GPa, "P_err_GPa": s.P_err_GPa,
            "fit_err_GPa": dP * smp_sp["r1_err_nm"].median(),
            "repeat_std_GPa": dP * math.sqrt((multi["r1_std_nm"] ** 2).mean()) if len(multi) else math.nan,
            "sample_ruby_std_GPa": dP * smp["r1_mean_nm"].std(ddof=1),
            "reference_ruby_std_GPa": dP0 * ref["r1_mean_nm"].std(ddof=1),
            "sample_part_GPa": dP * s.r1_sample_err_nm,
            "reference_part_GPa": dP0 * s.lambda0_err_nm,
            "single_pair_std_GPa": pair_std,
            "naive_pairwise_err_GPa": pair_std / math.sqrt(len(c)) if len(c) > 1 else math.nan,
            "n_reference_rubies": len(ref), "n_sample_rubies": len(smp), "n_pairs": len(c),
            "r1_fwhm_reference_nm": ref_sp["r1_fwhm_nm"].mean(),
            "r1_fwhm_reference_std_nm": ref_sp["r1_fwhm_nm"].std(ddof=1),
            "r1_fwhm_sample_nm": smp_sp["r1_fwhm_nm"].mean(),
            "r1_fwhm_sample_std_nm": smp_sp["r1_fwhm_nm"].std(ddof=1),
        })
    return pd.DataFrame(rows)


def _same_T(series, T):
    return series.isna() if pd.isna(T) else series == T


def _lookup(lam0, T):
    for key, val in lam0.items():
        if (pd.isna(key) and pd.isna(T)) or key == T:
            return val
    return (math.nan, math.nan, 0)


# ----------------------------------------------------------------------------------------
# Plots (matplotlib is only imported when plotting)
# ----------------------------------------------------------------------------------------

COLORS = {"data": "#8a8984", "fit": "#2a78d6", "R1": "#eb6834", "R2": "#1baf7a",
          "background": "#52514e", "grid": "#e6e5e0", "text": "#0b0b0b", "muted": "#52514e"}


def _style_axes(ax):
    ax.grid(True, color=COLORS["grid"], linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(COLORS["muted"])
    ax.tick_params(colors=COLORS["muted"], labelsize=9)


def plot_fit(fit, path=None, title=""):
    """Data, fit, R1/R2/background components and residuals. Saves to `path` if given."""
    from matplotlib.figure import Figure

    fig = Figure(figsize=(8, 6), dpi=110, layout="constrained")
    ax1, ax2 = fig.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    xf = np.linspace(fit.x[0], fit.x[-1], 1500)
    comps = fit.components(xf)
    p, e = fit.params, fit.errors
    ax1.plot(fit.x, fit.y, "o", ms=2.5, color=COLORS["data"], label="data")
    ax1.plot(xf, comps["R1"] + comps["background"], "--", lw=1.2, color=COLORS["R1"],
             label=f"R1  {p['r1_center']:.4f} ± {e['r1_center']:.4f} nm")
    ax1.plot(xf, comps["R2"] + comps["background"], "--", lw=1.2, color=COLORS["R2"],
             label=f"R2  {p['r2_center']:.4f} ± {e['r2_center']:.4f} nm")
    ax1.plot(xf, comps["background"], ":", lw=1.2, color=COLORS["background"], label="background")
    ax1.plot(xf, sum(comps.values()), lw=1.5, color=COLORS["fit"], label="fit")
    ax1.set_ylabel("Intensity (counts)")
    ax1.legend(frameon=False, fontsize=9, loc="upper left")
    status = "OK" if fit.ok else "CHECK: " + ", ".join(fit.flags)
    ax1.set_title(f"{title}\n{PROFILES[fit.profile]}   $r^2$ = {fit.r_squared:.5f}   {status}", fontsize=10, loc="left")
    ax2.plot(fit.x, fit.y - fit.model(), lw=0.8, color=COLORS["muted"])
    ax2.axhline(0, color=COLORS["text"], lw=0.6)
    for peak, color in (("r1", COLORS["R1"]), ("r2", COLORS["R2"])):   # mark the peak centres
        ax2.axvline(p[f"{peak}_center"], color=color, lw=1, ls="--")
        ax2.annotate(peak.upper(), (p[f"{peak}_center"], 1), xycoords=("data", "axes fraction"),
                     xytext=(3, -2), textcoords="offset points", va="top", fontsize=8, color=COLORS["muted"])
    ax2.set_xlabel("Wavelength (nm)")
    ax2.set_ylabel("Residual")
    for ax in (ax1, ax2):
        _style_axes(ax)
    if path:
        fig.savefig(path)
    return fig


def plot_summary(per_ruby, summary, path=None):
    """R1 positions and pressures per ruby, grouped by temperature."""
    from matplotlib.figure import Figure

    temps = list(summary["temperature_K"])
    labels = ["no T" if pd.isna(T) else f"{T:g} K" for T in temps]
    xpos = {i: i for i in range(len(temps))}
    fig = Figure(figsize=(10, 4.5), dpi=110, layout="constrained")
    ax1, ax2 = fig.subplots(1, 2)

    for i, T in enumerate(temps):
        grp = per_ruby[_same_T(per_ruby["temperature_K"], T)]
        for side, color, dx in (("reference", COLORS["fit"], -0.12), ("sample", COLORS["R1"], 0.12)):
            g = grp[grp["side"] == side]
            ax1.errorbar(np.full(len(g), xpos[i] + dx), g["r1_mean_nm"], yerr=g["r1_err_nm"], fmt="o",
                         ms=5, color=color, capsize=2, lw=1,
                         label=("reference rubies" if side == "reference" else "sample rubies") if i == 0 else None)
        smp = grp[grp["side"] == "sample"]
        jitter = np.linspace(-0.15, 0.15, len(smp)) if len(smp) > 1 else np.zeros(len(smp))
        ax2.errorbar(xpos[i] + jitter, smp["P_GPa"], yerr=smp["P_err_GPa"], fmt="o", ms=4,
                     color=COLORS["data"], capsize=2, lw=1, label="individual rubies" if i == 0 else None)
        row = summary.iloc[i]
        if np.isfinite(row["P_GPa"]):
            ax2.errorbar([xpos[i] + 0.3], [row["P_GPa"]], yerr=[row["P_err_GPa"]], fmt="D", ms=7,
                         color=COLORS["fit"], capsize=3, lw=1.5, label="mean ± error" if i == 0 else None)
            ax2.annotate(f"{1000 * row['P_GPa']:.0f} ± {1000 * row['P_err_GPa']:.0f} MPa",
                         (xpos[i] + 0.38, row["P_GPa"]), fontsize=9, color=COLORS["text"], va="center")

    for ax, ylabel in ((ax1, "R1 wavelength (nm)"), (ax2, "Pressure (GPa)")):
        ax.set_xticks(range(len(temps)), labels)
        ax.set_xlim(-0.6, len(temps) - 0.4 + 0.5 * (ax is ax2))
        ax.set_ylabel(ylabel)
        ax.legend(frameon=False, fontsize=9)
        _style_axes(ax)
    ax1.set_title("R1 position per ruby", fontsize=10, loc="left")
    ax2.set_title("Pressure per ruby", fontsize=10, loc="left")
    if path:
        fig.savefig(path)
    return fig


def _temp_label(T):
    return "no T" if pd.isna(T) else f"{T:g} K"


def _temperature_colors(n):
    """Sequential blue ramp for ordered temperatures."""
    from matplotlib.colors import LinearSegmentedColormap

    ramp = LinearSegmentedColormap.from_list("temps", ["#86b6ef", "#104281"])
    return [ramp(v) for v in (np.linspace(0, 1, n) if n > 1 else [0.6])]


def plot_cross(cross, per_ruby, path=None):
    """Heat map of the pressure of every reference x sample ruby pairing, one panel per temperature."""
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.figure import Figure

    panels, vmax = [], 0.0
    for T, c in cross.groupby("temperature_K", dropna=False, sort=True):
        M = c.pivot(index="reference_ruby", columns="sample_ruby", values="P_GPa") * 1000
        mean = float(np.nanmean(M.to_numpy()))
        vmax = max(vmax, float(np.nanmax(np.abs(M.to_numpy() - mean))))
        flagged = per_ruby[_same_T(per_ruby["temperature_K"], T) & per_ruby["outlier"]]
        panels.append((T, M, mean, set(zip(flagged["side"], flagged["ruby"]))))
    cmap = LinearSegmentedColormap.from_list("diverging", [COLORS["fit"], "#f0efec", "#e34948"])
    cmap.set_bad("#ffffff")
    fig = Figure(figsize=(max(4.4 * len(panels), 6.0), 4.6), dpi=110, layout="constrained")
    axes = np.atleast_1d(fig.subplots(1, len(panels)))
    image = None
    for ax, (T, M, mean, flagged) in zip(axes, panels):
        nr, nc = M.shape
        vals = np.full((nr + 1, nc + 1), np.nan)            # last row / column hold the means
        vals[:nr, :nc] = M.to_numpy()
        vals[:nr, nc] = np.nanmean(M.to_numpy(), axis=1)
        vals[nr, :nc] = np.nanmean(M.to_numpy(), axis=0)
        vals[nr, nc] = mean
        dev = np.full_like(vals, np.nan)
        dev[:nr, :nc] = vals[:nr, :nc] - mean
        image = ax.imshow(dev, cmap=cmap, vmin=-(vmax or 1), vmax=vmax or 1, aspect="auto")
        for i in range(nr + 1):
            for j in range(nc + 1):
                if np.isfinite(vals[i, j]):
                    ax.text(j, i, f"{vals[i, j]:.0f}", ha="center", va="center", fontsize=9, color=COLORS["text"],
                            fontweight="bold" if (i == nr or j == nc) else "normal")
        ax.axhline(nr - 0.5, color=COLORS["muted"], lw=1)
        ax.axvline(nc - 0.5, color=COLORS["muted"], lw=1)
        ax.set_xticks(range(nc + 1), [f"{s}*" if ("sample", s) in flagged else str(s) for s in M.columns] + ["mean"])
        ax.set_yticks(range(nr + 1), [f"{r}*" if ("reference", r) in flagged else str(r) for r in M.index] + ["mean"])
        ax.set_xlabel("Sample ruby")
        ax.set_ylabel("Reference ruby")
        ax.set_title(f"{_temp_label(T)}: mean of all pairings {mean:.0f} MPa", fontsize=10, loc="left")
        ax.tick_params(length=0, labelsize=9, colors=COLORS["muted"])
        for spine in ax.spines.values():
            spine.set_visible(False)
    cbar = fig.colorbar(image, ax=list(axes), shrink=0.85)
    cbar.set_label("Pairing pressure − mean (MPa)")
    note = "Pressure (MPa) for every reference × sample ruby pairing"
    if any(flagged for *_, flagged in panels):
        note += "   * outlier ruby (outside the 99% range of the other rubies)"
    fig.suptitle(note, fontsize=10, x=0.01, ha="left")
    if path:
        fig.savefig(path)
    return fig


BUDGET_ITEMS = [
    ("fit_err_GPa", "Fit only (one spectrum)"),
    ("repeat_std_GPa", "Repeat spectra of one ruby"),
    ("sample_ruby_std_GPa", "Between sample rubies"),
    ("reference_ruby_std_GPa", "Between reference rubies"),
    ("single_pair_std_GPa", "One reference + one sample ruby"),
    ("sample_part_GPa", "Mean P: from sample rubies"),
    ("reference_part_GPa", "Mean P: from reference rubies"),
    ("P_err_GPa", "Mean P: total (reported ±)"),
    ("naive_pairwise_err_GPa", "Naive: pairings as independent"),
]


def plot_error_budget(cross, budget, path=None):
    """(a) pressure of every pairing with the mean and its uncertainty; (b) uncertainty budget."""
    from matplotlib.figure import Figure

    fig = Figure(figsize=(12, 5.4), dpi=110, layout="constrained")
    ax1, ax2 = fig.subplots(1, 2, gridspec_kw={"width_ratios": [1, 1.35]})
    colors = _temperature_colors(len(budget))
    for i, b in enumerate(budget.itertuples()):
        c = cross[_same_T(cross["temperature_K"], b.temperature_K)]
        y = c["P_GPa"].to_numpy() * 1000
        jitter = np.linspace(-0.2, 0.2, len(y)) if len(y) > 1 else np.zeros(len(y))
        ax1.plot(i - 0.12 + jitter, y, "o", ms=4, color=COLORS["data"], alpha=0.8,
                 label="each pairing" if i == 0 else None)
        if np.isfinite(b.P_GPa):
            if np.isfinite(b.single_pair_std_GPa):
                ax1.errorbar([i + 0.22], [1000 * b.P_GPa], yerr=[1000 * b.single_pair_std_GPa], fmt="none",
                             ecolor=COLORS["R1"], elinewidth=1.5, capsize=4,
                             label="± spread of single pairings" if i == 0 else None)
            ax1.errorbar([i + 0.32], [1000 * b.P_GPa], yerr=[1000 * b.P_err_GPa], fmt="D", ms=7,
                         color=COLORS["fit"], capsize=4, lw=1.5, label="mean ± reported uncertainty" if i == 0 else None)
            ax1.annotate(f"{1000 * b.P_GPa:.0f} ± {1000 * b.P_err_GPa:.0f}", (i + 0.42, 1000 * b.P_GPa),
                         fontsize=9, va="center", color=COLORS["text"])
    ax1.set_xticks(range(len(budget)), [_temp_label(T) for T in budget["temperature_K"]])
    ax1.set_xlim(-0.6, len(budget) - 0.4 + 0.5)
    ax1.set_ylabel("Pressure (MPa)")
    ax1.set_title("Pressure from every reference × sample pairing", fontsize=10, loc="left")
    ax1.legend(frameon=False, fontsize=9)

    labels = [label for _, label in BUDGET_ITEMS]
    y = np.arange(len(labels))
    h = 0.8 / max(len(budget), 1)
    for k, b in enumerate(budget.itertuples()):
        vals = np.array([getattr(b, key) for key, _ in BUDGET_ITEMS], dtype=float) * 1000
        pos = y - 0.4 + h * (k + 0.5)
        ax2.barh(pos, vals, height=h * 0.88, color=colors[k], label=_temp_label(b.temperature_K))
        for yy, v in zip(pos, vals):
            if np.isfinite(v):
                ax2.text(v, yy, f" {v:.1f}", va="center", fontsize=7.5, color=COLORS["muted"])
    ax2.set_yticks(y, labels)
    ax2.invert_yaxis()
    for sep in (4.5, 7.5):
        ax2.axhline(sep, color=COLORS["grid"], lw=1)
    ax2.set_xlabel("1σ (MPa)")
    ax2.set_title("Uncertainty budget: single measurements (top) → mean pressure (bottom)", fontsize=10, loc="left")
    ax2.legend(frameon=False, fontsize=9, loc="lower right")
    ax2.margins(x=0.12)
    for ax in (ax1, ax2):
        _style_axes(ax)
    ax2.grid(False, axis="y")
    if path:
        fig.savefig(path)
    return fig


SPREAD_ITEMS = [
    ("reference_ruby_std_GPa", "Between reference rubies", COLORS["fit"]),
    ("sample_ruby_std_GPa", "Between sample rubies", COLORS["R1"]),
    ("repeat_std_GPa", "Repeat spectra of one ruby", COLORS["R2"]),
    ("fit_err_GPa", "Fit only (one spectrum)", "#eda100"),
]


def plot_spread(budget, path=None):
    """Scatter of the R1 measurements and R1 line width as a function of temperature and pressure."""
    from matplotlib.figure import Figure

    fig = Figure(figsize=(11, 7.6), dpi=110, layout="constrained")
    axes = fig.subplots(2, 2)
    xs = [(budget["temperature_K"].to_numpy(float), "Temperature (K)"),
          (budget["P_GPa"].to_numpy(float) * 1000, "Pressure (MPa)")]
    for col, (x, xlabel) in enumerate(xs):
        order = np.argsort(x)
        ax = axes[0, col]
        for key, label, color in SPREAD_ITEMS:
            ax.plot(x[order], budget[key].to_numpy(float)[order] * 1000, "-o", color=color, lw=1.5, ms=6, label=label)
        ax.set_ylim(bottom=0)
        ax.set_ylabel("1σ spread (MPa)")
        ax.set_xlabel(xlabel)
        ax.set_title(f"Scatter of R1 measurements vs {xlabel.split()[0].lower()}", fontsize=10, loc="left")
        ax = axes[1, col]
        for side, label, color in (("reference", "Reference rubies", COLORS["fit"]), ("sample", "Sample rubies", COLORS["R1"])):
            ax.errorbar(x[order], budget[f"r1_fwhm_{side}_nm"].to_numpy(float)[order],
                        yerr=budget[f"r1_fwhm_{side}_std_nm"].to_numpy(float)[order],
                        fmt="-o", color=color, lw=1.5, ms=6, capsize=3, label=label)
        ax.set_ylabel("R1 FWHM (nm), mean ± sd")
        ax.set_xlabel(xlabel)
        ax.set_title(f"R1 line width vs {xlabel.split()[0].lower()}", fontsize=10, loc="left")
    axes[0, 0].legend(frameon=False, fontsize=9)
    axes[1, 0].legend(frameon=False, fontsize=9)
    for ax in axes.ravel():
        _style_axes(ax)
    if path:
        fig.savefig(path)
    return fig


# ----------------------------------------------------------------------------------------
# Command line
# ----------------------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description="Fit ruby R1/R2 spectra and compute pressures.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="files, folders or glob patterns")
    ap.add_argument("-o", "--out", default="rubyfit_results", help="output folder")
    ap.add_argument("--pattern", default="*.csv", help="file pattern used inside folders")
    ap.add_argument("--window", nargs=2, type=float, default=DEFAULT_WINDOW, metavar=("LO", "HI"),
                    help="fit window relative to the R1 peak (nm)")
    ap.add_argument("--full-range", action="store_true", help="fit the whole spectrum (as the notebook)")
    ap.add_argument("--search", nargs=2, type=float, metavar=("LO", "HI"),
                    help="only look for R1 between these wavelengths (nm)")
    ap.add_argument("--min-r2", type=float, default=DEFAULT_MIN_R2, help="flag fits below this r-squared")
    ap.add_argument("--profile", choices=sorted(PROFILES), default=DEFAULT_PROFILE,
                    help="peak line shape (split_pearson7: asymmetric, center = peak maximum)")
    ap.add_argument("--scale", choices=sorted(SCALES), default="mao1986", help="ruby pressure scale")
    ap.add_argument("--ref-wavelength", type=float,
                    help="fixed reference R1 (nm), e.g. 694.22; otherwise reference rubies are used")
    ap.add_argument("--ref-err", type=float, default=0.02, help="uncertainty of --ref-wavelength (nm)")
    ap.add_argument("--ref-keys", default=",".join(DEFAULT_REF_KEYS),
                    help="file-name words marking reference rubies")
    ap.add_argument("--sample-keys", default=",".join(DEFAULT_SAMPLE_KEYS),
                    help="file-name words marking sample rubies")
    ap.add_argument("--plots", choices=("none", "flagged", "all"), default="flagged",
                    help="which per-file fit plots to save")
    ap.add_argument("-j", "--jobs", type=int, default=1, help="parallel processes")
    args = ap.parse_args(argv)

    paths = find_files(args.inputs, args.pattern)
    if not paths:
        sys.exit("no input files found")
    os.makedirs(args.out, exist_ok=True)

    t0 = time.perf_counter()
    fits = fit_files(paths, jobs=args.jobs,
                     plot_dir=os.path.join(args.out, "plots") if args.plots != "none" else None,
                     plot_mode=args.plots,
                     ref_keys=tuple(k.strip() for k in args.ref_keys.split(",") if k.strip()),
                     sample_keys=tuple(k.strip() for k in args.sample_keys.split(",") if k.strip()),
                     window=None if args.full_range else tuple(args.window),
                     search=tuple(args.search) if args.search else None,
                     min_r2=args.min_r2, profile=args.profile)
    t_fit = time.perf_counter() - t0
    fits, per_ruby, summary = pressure_tables(fits, args.scale, args.ref_wavelength, args.ref_err)

    # the cross matrix only makes sense with reference rubies (dual-ruby mode)
    cross = cross_matrix(per_ruby if args.ref_wavelength is None else per_ruby.iloc[:0], args.scale)
    budget = error_budget(fits, per_ruby, summary, cross, args.scale, args.ref_wavelength)

    fits.drop(columns="path").to_csv(os.path.join(args.out, "fits.csv"), index=False)
    per_ruby.to_csv(os.path.join(args.out, "pressure_by_ruby.csv"), index=False)
    summary.to_csv(os.path.join(args.out, "pressure_summary.csv"), index=False)
    budget.to_csv(os.path.join(args.out, "error_budget.csv"), index=False)
    if len(cross):
        cross.to_csv(os.path.join(args.out, "pressure_cross.csv"), index=False)
        plot_cross(cross, per_ruby, os.path.join(args.out, "cross_matrix.png"))
    if len(summary):
        plot_summary(per_ruby, summary, os.path.join(args.out, "summary.png"))
        plot_error_budget(cross, budget, os.path.join(args.out, "error_budget.png"))
        plot_spread(budget, os.path.join(args.out, "spread_vs_T_P.png"))

    bad = fits[~fits["ok"].astype(bool)]
    print(f"Fitted {len(fits)} files in {t_fit:.2f} s ({len(fits) - len(bad)} ok, {len(bad)} flagged)")
    for _, r in bad.iterrows():
        print(f"  flagged: {r['file']}: {r['flags']}")
    print(f"Line shape: {PROFILES[args.profile]}; pressure scale: {SCALES[args.scale]}")
    for _, r in summary.iterrows():
        T = "no T" if pd.isna(r["temperature_K"]) else f"{r['temperature_K']:g} K"
        if np.isfinite(r["P_GPa"]):
            spread = (f", spread between rubies {1000 * r['P_ruby_std_GPa']:.0f} MPa"
                      if np.isfinite(r["P_ruby_std_GPa"]) else "")
            source = (f"from {int(r['n_reference_rubies'])} reference rubies" if args.ref_wavelength is None
                      else "fixed")
            print(f"  {T:>7}: P = {r['P_GPa']:.4f} ± {r['P_err_GPa']:.4f} GPa "
                  f"({1000 * r['P_GPa']:.0f} ± {1000 * r['P_err_GPa']:.0f} MPa){spread}; "
                  f"lambda0 = {r['lambda0_nm']:.4f} nm ({source}), {int(r['n_sample_rubies'])} sample rubies")
        else:
            print(f"  {T:>7}: no pressure (missing sample or reference rubies at this temperature)")
    for b in budget.itertuples():
        if np.isfinite(b.P_err_GPa):
            print(f"  {_temp_label(b.temperature_K):>7}  uncertainty: reference rubies ±{1000 * b.reference_part_GPa:.0f}, "
                  f"sample rubies ±{1000 * b.sample_part_GPa:.0f} MPa; one ref + one sample ruby ±"
                  f"{1000 * b.single_pair_std_GPa:.0f}, one spectrum (fit only) ±{1000 * b.fit_err_GPa:.1f} MPa")
    for r in per_ruby[per_ruby["outlier"]].itertuples():
        print(f"  ⚠ outlier: {r.side or 'sample'} ruby {r.ruby} at {_temp_label(r.temperature_K)}: R1 "
              f"{r.r1_dev_from_others_nm:+.4f} nm from the other {r.side or 'sample'} rubies "
              f"(t = {r.outlier_t:.1f}, 99% limit {r.outlier_t_limit:.1f})")
    if len(cross):
        print("Cross-check, pressure (MPa) for every reference (rows) × sample (columns) ruby:")
        for T, c in cross.groupby("temperature_K", dropna=False, sort=True):
            M = c.pivot(index="reference_ruby", columns="sample_ruby", values="P_GPa") * 1000
            M.index.name, M.columns.name = "ref", "sample"
            print(f"  {_temp_label(T)}")
            print("    " + M.to_string(float_format=lambda v: f"{v:.0f}").replace("\n", "\n    "))
    print(f"Results written to {os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()
