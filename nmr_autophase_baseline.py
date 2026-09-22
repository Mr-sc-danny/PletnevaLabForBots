"""
nmr_autophase_baseline.py

Automated phase correction + baseline correction for 1D NMR spectra exported
from TopSpin as ASCII text (real-only or real+imaginary), with the
optimization deliberately focused on a negative-ppm region and an extra
emphasis window around a chosen chemical shift (default: -3 ppm).

No third-party packages beyond numpy / scipy / pandas / matplotlib are used.

------------------------------------------------------------------------
WHAT IT DOES
------------------------------------------------------------------------
1. Parses a TopSpin ASCII export (the "# LEFT = ... RIGHT = ..." /
   "# SIZE = ..." header format), handling both real-only traces and
   real+imaginary complex traces (e.g. "12345.6-789.0i"), including
   scientific notation.

2. Phase correction:
   - Applies a single global zero-order (PH0) + first-order (PH1) phase
     to the WHOLE spectrum (phase is a property of the whole spectrum,
     not just the region you care about).
   - Optimizes PH0/PH1 by minimizing an ACME-style cost (entropy of the
     first derivative + a penalty for negative-going signal), evaluated
     ONLY on the focus region (ppm < focus_max, default 0 ppm).
   - Adds an *extra* negativity penalty inside a Gaussian-weighted
     emphasis window around `emphasis_center` (default -3 ppm), so the
     optimizer works especially hard to make that specific region look
     clean and absorptive, even if it has to compromise slightly
     elsewhere in the focus region.
   - Multi-start Nelder-Mead search to avoid local minima.

3. Baseline correction (focus region only):
   - Fits a smooth baseline with Asymmetric Least Squares (ALS), which
     naturally sits at/below broad humps and peak tails without eating
     into sharp peaks.
   - ALS's asymmetric penalty is known to bias the fitted baseline
     slightly low. This script auto-detects peaks (scipy.signal.find_peaks
     on the initial ALS residual) and estimates a constant, peak-free
     residual offset that's added back so the *quiet* regions genuinely
     center on zero (not just "close to it").

4. Outputs:
   - A CSV with ppm, raw real/imag, phased real/imag, fitted baseline,
     and final phased+baseline-corrected intensity (only for the
     focus region, since that's what's being corrected).
   - Diagnostic PNGs: full-spectrum before/after phasing, and a
     zoomed/baseline-corrected view of the focus region.

------------------------------------------------------------------------
USAGE
------------------------------------------------------------------------
    python nmr_autophase_baseline.py input.txt --outdir ./out

Optional flags (all have sensible defaults):
    --focus-max        ppm boundary of the focus region (default 0.0;
                        the optimizer only "looks at" ppm < focus_max)
    --emphasis-center   ppm center of the extra-weighted window (default -3.0)
    --emphasis-width     half-width in ppm of the emphasis window (default 0.3)
    --emphasis-gain      how much extra negativity penalty in that window
                          relative to the rest of the focus region (default 4.0)
    --als-lam             ALS smoothness parameter (default 5e4)
    --als-p                ALS asymmetry parameter (default 1e-3)

You can also `import` this file and call `process_spectrum(...)` directly.
"""

import argparse
import os
import re

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.optimize import minimize
from scipy.signal import find_peaks, peak_widths

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------
# 1. Parsing
# ----------------------------------------------------------------------

_NUM = r'[+-]?\d+\.?\d*(?:[Ee][+-]?\d+)?'
_COMPLEX_RE = re.compile(rf'^({_NUM})({_NUM})i$')


def load_topspin_ascii(path):
    """
    Load a TopSpin ASCII export (real-only or real+imaginary).

    Returns
    -------
    spec : complex ndarray
        Complex spectrum (imag=0 if the file was real-only).
    ppm : ndarray
        ppm axis, same length as spec, ordered exactly as the file
        (left ppm value first).
    """
    with open(path) as f:
        lines = f.readlines()

    left_ppm = right_ppm = size = None
    for l in lines:
        if l.startswith('# LEFT') or ('LEFT' in l and 'RIGHT' in l):
            m = re.search(rf'LEFT\s*=\s*({_NUM})\s*ppm\.\s*RIGHT\s*=\s*({_NUM})\s*ppm', l)
            if m:
                left_ppm, right_ppm = float(m.group(1)), float(m.group(2))
        if l.startswith('# SIZE'):
            m = re.search(r'SIZE\s*=\s*(\d+)', l)
            if m:
                size = int(m.group(1))

    if left_ppm is None or right_ppm is None or size is None:
        raise ValueError("Could not find LEFT/RIGHT/SIZE header in file: " + path)

    data_lines = [l.strip() for l in lines if l.strip() and not l.startswith('#')]
    if len(data_lines) != size:
        # Not fatal -- just trust what's actually there.
        size = len(data_lines)

    re_vals = np.zeros(size)
    im_vals = np.zeros(size)
    is_complex = _COMPLEX_RE.match(data_lines[0]) is not None

    for i, dl in enumerate(data_lines):
        if is_complex:
            m = _COMPLEX_RE.match(dl)
            if m is None:
                raise ValueError(f"Line {i} does not match expected complex format: {dl!r}")
            re_vals[i] = float(m.group(1))
            im_vals[i] = float(m.group(2))
        else:
            re_vals[i] = float(dl)

    spec = re_vals + 1j * im_vals
    ppm = np.linspace(left_ppm, right_ppm, size)
    return spec, ppm


# ----------------------------------------------------------------------
# 2. Phase correction
# ----------------------------------------------------------------------

def apply_phase(spec, phc0_deg, phc1_deg):
    """Apply a global zero+first order phase correction (pivot = center)."""
    N = len(spec)
    n = np.arange(N)
    # ph is the correction angle at each point n. phc0 shifts every point by
    # the same amount; phc1 adds a term that grows linearly with distance
    # from the center point (n - N/2), which is why it's called "first
    # order" -- it's the linear (in index/frequency) part of the phase error.
    ph = np.deg2rad(phc0_deg) + np.deg2rad(phc1_deg) * (n - N / 2) / N
    # Multiplying a complex number by exp(i*ph) rotates it by angle ph
    # without changing its magnitude -- this is the actual "phasing" step.
    # Done correctly, it rotates each point so the real part is all
    # absorptive (Lorentzian) and the dispersive character moves into the
    # (normally discarded) imaginary part.
    return spec * np.exp(1j * ph)


def _autophase_cost(params, spec, ppm, focus_mask, emphasis_mask,
                     gamma=8e-5, emphasis_gain=4.0):
    phc0, phc1 = params
    # Try this candidate phase and look only at the region we care about.
    corrected = apply_phase(spec, phc0, phc1)
    real = corrected[focus_mask].real

    # --- Term 1: entropy of the first derivative ("ACME" trick) ---
    # A well-phased spectrum is mostly flat baseline with a few sharp
    # jumps at peak edges -- i.e. the derivative is concentrated in a few
    # large values. A badly-phased (dispersive/wavy) spectrum smears
    # energy across many small derivative values everywhere. Normalizing
    # |d| into a probability distribution `p` and computing its Shannon
    # entropy gives a single number that's LOW when the signal is
    # concentrated (good phase) and HIGH when it's smeared (bad phase).
    d = np.diff(real)
    d_abs = np.abs(d)
    s = np.sum(d_abs)
    if s < 1e-9:
        return 1e6  # degenerate/all-zero case -- heavily penalize
    p = d_abs / s
    p = p[p > 1e-12]  # avoid log(0)
    entropy = -np.sum(p * np.log(p))

    # --- Term 2: negativity penalty ---
    # Real absorptive NMR peaks point up, not down. Squaring and summing
    # the negative-going part of the signal punishes dispersive
    # (S-shaped, partially negative) lineshapes anywhere in the focus
    # region. This term alone would not find good phase (a flat
    # all-zero spectrum would "win"), which is why it's combined with
    # the entropy term above rather than used by itself.
    neg = np.clip(real, None, 0)  # keeps negative values, zeroes out positive ones
    penalty = gamma * np.sum(neg ** 2) / len(real)

    # --- Term 3: extra-weighted negativity penalty near emphasis_center ---
    # Identical idea to Term 2, but computed only inside the emphasis
    # window and multiplied by emphasis_gain. This is what makes the
    # optimizer prioritize a clean, absorptive lineshape specifically
    # around e.g. -3 ppm, even at a small cost to entropy/penalty
    # elsewhere in the focus region.
    emphasis_penalty = 0.0
    if emphasis_mask is not None and emphasis_mask.sum() > 5:
        real_emph = corrected[emphasis_mask].real
        neg_emph = np.clip(real_emph, None, 0)
        emphasis_penalty = emphasis_gain * gamma * np.sum(neg_emph ** 2) / len(real_emph)

    # Lower total cost = better phase. minimize() will search (phc0, phc1)
    # to push this number as low as possible.
    return entropy + penalty + emphasis_penalty


def autophase(spec, ppm, focus_max=0.0, emphasis_center=-3.0,
              emphasis_width=0.3, emphasis_gain=4.0, gamma=8e-5,
              p0_grid=None, p1_grid=None):
    """
    Find the global (PH0, PH1) that best phases the focus region
    (ppm < focus_max), with extra weight on the emphasis window.
    """
    focus_mask = ppm < focus_max
    if focus_mask.sum() < 20:
        raise ValueError("Focus region has too few points -- check focus_max.")

    emphasis_mask = (ppm > emphasis_center - emphasis_width) & \
                    (ppm < emphasis_center + emphasis_width)

    if p0_grid is None:
        p0_grid = np.linspace(-180, 180, 13)
    if p1_grid is None:
        p1_grid = np.linspace(-720, 720, 9)

    best = None
    for p0 in p0_grid:
        for p1 in p1_grid:
            res = minimize(
                _autophase_cost, x0=[p0, p1],
                args=(spec, ppm, focus_mask, emphasis_mask, gamma, emphasis_gain),
                method='Nelder-Mead',
                options={'xatol': 1e-3, 'fatol': 1e-8, 'maxiter': 2000},
            )
            if best is None or res.fun < best.fun:
                best = res

    phc0, phc1 = best.x
    return phc0, phc1, best.fun


# ----------------------------------------------------------------------
# 3. Baseline correction
# ----------------------------------------------------------------------

def baseline_als(y, lam=5e4, p=1e-3, niter=15):
    """Asymmetric Least Squares baseline (Eilers)."""
    y = np.asarray(y, dtype=np.float64)
    L = len(y)
    # D encodes a discrete second-derivative operator. D.dot(D.T), scaled
    # by `lam`, becomes a "roughness penalty" matrix: large values force
    # the fitted baseline z to have a small second derivative, i.e. to be
    # smooth/gently curving rather than following every wiggle in y.
    # Bigger lam -> smoother (less peak-sensitive) baseline.
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2))
    D = lam * D.dot(D.transpose())
    w = np.ones(L)          # per-point fit weights, updated each iteration
    W = sparse.spdiags(w, 0, L, L)
    z = y.copy()
    for _ in range(niter):
        # Solve the weighted, smoothness-penalized least-squares system
        # (W + D) z = W y  for the baseline z.
        W.setdiag(w)
        Z = (W + D).tocsc()
        z = spsolve(Z, w * y)
        # Re-weight: points where y sits ABOVE the current baseline guess
        # (likely peaks) get a tiny weight `p`, so they barely pull the
        # next fit upward. Points where y sits BELOW get a large weight
        # `1-p`, so the baseline is pulled down to hug the valleys.
        # Iterating this converges to a curve that tracks the low
        # envelope between peaks while mostly ignoring the peaks
        # themselves. Smaller p -> baseline sits closer to the valley
        # floor (better peak preservation, but more prone to sitting
        # systematically low -- see the offset correction in
        # baseline_correct() below).
        w = p * (y > z) + (1 - p) * (y < z)
    return z


def detect_peak_mask(residual, prominence_sigma=3.0):
    """
    Auto-detect peak regions in a (roughly baseline-flat) residual trace,
    so we know which points are safe to use for estimating a flat-region
    offset. Returns a boolean mask, True = inside a detected peak.
    """
    noise = np.std(residual[np.abs(residual) < np.percentile(np.abs(residual), 60)])
    if noise <= 0:
        noise = np.std(residual) + 1e-9
    peaks, props = find_peaks(residual, prominence=prominence_sigma * noise)
    mask = np.zeros(len(residual), dtype=bool)
    if len(peaks) == 0:
        return mask
    widths, _, left_ips, right_ips = peak_widths(residual, peaks, rel_height=0.9)
    for lip, rip in zip(left_ips, right_ips):
        lo = max(0, int(np.floor(lip)))
        hi = min(len(residual), int(np.ceil(rip)) + 1)
        mask[lo:hi] = True
    return mask


def baseline_correct(spec_real, lam=5e4, p=1e-3, niter=15):
    """
    Fit + apply ALS baseline correction, with the auto peak-free-offset
    recentering step described in the module docstring.

    Returns (baseline, corrected).
    """
    baseline = baseline_als(spec_real, lam=lam, p=p, niter=niter)
    residual = spec_real - baseline

    # ALS's small p means the fitted baseline tracks the *valley floor*
    # of the noise, not its mean -- so after subtraction, the "empty"
    # (peak-free) regions sit systematically ABOVE zero rather than
    # centered on it. Fix: look only at points NOT inside a detected
    # peak, take their median (should be ~0 if unbiased), and shift the
    # baseline up by that amount. This recenters the noise floor to a
    # true zero without touching peak heights.
    peak_mask = detect_peak_mask(residual)
    quiet = residual[~peak_mask]
    offset = np.median(quiet) if len(quiet) > 20 else 0.0

    baseline_final = baseline + offset
    corrected = spec_real - baseline_final
    return baseline_final, corrected


# ----------------------------------------------------------------------
# 4. End-to-end driver
# ----------------------------------------------------------------------

def process_spectrum(path, outdir, focus_max=0.0, emphasis_center=-3.0,
                      emphasis_width=0.3, emphasis_gain=4.0,
                      als_lam=5e4, als_p=1e-3):
    os.makedirs(outdir, exist_ok=True)
    base = os.path.splitext(os.path.basename(path))[0]

    spec, ppm = load_topspin_ascii(path)

    print(f"[{base}] loaded {len(spec)} points, ppm {ppm[0]:.3f} to {ppm[-1]:.3f}")

    phc0, phc1, cost = autophase(
        spec, ppm, focus_max=focus_max,
        emphasis_center=emphasis_center, emphasis_width=emphasis_width,
        emphasis_gain=emphasis_gain,
    )
    print(f"[{base}] best phase: PH0={phc0:.2f} deg, PH1={phc1:.2f} deg (cost={cost:.4f})")

    corrected = apply_phase(spec, phc0, phc1)

    focus_mask = ppm < focus_max
    focus_ppm = ppm[focus_mask]
    focus_real = corrected[focus_mask].real
    focus_imag = corrected[focus_mask].imag

    baseline_final, final_corrected = baseline_correct(
        focus_real, lam=als_lam, p=als_p,
    )

    # ---- Save CSV (intermediate data) ----
    df = pd.DataFrame({
        "ppm": focus_ppm,
        "real_raw": spec[focus_mask].real,
        "imag_raw": spec[focus_mask].imag,
        "real_phased": focus_real,
        "imag_phased": focus_imag,
        "baseline_fit": baseline_final,
        "real_phased_baseline_corrected": final_corrected,
    })
    csv_path = os.path.join(outdir, f"{base}_phased_baseline_corrected.csv")
    df.to_csv(csv_path, index=False)
    print(f"[{base}] wrote {csv_path}")

    # ---- Diagnostic plots ----
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    axes[0].plot(ppm, spec.real, lw=0.6, color='tab:red')
    axes[0].axhline(0, color='gray', lw=0.5)
    axes[0].axvline(focus_max, color='k', lw=0.5, ls='--')
    axes[0].set_title(f"{base}: BEFORE phasing (full spectrum)")
    axes[1].plot(ppm, corrected.real, lw=0.6, color='tab:blue')
    axes[1].axhline(0, color='gray', lw=0.5)
    axes[1].axvline(focus_max, color='k', lw=0.5, ls='--')
    axes[1].set_title(f"PH0={phc0:.1f} deg, PH1={phc1:.1f} deg  |  AFTER phasing")
    axes[1].set_xlabel("ppm")
    axes[1].invert_xaxis()
    plt.tight_layout()
    full_png = os.path.join(outdir, f"{base}_full_spectrum_before_after.png")
    plt.savefig(full_png, dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(focus_ppm, final_corrected, lw=0.6, color='tab:green')
    ax.axhline(0, color='gray', lw=0.6)
    ax.axvline(emphasis_center, color='purple', lw=0.8, ls=':',
               label=f"emphasis center ({emphasis_center} ppm)")
    ax.invert_xaxis()
    ax.set_xlabel("ppm")
    ax.set_title(f"{base}: focus region, phased + baseline-corrected")
    ax.legend(fontsize=8)
    plt.tight_layout()
    focus_png = os.path.join(outdir, f"{base}_focus_region_corrected.png")
    plt.savefig(focus_png, dpi=140)
    plt.close(fig)

    print(f"[{base}] wrote {full_png}")
    print(f"[{base}] wrote {focus_png}")

    return {
        "phc0": phc0, "phc1": phc1, "cost": cost,
        "csv_path": csv_path, "full_png": full_png, "focus_png": focus_png,
        "ppm": focus_ppm, "corrected": final_corrected,
    }


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input", help="TopSpin ASCII export (.txt)")
    ap.add_argument("--outdir", default="./nmr_out")
    ap.add_argument("--focus-max", type=float, default=0.0)
    ap.add_argument("--emphasis-center", type=float, default=-3.0)
    ap.add_argument("--emphasis-width", type=float, default=0.3)
    ap.add_argument("--emphasis-gain", type=float, default=4.0)
    ap.add_argument("--als-lam", type=float, default=5e4)
    ap.add_argument("--als-p", type=float, default=1e-3)
    args = ap.parse_args()

    process_spectrum(
        args.input, args.outdir,
        focus_max=args.focus_max,
        emphasis_center=args.emphasis_center,
        emphasis_width=args.emphasis_width,
        emphasis_gain=args.emphasis_gain,
        als_lam=args.als_lam,
        als_p=args.als_p,
    )


if __name__ == "__main__":
    main()
