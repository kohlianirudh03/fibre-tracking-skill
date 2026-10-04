"""Image quality for fibre volumes: true CNR, contrast transfer, regime check.

The conventional contrast-to-noise ratio estimates noise from the variance
within each phase. In a fibre composite the dominant contribution to that
variance is the material's own texture -- real, reproducible structure, not
detector noise. It read a flat 2.6-2.9 across ten reconstructions spanning a
20x range of voxel size and two acquisition protocols, and concealed a 29.8%
protocol improvement that it reported as exactly zero.

Noise measured from adjacent-slice differences isolates the noise instead: over
one slice spacing the microstructure is nearly identical while the noise is
independent.

  python quality.py --dir path/to/slices --voxel 0.98 --fibre-d 7.0
"""
import argparse
import glob
import os
import numpy as np
from scipy import ndimage
from skimage import filters

import fibretrack as FT


def contrast_and_noise(slices, mask_threshold=None, fibre_d_vox=None,
                       smooth=3.0):
    """Contrast, noise SD, apparent CNR and true CNR from a slice stack.

    slices : list/array of 2D arrays, consecutive in z and evenly spaced
    """
    sl = [np.asarray(s, np.float32) for s in slices]
    if len(sl) < 2:
        raise ValueError("need at least two consecutive slices for the "
                         "adjacent-slice noise estimate")

    contrasts, apparent = [], []
    for s in sl:
        if mask_threshold is not None and fibre_d_vox is not None:
            m = FT.specimen_mask(s, mask_threshold, fibre_d_vox, smooth=smooth)
        else:
            m = np.ones(s.shape, bool)
        if m.sum() < 400:
            continue
        t = filters.threshold_otsu(s[m])
        hi, lo = s[m & (s > t)], s[m & (s <= t)]
        if hi.size < 50 or lo.size < 50:
            continue
        contrasts.append(hi.mean() - lo.mean())
        # the conventional metric, kept only so the difference is visible
        apparent.append((hi.mean() - lo.mean()) /
                        np.sqrt((hi.var() + lo.var()) / 2.0))
    if not contrasts:
        raise ValueError("no slice yielded a usable fibre/matrix split; check "
                         "the mask threshold")

    diffs = [np.std(sl[i + 1] - sl[i]) / np.sqrt(2.0) for i in range(len(sl) - 1)]
    noise_sd = float(np.mean(diffs))
    contrast = float(np.mean(contrasts))
    return {"contrast": contrast,
            "noise_sd": noise_sd,
            "apparent_cnr": float(np.mean(apparent)),
            "true_cnr": contrast / noise_sd if noise_sd > 0 else float("inf"),
            "n_slices": len(sl)}


def contrast_transfer(per_voxel_contrast):
    """Fraction of finest-resolution contrast surviving at each voxel size.

    per_voxel_contrast : {voxel_um: contrast}
    This curve is the target the phantom's blur is fitted to. It is measured on
    real data and is independent of any orientation statistic, which is what
    keeps the later bias correction non-circular.
    """
    if not per_voxel_contrast:
        return {}
    ref = per_voxel_contrast[min(per_voxel_contrast)]
    return {v: c / ref for v, c in sorted(per_voxel_contrast.items())}


def _load(directory, n, start_frac=0.45):
    import tifffile
    files = sorted(glob.glob(os.path.join(directory, "*.tif*")))
    if not files:
        raise SystemExit("no TIFFs in %s" % directory)
    i0 = int(len(files) * start_frac)
    return [tifffile.imread(f) for f in files[i0:i0 + n]], len(files)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="directory of TIFF slices")
    ap.add_argument("--voxel", type=float, required=True, help="um")
    ap.add_argument("--fibre-d", type=float, required=True, help="um")
    ap.add_argument("--threshold", type=float,
                    help="specimen threshold (grey level); omit to use the "
                         "whole frame, valid only if there is no air in it")
    ap.add_argument("--n", type=int, default=8, help="consecutive slices")
    a = ap.parse_args()

    g = FT.regime(a.voxel, a.fibre_d)
    print("r = %.3f  (%s)" % (g["r"], g["regime"]))
    print("supports: %s\n" % g["supports"])

    slices, total = _load(a.dir, a.n)
    print("read %d of %d slices from %s" % (len(slices), total, a.dir))
    q = contrast_and_noise(slices, a.threshold, a.fibre_d / a.voxel)

    print("\ncontrast            %10.1f" % q["contrast"])
    print("noise SD            %10.1f   (adjacent-slice difference)" % q["noise_sd"])
    print("apparent CNR        %10.2f   <- conventional, counts texture as noise"
          % q["apparent_cnr"])
    print("true CNR            %10.2f   <- use this one" % q["true_cnr"])
    ratio = q["true_cnr"] / max(q["apparent_cnr"], 1e-9)
    print("\nthe conventional metric understates CNR by %.1fx here." % ratio)
    if not g["per_fibre_valid"]:
        print("\nr > 0.3: do not emit a per-fibre table from this scan. Write a")
        print("voxel-wise orientation field, and correct the mean angle against")
        print("a calibrated phantom (calibrate.py) before quoting it.")
