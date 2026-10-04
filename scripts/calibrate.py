"""Image-formation calibration and resolution-bias correction.

Needed whenever r > 0.3, where partial-volume averaging inflates the measured
misalignment. On one ten-scan series the measured mean angle rose from 1.00 deg
at r = 0.05 to 22.74 deg at r = 1.00 on an unchanged specimen; correcting each
scan against a phantom calibrated to that scan's own contrast transfer and CNR
collapsed the 23-fold spread to 3-fold.

Two properties must match at every voxel size: how much contrast survives, and
how much noise sits on it. Three attempts were needed.

  1. Supplying the real absolute noise.  Phantom contrast is fixed by
     construction while real contrast falls with coarsening, so coarse
     phantoms came out up to 13x too clean.
  2. Matching CNR using the real contrast.  The measured contrast at coarse
     voxels is what survived partial-volume averaging, not a material
     property. Dividing by it double-counts the loss: 5.3x too noisy.
  3. Fit blur to the contrast-transfer curve, then set noise from the
     PHANTOM's own observed contrast and the real true CNR.  This matches.

Blur is fitted to contrast transfer and noise to CNR. Neither target is the
orientation angle. Tuning against the quantity being corrected would make the
bias estimate circular and the correction worthless.

  python calibrate.py --fibre-d 5.0 --voxels 0.25 0.5 1 2.5 5 \
      --real-contrast 4140 3840 2340 650 461 --real-cnr 10.4 13.8 10.1 7.2 2.9
"""
import argparse
import json
import numpy as np
from skimage import filters

import phantom as PH

BLURS = [1.0, 1.4, 1.8, 2.2, 2.8, 3.4, 4.2]


def observed_contrast(vol, n=8):
    """Fibre-matrix contrast actually present in a generated volume."""
    out = []
    for z in np.linspace(vol.shape[0] * 0.3, vol.shape[0] * 0.7, n).astype(int):
        sl = vol[z].astype(np.float32)
        t = filters.threshold_otsu(sl)
        hi, lo = sl[sl > t], sl[sl <= t]
        if hi.size < 50 or lo.size < 50:
            continue
        out.append(hi.mean() - lo.mean())
    return float(np.mean(out)) if out else float("nan")


def _make(voxel, fibre_d, blur, size=None, **kw):
    n = size or max(80, int(round(100.0 / voxel)))
    vol, _ = PH.build(size=n, voxel=voxel, fibre_d=fibre_d, psf_vox=blur,
                      cnr=0, seed=77, **kw)
    return observed_contrast(vol)


def fit_blur(voxels, real_contrast, fibre_d, **kw):
    """Blur per voxel size reproducing the measured contrast-transfer curve.

    A single blur value is not sufficient: a phantom retains more of its
    fine-resolution contrast than a real system does, so one voxel of blur
    under-models the real point spread at coarse sampling.
    """
    ref_real = real_contrast[min(real_contrast)]
    base = _make(min(voxels), fibre_d, 1.0, **kw)
    out = {}
    for v in voxels:
        target = real_contrast[v] / ref_real
        best = None
        for b in BLURS:
            c = _make(v, fibre_d, b, **kw)
            frac = c / base
            err = abs(np.log(max(frac, 1e-6) / max(target, 1e-6)))
            if best is None or err < best[0]:
                best = (err, b, c, frac)
        _, blur, contrast, frac = best
        out[v] = {"voxel_um": v, "blur_voxels": blur,
                  "phantom_contrast": contrast, "phantom_frac": frac,
                  "real_frac": target}
    return out


def set_noise(cal, real_cnr):
    """Noise SD from the PHANTOM's own contrast and the real true CNR.

    Using the phantom's observed contrast rather than the real one is what
    avoids double-counting the partial-volume loss.
    """
    for v, d in cal.items():
        d["real_cnr"] = real_cnr[v]
        d["noise_sd"] = d["phantom_contrast"] / max(real_cnr[v], 1e-9)
        d["achieved_cnr"] = d["phantom_contrast"] / max(d["noise_sd"], 1e-9)
    return cal


def bias_series(cal, fibre_d, misalignments=(0.5, 1.0, 2.0, 4.0, 8.0),
                size=280, measure_fn=None, **kw):
    """Measured minus true mean theta, per voxel size.

    measure_fn(vol, voxel_um, fibre_d) -> mean theta in degrees. Supply your
    own (structure tensor, or track-and-fit); there is no single right choice,
    and the bias belongs to the measurement you actually use.
    """
    if measure_fn is None:
        raise ValueError("supply measure_fn: the bias is specific to the "
                         "measurement you will apply it to")
    rows = []
    for mis in misalignments:
        for v, d in sorted(cal.items()):
            vol, truth = PH.build(size=size, voxel=v, fibre_d=fibre_d,
                                  misalign_deg=mis, psf_vox=d["blur_voxels"],
                                  cnr=d["achieved_cnr"], seed=7, **kw)
            a = truth["axes_zyx"]
            a = a.copy()
            a[a[:, 0] < 0] *= -1
            true_mean = float(np.degrees(np.arccos(np.clip(a[:, 0], -1, 1))).mean())
            meas = float(measure_fn(vol, v, fibre_d))
            rows.append({"misalign_requested": mis, "voxel_um": v,
                         "r": v / fibre_d, "true_mean": true_mean,
                         "meas_mean": meas, "bias": meas - true_mean})
    return rows


def correct(measured_by_voxel, bias_by_voxel):
    """Apply the measured bias. Report raw and corrected side by side, always."""
    return {v: {"raw": m, "corrected": m - bias_by_voxel.get(v, 0.0),
                "bias": bias_by_voxel.get(v, 0.0)}
            for v, m in measured_by_voxel.items()}


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fibre-d", type=float, required=True, help="um")
    ap.add_argument("--voxels", type=float, nargs="+", required=True)
    ap.add_argument("--real-contrast", type=float, nargs="+", required=True,
                    help="measured contrast at each voxel size (quality.py)")
    ap.add_argument("--real-cnr", type=float, nargs="+", required=True,
                    help="measured TRUE CNR at each voxel size (quality.py)")
    ap.add_argument("--vf", type=float, default=0.52)
    ap.add_argument("--out", default="calibration.json")
    a = ap.parse_args()

    if not (len(a.voxels) == len(a.real_contrast) == len(a.real_cnr)):
        raise SystemExit("--voxels, --real-contrast and --real-cnr must match "
                         "in length")

    rc = dict(zip(a.voxels, a.real_contrast))
    cn = dict(zip(a.voxels, a.real_cnr))
    print("fitting blur to the contrast-transfer curve ...")
    cal = fit_blur(a.voxels, rc, a.fibre_d, volume_fraction=a.vf)
    cal = set_noise(cal, cn)

    print("\n%8s %11s %11s %9s %10s %10s"
          % ("voxel", "real frac", "phantom", "blur", "noise SD", "CNR"))
    for v in sorted(cal):
        d = cal[v]
        print("%8.2f %11.3f %11.3f %9.1f %10.0f %10.2f"
              % (v, d["real_frac"], d["phantom_frac"], d["blur_voxels"],
                 d["noise_sd"], d["achieved_cnr"]))

    json.dump({str(k): v for k, v in cal.items()}, open(a.out, "w"), indent=1)
    print("\nwrote %s" % a.out)
    print("\nNext: call bias_series() with the measurement function you will")
    print("actually apply, then correct() the real numbers. Report raw and")
    print("corrected side by side -- never corrected alone.")
