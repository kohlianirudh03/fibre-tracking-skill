"""Ground-truth fibre phantom, and detector scoring against it.

Use this before reporting any fibre-tracking result, positive or negative. A
sweep over one method's parameters maps that method's ceiling and nothing more;
consistency across many settings of a broken method reads exactly like a real
physical limit.

  python phantom.py --fibre-d 7.0 --voxel 0.98 --spacing 8.8 --vf 0.55 \
                    --psf 2.0 --cnr 5 --score

Generate each resolution independently at its own voxel size. Downsampling a
fine volume averages noise as the cube root of the binning factor and produces
something far cleaner than a real coarse scan.
"""
import argparse
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

import fibretrack as FT


# --------------------------------------------------------------- geometry
def hex_centres(n_xy, voxel, diameter, volume_fraction, rng, margin=1.5):
    """Disordered hexagonal packing at a target area fraction.

    A real laminate or tow is a perturbed hexagonal array, not a random
    scatter: random sequential packing saturates near 54% and cannot reach
    the fractions these materials actually reach. Solve the pitch from the
    requested fraction, then jitter.
    """
    r = diameter / 2.0
    pitch = r * np.sqrt(2.0 * np.pi / (np.sqrt(3.0) * max(volume_fraction, 1e-3)))
    row_h = pitch * np.sqrt(3.0) / 2.0
    pad = margin * diameter
    extent = n_xy * voxel
    ys = np.arange(pad, extent - pad, row_h)
    pts = []
    for i, y in enumerate(ys):
        off = (pitch / 2.0) if (i % 2) else 0.0
        xs = np.arange(pad + off, extent - pad, pitch)
        for x in xs:
            pts.append((y, x))
    pts = np.array(pts, float)
    if len(pts) == 0:
        return pts
    pts += rng.normal(0.0, 0.10 * pitch, pts.shape)
    return pts


def sample_axes(n, misalign_deg, aniso, phi_deg, rng):
    """Unit axes near +z, optionally confined to a preferred plane.

    aniso 0 gives an isotropic cone about z; aniso -> 1 confines the tilt to
    the plane containing z and phi_deg. Azimuth is sampled modulo pi so the
    axes are undirected.
    """
    theta = np.abs(rng.normal(0.0, np.radians(misalign_deg), size=n))
    if aniso <= 0:
        phi = rng.uniform(0, 2 * np.pi, size=n)
    else:
        spread = max(1e-3, 1.0 - aniso) * np.pi / 2.0
        phi = (np.radians(phi_deg) + rng.normal(0.0, spread, size=n)
               + np.pi * rng.integers(0, 2, size=n))
    d = np.stack([np.cos(theta),
                  np.sin(theta) * np.cos(phi),
                  np.sin(theta) * np.sin(phi)], axis=1)      # (z, y, x)
    return d / np.linalg.norm(d, axis=1, keepdims=True)


def build(size=256, voxel=1.0, fibre_d=7.0, fibre_d_sd=0.0, volume_fraction=0.55,
          misalign_deg=1.0, aniso=0.0, phi_deg=100.0, psf_vox=1.0, cnr=8.0,
          noise_corr=1.15, fibre_level=29000, matrix_level=23000, seed=0):
    """Return (volume uint16, truth dict with exact centres and axes)."""
    rng = np.random.default_rng(seed)
    centres = hex_centres(size, voxel, fibre_d, volume_fraction, rng)
    if len(centres) == 0:
        raise ValueError("no fibres fit: check size, voxel and fibre diameter")
    axes = sample_axes(len(centres), misalign_deg, aniso, phi_deg, rng)
    radii = (fibre_d / 2.0) + rng.normal(0, fibre_d_sd / 2.0, len(centres))
    radii = np.clip(radii, 0.2 * fibre_d, None)

    mask = np.zeros((size, size, size), bool)
    yy, xx = np.mgrid[0:size, 0:size].astype(np.float32) * voxel
    zc = (size / 2.0) * voxel
    for k in range(size):
        dz = k * voxel - zc
        # each fibre's centre drifts with z according to its own axis
        cy = centres[:, 0] + axes[:, 1] / axes[:, 0] * dz
        cx = centres[:, 1] + axes[:, 2] / axes[:, 0] * dz
        sl = np.zeros((size, size), bool)
        for (y0, x0, rr) in zip(cy, cx, radii):
            if not (-rr <= y0 <= size * voxel + rr and -rr <= x0 <= size * voxel + rr):
                continue
            j0 = max(0, int((y0 - rr) / voxel)); j1 = min(size, int((y0 + rr) / voxel) + 2)
            i0 = max(0, int((x0 - rr) / voxel)); i1 = min(size, int((x0 + rr) / voxel) + 2)
            if j1 <= j0 or i1 <= i0:
                continue
            sub = ((yy[j0:j1, i0:i1] - y0) ** 2 + (xx[j0:j1, i0:i1] - x0) ** 2) <= rr ** 2
            sl[j0:j1, i0:i1] |= sub
        mask[k] = sl

    vol = np.where(mask, float(fibre_level), float(matrix_level)).astype(np.float32)
    if psf_vox > 0:
        vol = ndimage.gaussian_filter(vol, psf_vox)

    contrast = float(fibre_level - matrix_level)
    if cnr and cnr > 0:
        # correlated noise: projection noise emerges spatially correlated after
        # reconstruction, and the structure tensor responds to the noise power
        # spectrum, not only its variance. White noise understates the error.
        n = rng.normal(0.0, 1.0, vol.shape).astype(np.float32)
        if noise_corr > 0:
            n = ndimage.gaussian_filter(n, noise_corr)
            n /= (n.std() + 1e-9)
        vol = vol + n * (contrast / float(cnr))

    vol = np.clip(vol, 0, 65535).astype(np.uint16)
    truth = {"centres_yx_um": centres, "axes_zyx": axes, "radii_um": radii,
             "voxel_um": voxel, "fibre_d_um": fibre_d, "n_fibres": len(centres),
             "z_centre_um": zc, "size": size}
    return vol, truth


# ---------------------------------------------------------------- scoring
def truth_centres_at(truth, z_index):
    """Where each fibre's centre sits in slice z_index, in um."""
    dz = z_index * truth["voxel_um"] - truth["z_centre_um"]
    a = truth["axes_zyx"]
    return np.stack([truth["centres_yx_um"][:, 0] + a[:, 1] / a[:, 0] * dz,
                     truth["centres_yx_um"][:, 1] + a[:, 2] / a[:, 0] * dz], axis=1)


def score_slice(detected_yx_um, truth_yx_um, tol_um, field_um=None, edge_um=0.0):
    """Recall and precision for one slice, one-to-one nearest matching.

    Interior scoring only when edge_um > 0, so boundary truncation is not
    counted as a detection failure.
    """
    t = np.asarray(truth_yx_um, float)
    d = np.asarray(detected_yx_um, float)
    if edge_um > 0 and field_um is not None:
        keep = ((t[:, 0] > edge_um) & (t[:, 0] < field_um - edge_um) &
                (t[:, 1] > edge_um) & (t[:, 1] < field_um - edge_um))
        t = t[keep]
        if len(d):
            kd = ((d[:, 0] > edge_um) & (d[:, 0] < field_um - edge_um) &
                  (d[:, 1] > edge_um) & (d[:, 1] < field_um - edge_um))
            d = d[kd]
    if len(t) == 0 or len(d) == 0:
        return {"n_true": int(len(t)), "n_det": int(len(d)), "tp": 0,
                "recall": 0.0, "precision": 0.0}
    dist, idx = cKDTree(d).query(t, k=1)
    taken = {}
    for i, (dd, jj) in enumerate(zip(dist, idx)):
        if dd > tol_um:
            continue
        if jj in taken and taken[jj][1] <= dd:
            continue
        taken[jj] = (i, dd)
    tp = len(taken)
    return {"n_true": int(len(t)), "n_det": int(len(d)), "tp": int(tp),
            "recall": tp / len(t), "precision": tp / len(d),
            "f1": 2 * tp / (len(t) + len(d))}


def score_detectors(vol, truth, detectors=("watershed", "log"), n_slices=8,
                    tol_frac=0.5):
    """Score each detector per slice against truth. This is the whole point."""
    vox = truth["voxel_um"]
    D_vox = truth["fibre_d_um"] / vox
    field = truth["size"] * vox
    tol = tol_frac * truth["fibre_d_um"]
    zs = np.linspace(truth["size"] * 0.3, truth["size"] * 0.7, n_slices).astype(int)
    out = {}
    for name in detectors:
        fn = FT.DETECTORS[name]
        rec, prec = [], []
        for z in zs:
            pos_vox, _ = fn(vol[z], D_vox)
            pos_um = pos_vox * vox if len(pos_vox) else np.zeros((0, 2))
            s = score_slice(pos_um, truth_centres_at(truth, z), tol,
                            field_um=field, edge_um=1.5 * truth["fibre_d_um"])
            rec.append(s["recall"]); prec.append(s["precision"])
        out[name] = {"recall": float(np.mean(rec)),
                     "precision": float(np.mean(prec)),
                     "f1": float(2 * np.mean(rec) * np.mean(prec) /
                                 max(np.mean(rec) + np.mean(prec), 1e-9))}
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--size", type=int, default=192)
    ap.add_argument("--voxel", type=float, default=1.0, help="um")
    ap.add_argument("--fibre-d", type=float, default=7.0, help="um")
    ap.add_argument("--vf", type=float, default=0.55, help="volume fraction")
    ap.add_argument("--misalign-deg", type=float, default=1.0)
    ap.add_argument("--aniso", type=float, default=0.0)
    ap.add_argument("--psf", type=float, default=1.0, help="blur, voxels")
    ap.add_argument("--cnr", type=float, default=8.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--score", action="store_true",
                    help="score both detectors against truth")
    ap.add_argument("--out", help="optional .tif to write")
    a = ap.parse_args()

    g = FT.regime(a.voxel, a.fibre_d)
    print("r = %.3f (%s), fibre = %.1f voxels"
          % (g["r"], g["regime"], a.fibre_d / a.voxel))

    vol, truth = build(size=a.size, voxel=a.voxel, fibre_d=a.fibre_d,
                       volume_fraction=a.vf, misalign_deg=a.misalign_deg,
                       aniso=a.aniso, psf_vox=a.psf, cnr=a.cnr, seed=a.seed)
    print("built %d fibres in %d^3" % (truth["n_fibres"], a.size))

    if a.out:
        import tifffile
        tifffile.imwrite(a.out, vol, imagej=True,
                         metadata={"spacing": a.voxel, "unit": "um", "axes": "ZYX"})
        print("wrote", a.out)

    if a.score:
        print("\nper-slice detection against truth:")
        print("%-12s %8s %10s %7s" % ("detector", "recall", "precision", "F1"))
        for name, s in score_detectors(vol, truth).items():
            print("%-12s %8.3f %10.3f %7.3f"
                  % (name, s["recall"], s["precision"], s["f1"]))
        print("\nIf neither reaches ~1.0 on clean data (--psf 0 --cnr 0), the")
        print("detector is misconfigured, not the data unsuitable. Check that")
        print("before concluding anything about the real scan.")
