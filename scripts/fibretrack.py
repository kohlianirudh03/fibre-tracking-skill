"""Fibre detection, linking, axis fitting and cross-scan matching.

Both detectors live here. Which one you pick dominates the result, so choose it
by scoring against a phantom (phantom.py), never by looking at the output.

  detect_watershed  resolved laminates, r <= 0.3, clean fibre/matrix contrast
  detect_log        tightly packed tows, fibres nearly touching

Importable with no side effects. Angle conventions are stated in
angles_from_axis and are the ones used throughout.
"""
import argparse
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree
from skimage import feature, filters, measure, segmentation

__all__ = ["specimen_mask", "detect_watershed", "detect_log", "link_frames",
           "fit_axes", "match_axes", "angles_from_axis", "azimuth_stats",
           "orientation_tensor", "regime"]


# ------------------------------------------------------------------ regime
def regime(voxel_um, fibre_d_um):
    """Resolution ratio and what it supports. Quote r beside every angle."""
    r = float(voxel_um) / float(fibre_d_um)
    if r <= 0.3:
        name, supports = "resolved", "per-fibre orientation, diameter, tracking"
    elif r <= 3.0:
        name, supports = "marginal", "bulk orientation tensor, corrected only"
    else:
        name, supports = "subvoxel", "volume fraction and gross architecture only"
    return {"r": r, "regime": name, "supports": supports,
            "per_fibre_valid": r <= 0.3}


# -------------------------------------------------------------------- mask
def specimen_mask(sl, threshold, fibre_d_vox, smooth=3.0, ndim=2):
    """Specimen interior, eroded clear of the surface.

    Smooth BEFORE thresholding: the mask describes geometry, not texture.
    Thresholding raw data gives a speckled mask that the erosion destroys.
    """
    m = ndimage.gaussian_filter(np.asarray(sl, np.float32), smooth) > threshold
    if m.mean() < 1e-4:
        return np.zeros(m.shape, bool)
    m = ndimage.binary_fill_holes(m)
    footprint = np.ones((3,) * ndim)
    iters = int(max(4, 2 * fibre_d_vox))
    return ndimage.binary_erosion(m, footprint, iterations=iters, border_value=0)


# --------------------------------------------------------------- detectors
def detect_watershed(sl, fibre_d_vox, mask=None, min_area=None):
    """Watershed on the distance transform. For resolved laminates.

    Runs in 2D on one slice, deliberately: a 3D watershed fragments fibres
    along their length (4975 objects where ~51 fibres exist).

    Returns (centroids Nx2 in (row, col), equivalent diameters in voxels).
    """
    sl = np.asarray(sl, np.float32)
    D = float(fibre_d_vox)
    if mask is None:
        mask = np.ones(sl.shape, bool)
    if mask.sum() < 400:
        return np.zeros((0, 2)), np.zeros(0)

    sm = ndimage.gaussian_filter(sl, max(0.6, D / 13.0))
    # Otsu on specimen voxels only. Over the whole frame it finds the
    # air/specimen step and fibre and matrix land on the same side of it.
    thr = filters.threshold_otsu(sm[mask])
    fib = ndimage.binary_opening((sm > thr) & mask, np.ones((3, 3)))
    if fib.sum() < 100:
        return np.zeros((0, 2)), np.zeros(0)

    dist = ndimage.distance_transform_edt(fib)
    seeds = feature.peak_local_max(dist, min_distance=int(max(2, round(0.5 * D))),
                                   labels=fib)
    if len(seeds) < 5:
        return np.zeros((0, 2)), np.zeros(0)
    mk = np.zeros(fib.shape, np.int32)
    mk[tuple(seeds.T)] = np.arange(1, len(seeds) + 1)
    lab = segmentation.watershed(-dist, mk, mask=fib)

    amin = min_area if min_area is not None else max(9, 0.15 * np.pi * (D / 2) ** 2)
    pos, dia = [], []
    for p in measure.regionprops(lab):
        if p.area < amin:
            continue
        pos.append(p.centroid)
        dia.append(2.0 * np.sqrt(p.area / np.pi))
    if not pos:
        return np.zeros((0, 2)), np.zeros(0)
    return np.array(pos), np.array(dia)


def detect_log(sl, fibre_d_vox, mask=None, min_distance=3, thresh_sd=0.5):
    """Laplacian-of-Gaussian blob detection. For tightly packed tows.

    Fibres are bright, so the NEGATED Laplacian peaks at fibre centres, and
    sigma = half the fibre radius is the blob-matched scale.

    Scored on a 1K-tow phantom (7.1 px fibres, 1.9 px gaps): recall 1.00 /
    precision 0.95 clean, 0.89 / 0.90 at PSF 2 px and CNR 5. Difference-of-
    Gaussians on the same phantom scored 0.55, and produced a false
    "untrackable" verdict on real data.

    Returns (centroids Nx2 in (row, col), diameters in voxels — nominal, since
    LoG locates centres rather than measuring extent).
    """
    sl = np.asarray(sl, np.float32)
    r_vox = float(fibre_d_vox) / 2.0
    resp = -ndimage.gaussian_laplace(sl, sigma=0.5 * r_vox)
    if mask is not None:
        resp = np.where(mask, resp, 0.0)
    sd = resp.std()
    if sd <= 0:
        return np.zeros((0, 2)), np.zeros(0)
    peaks = feature.peak_local_max(resp, min_distance=int(min_distance),
                                   threshold_abs=thresh_sd * sd)
    if len(peaks) == 0:
        return np.zeros((0, 2)), np.zeros(0)
    return peaks.astype(float), np.full(len(peaks), float(fibre_d_vox))


DETECTORS = {"watershed": detect_watershed, "log": detect_log}


# ----------------------------------------------------------------- linking
def link_frames(frames, tol, max_gap=4):
    """Link per-slice detections into tracks, tolerating short gaps.

    frames : list of (z_um, positions Nx2 in um, diameters N in um)
    tol    : transverse link tolerance in um, typically 0.6 * fibre diameter

    The gap tolerance is not cosmetic. Without it a single bad slice retires
    the track and the fibre is lost, which was the dominant source of false
    negatives.
    """
    tracks, active, nid = {}, {}, 0
    for fi, (z_um, pos, dia) in enumerate(frames):
        if len(pos) == 0:
            continue
        live = [(i, p) for i, (p, last) in active.items() if fi - last <= max_gap]
        if not live:
            for p, d in zip(pos, dia):
                tracks[nid] = dict(z=[z_um], y=[p[0]], x=[p[1]], d=[d])
                active[nid] = (p, fi)
                nid += 1
            continue

        ids = [i for i, _ in live]
        tree = cKDTree(np.array([p for _, p in live]))
        dist, idx = tree.query(pos, k=1)
        claimed = {}
        for j, (dd, ii) in enumerate(zip(dist, idx)):
            if dd > tol:
                continue
            tid = ids[ii]
            if tid in claimed and dist[claimed[tid]] <= dd:
                continue
            claimed[tid] = j
        used = set()
        for tid, j in claimed.items():
            tracks[tid]["z"].append(z_um)
            tracks[tid]["y"].append(pos[j][0])
            tracks[tid]["x"].append(pos[j][1])
            tracks[tid]["d"].append(dia[j])
            active[tid] = (pos[j], fi)
            used.add(j)
        for j in range(len(pos)):
            if j in used:
                continue
            tracks[nid] = dict(z=[z_um], y=[pos[j][0]], x=[pos[j][1]], d=[dia[j]])
            active[nid] = (pos[j], fi)
            nid += 1
    return tracks


def fit_axes(tracks, n_slices, min_frac=0.5):
    """Total-least-squares straight axis per track.

    Suits fibres straight over the traced span. For wavy or woven material,
    fit per segment and report curvature instead.
    """
    out = []
    need = max(8, int(min_frac * n_slices))
    for tid, t in tracks.items():
        z = np.asarray(t["z"], float)
        y = np.asarray(t["y"], float)
        x = np.asarray(t["x"], float)
        if len(z) < need:
            continue
        zc = z - z.mean()
        denom = (zc ** 2).sum()
        if denom <= 0:
            continue
        sy = (zc * (y - y.mean())).sum() / denom
        sx = (zc * (x - x.mean())).sum() / denom
        v = np.array([1.0, sy, sx])
        v /= np.linalg.norm(v)
        theta, phi, phi_ax = angles_from_axis(v)
        out.append(dict(id=int(tid), n=int(len(z)),
                        y=float(y.mean()), x=float(x.mean()),
                        dia=float(np.median(t["d"])),
                        theta=theta, phi=phi, phi_axial=phi_ax,
                        axis=[float(a) for a in v],
                        span_um=float(z.max() - z.min())))
    return out


def match_axes(A, B, offset_yx, tol):
    """Pair the same fibre between two scans after a transverse offset.

    tol below half a fibre spacing (0.45 * D) makes a wrong pairing
    geometrically impossible rather than merely unlikely. Axial offset is not
    needed: for straight fibres the fitted axis is unchanged by a z shift.
    """
    if not A or not B:
        return []
    pa = np.array([[a["y"] + offset_yx[0], a["x"] + offset_yx[1]] for a in A])
    pb = np.array([[b["y"], b["x"]] for b in B])
    dist, idx = cKDTree(pb).query(pa, k=1)
    taken = {}
    for i, (d, j) in enumerate(zip(dist, idx)):
        if d > tol:
            continue
        if j in taken and taken[j][1] <= d:
            continue
        taken[j] = (i, d)
    return [dict(sep=float(d), a=A[i], b=B[j]) for j, (i, d) in taken.items()]


# ------------------------------------------------------------------ angles
def angles_from_axis(v):
    """(theta, phi, phi_axial) in degrees from a unit axis in (z, y, x).

    theta     polar misalignment from +z, 0-90; folded to the +z hemisphere
              because a fibre axis has no head or tail
    phi       azimuth atan2(dy, dx) from +x, 0-360
    phi_axial phi folded to 0-180; use THIS for circular statistics
    """
    v = np.asarray(v, float)
    if v[0] < 0:
        v = -v
    theta = float(np.degrees(np.arccos(np.clip(v[0], -1, 1))))
    phi = float(np.degrees(np.arctan2(v[1], v[2])) % 360.0)
    return theta, phi, phi % 180.0


def azimuth_stats(theta_deg, phi_axial_deg, weight_by_theta=True):
    """Mean lean plane and its concentration.

    Two corrections, both required. Doubled-angle because phi and phi+180 are
    the same axis, so an arithmetic mean is undefined. Weighted by theta
    because a near-axial fibre has a transverse component near zero, making its
    phi numerically noisy and physically uninformative.

    R is 1 for a single plane, 0 for isotropic.
    """
    th = np.asarray(theta_deg, float)
    pa = np.asarray(phi_axial_deg, float)
    w = th if weight_by_theta else np.ones_like(th)
    if w.sum() <= 0:
        return {"phi_mean_deg": float("nan"), "R": 0.0, "n": int(len(th))}
    c = np.sum(w * np.exp(2j * np.radians(pa))) / np.sum(w)
    return {"phi_mean_deg": float((np.degrees(np.angle(c)) / 2) % 180),
            "R": float(abs(c)), "n": int(len(th)),
            "weighted_by_theta": bool(weight_by_theta)}


def orientation_tensor(axes_zyx):
    """Advani-Tucker A2, its eigenvalues, FA, and the axis error.

    axis_err_deg is the angle between A2's principal direction and +z: a check
    that the nominal fibre direction is the measured one, and a fast way to
    catch a mis-stacked volume.
    """
    a = np.asarray(axes_zyx, float)
    a = a / np.linalg.norm(a, axis=1, keepdims=True)
    a[a[:, 0] < 0] *= -1
    A2 = (a.T @ a) / len(a)
    w, V = np.linalg.eigh(A2)
    order = np.argsort(w)[::-1]
    w = w[order]
    pr = V[:, order[0]]
    if pr[0] < 0:
        pr = -pr
    fa = float(np.sqrt(1.5 * np.sum((w - w.mean()) ** 2) / np.sum(w ** 2)))
    return {"A2": A2.tolist(), "A2_zz": float(A2[0, 0]),
            "eigenvalues": [float(x) for x in w], "FA": fa,
            "principal_zyx": [float(x) for x in pr],
            "axis_err_deg": float(np.degrees(np.arccos(min(1.0, abs(pr[0]))))),
            "n": int(len(a))}


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--voxel", type=float, required=True, help="voxel size, um")
    ap.add_argument("--fibre-d", type=float, required=True,
                    help="fibre diameter, um")
    a = ap.parse_args()
    g = regime(a.voxel, a.fibre_d)
    print("r = %.3f  (%s)" % (g["r"], g["regime"]))
    print("supports: %s" % g["supports"])
    print("per-fibre tables valid: %s" % ("yes" if g["per_fibre_valid"] else
                                          "NO - write a voxel-wise field instead"))
    print("\nsuggested detector: %s"
          % ("watershed" if g["per_fibre_valid"] else
             "log (or voxel-wise structure tensor only)"))
    print("score it against a phantom before trusting it (phantom.py)")
