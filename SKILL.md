---
name: fibre-tracking
description: Track and measure individual fibres in X-ray CT or microscopy volumes of fibre-reinforced composites — per-fibre orientation (polar theta, azimuth phi), diameter, volume fraction, orientation tensor, and tracking across scans or resolutions. Use for any request about fibre orientation, fibre misalignment, fibre waviness, fibre tracking, tow twist, fibre segmentation, orientation tensor, or fibre volume fraction from tomography or micrographs; also when asked to compare scans of the same specimen, judge whether a scan can resolve fibres at all, or decide a voxel size. Covers carbon and glass fibre, UD laminates, tows and woven architectures. Triggers on "fibre tracking", "fibre orientation", "misalignment angle", "can I track fibres in this scan", "A2 tensor", "fibre volume fraction", "CT of composite", "voxel size for fibres".
---

# Fibre tracking in tomographic volumes

Measuring per-fibre orientation is easy to do and easy to do wrong. Wrong
answers here look completely plausible: a smooth number, a tight distribution,
consistent across a parameter sweep. This skill exists because several such
answers were produced and published internally before being caught.

## The one rule

**Never report a fibre-tracking result — positive or negative — that has not
been scored against a phantom with known fibre positions.**

This is not a nicety. A 30-point parameter sweep, an axial-coherence
measurement and a structure-tensor test once agreed that a carbon tow "cannot
be fibre-tracked". All three were consistent and all three were wrong: the
detector was the limit, not the data. Swapping difference-of-Gaussians for
Laplacian-of-Gaussian took per-slice recall from 0.55 to 1.00, and 58% of
fibres then tracked.

A sweep over one method's parameters maps that method's ceiling, nothing more.
Consistency across many settings of a broken method is indistinguishable from
a real physical limit. Only ground truth separates *the data cannot support
this* from *my method cannot do this*. Report those two conclusions
separately, always.

See `reference/pitfalls.md` for the full list of traps, each of which changed
a number or withdrew a claim.

## Start here: is the scan even in range?

Compute the resolution ratio before anything else:

```
r = voxel_size / fibre_diameter
```

| r | Regime | What it supports |
|---|---|---|
| ≤ 0.3 | Resolved | Per-fibre orientation, diameter, tracking |
| 0.3 – 3.0 | Marginal | Bulk orientation tensor, and only with a bias correction |
| > 3.0 | Subvoxel | Volume fraction and gross architecture only |

Quote `r` beside every angle you report. A misalignment of 2° measured at
r = 0.05 and one measured at r = 0.5 are not the same quantity, and the second
is not a measurement of the material. In the marginal regime the mean angle is
the statistic most inflated by partial-volume averaging; prefer the
distribution's shape, or correct against a calibrated phantom.

## Workflow

1. **Establish `r`** and state the regime. If r > 0.3, say so before reporting
   any per-fibre number, and expect to need `calibrate.py`.
2. **Pick a detector.** This choice dominates the result. See
   `reference/detectors.md` — the short version is Laplacian-of-Gaussian for
   tightly packed tows, watershed-on-distance-transform for resolved
   laminates. Never accept the first detector's verdict.
3. **Build a phantom** matching the real geometry — feature size, spacing,
   packing, tilt, blur, noise — with `scripts/phantom.py`, and score the
   detector against it. Check it reaches a sensible ceiling on easy cases
   before trusting it on hard ones.
4. **Detect and link** with `scripts/fibretrack.py`. Detect per slice in 2D,
   then link centroids along z. Never run a 3D watershed (see pitfalls).
5. **Fit axes and aggregate.** Angles and circular statistics have two
   corrections that are easy to miss; `reference/parameters.md` gives both.
6. **Check image quality** with `scripts/quality.py`. The conventional CNR is
   misleading on textured material — it counts the microstructure as noise.
7. **Validate.** Score against the phantom, and where possible against a second
   independent method (structure tensor vs per-fibre tracking) or a second
   acquisition. Report recall and precision, not just a count.

## Scripts

All are importable modules with no side effects on import. Run any of them with
`--help`.

| Script | Purpose |
|---|---|
| `scripts/fibretrack.py` | Detection (both detectors), linking, axis fitting, cross-scan matching, angle conventions |
| `scripts/phantom.py` | Ground-truth phantom generator and detector scoring |
| `scripts/quality.py` | True CNR from adjacent slices, contrast transfer, regime check |
| `scripts/calibrate.py` | Three-stage image-formation calibration and resolution-bias correction |

Dependencies: numpy, scipy, scikit-image, tifffile. All are in the
`fiber-tracer` venv at `C:\Users\Anirudh\.venvs\fiber-tracer`.

## Reference

- `reference/detectors.md` — choosing and scoring a detector; the LoG/DoG result
- `reference/parameters.md` — every parameter, with the reasoning, and the angle
  and circular-statistics conventions
- `reference/pitfalls.md` — the traps, each with what it cost

## Prior work to build on, not repeat

- **AMC unidirectional laminate**, ten scans at five voxel sizes: the
  resolution-bias correction and its phantom calibration live in
  `C:\Users\Anirudh\fiber_work\run`. A full methodology document is at
  `Desktop\fibre_figures\AMC_methodology.docx`.
- **CF1K carbon tow**: trackable with LoG only; quality varies strongly between
  sub-scans, so check each one rather than generalising. Scripts in
  `C:\Users\Anirudh\fibre_tracking_cf1k\scripts`.

Both are ~7 µm or 5 µm fibre diameter at roughly 1 µm voxels, so their tuned
parameters are a reasonable starting point for similar material — but re-score
on a phantom before relying on them.
