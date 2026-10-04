# Choosing a detector

The detector dominates the result. More than linking, more than thresholding,
more than the parameter sweep you are about to run. Choose it by scoring
candidates against a phantom, never by inspecting output.

## The two that work

### Laplacian-of-Gaussian — tightly packed tows

Use when fibres nearly touch: gaps of a few voxels or less, packing above
~0.5, tow or bundle geometry.

```python
resp = -ndimage.gaussian_laplace(slice, sigma=0.5 * r_vox)   # r_vox = fibre RADIUS in voxels
peaks = feature.peak_local_max(resp, min_distance=3,
                               threshold_abs=0.5 * resp.std())
```

Scored on a phantom at 1K-tow packing (7.1 px fibres, 1.9 px gaps):

| Data | Recall | Precision |
|---|---|---|
| Clean | 1.00 | 0.95 |
| PSF 2 px, CNR 5 | 0.89 | 0.90 |

The sign matters: fibres are bright, so the *negated* Laplacian peaks at fibre
centres. `sigma = 0.5 * radius` is the blob-matched scale.

### Watershed on the distance transform — resolved laminates

Use when fibres are separable after thresholding: r ≤ 0.3 and a clean
fibre/matrix contrast.

```python
sm   = ndimage.gaussian_filter(slice, max(0.6, D_vox / 13))
fib  = ndimage.binary_opening(sm > threshold_otsu(sm[mask]) , np.ones((3, 3))) & mask
dist = ndimage.distance_transform_edt(fib)
seeds = feature.peak_local_max(dist, min_distance=max(2, round(0.5 * D_vox)), labels=fib)
labels = segmentation.watershed(-dist, markers_from(seeds), mask=fib)
```

Scored against phantom truth at r = 0.05: precision 0.994, recall 0.870, mean
angular error 0.26°. Recall degrades with r — 0.83 at r = 0.10, 0.59 at
r = 0.20 — and the failures are per-slice detection failures, not linking
failures.

## Head-to-head, scored on this skill's phantom

Run `phantom.py --score` yourself; these are the numbers it gives.

| Phantom | Detector | Recall | Precision |
|---|---|---|---|
| Tow, 7 µm at 0.98 µm, clean | watershed | 0.975 | 0.996 |
| | LoG | 0.987 | 1.000 |
| Tow, PSF 2 px, CNR 5 | watershed | 0.431 | 0.779 |
| | **LoG** | **0.841** | **0.899** |
| UD, 5 µm at 1 µm (r = 0.20), PSF 1.4, CNR 10 | watershed | 0.462 | 0.744 |
| | **LoG** | **0.842** | **0.998** |

Both are fine on clean data. Under realistic blur and noise the watershed
loses roughly half its recall and LoG keeps most of its own, because the
distance transform of a blurred, noisy binary mask has far less reliable peak
structure than the LoG response of the greyscale image does.

**This bears on a published-internally result.** The AMC laminate study used
the watershed throughout and reported recall falling to 0.59 at r = 0.20,
diagnosing the misses as per-slice segmentation failures rather than linking
failures. That diagnosis looks right, and this table suggests the fix: at
comparable conditions LoG recovers 0.84 where the watershed manages 0.46. That
comparison has **not** yet been run on the real AMC data or against that
study's own calibrated phantom, so treat it as a strong reason to test rather
than as a result.

The watershed keeps one advantage LoG does not have: it measures each
cross-section's area, so it yields a diameter and an area fraction. LoG locates
centres only. Where morphometry matters, either use the watershed at
r ≤ 0.1 where its recall is intact, or use LoG for orientation and the
watershed for diameter and say which produced which.

## The one that does not

**Difference-of-Gaussians** (σ = r/√2 and r·√2) scored 0.55 recall on the same
tow phantom where LoG scored 1.00, and 0.29 at PSF 2 / CNR 5. It produced a
false "this dataset is untrackable" verdict that survived a 30-point parameter
sweep.

DoG is a fine blob detector for sparse, well-separated features. At tow packing
the outer Gaussian spans neighbouring fibres and the response between two
fibres rivals the response at their centres.

## How to choose, concretely

1. Generate a phantom at the real geometry with `scripts/phantom.py`: fibre
   diameter, centre spacing, packing fraction, tilt, PSF width, CNR — all
   measured from the real data, not guessed.
2. Score every candidate detector against it for recall and precision per
   slice. Include the easy case (clean, no blur) to confirm each reaches a
   sensible ceiling; a detector that cannot hit ~1.0 on clean data is
   misconfigured, not unsuitable.
3. Degrade the phantom to the real PSF and CNR and re-score.
4. Only then run on real data.

Scoring is per-slice detection, deliberately. Linking is not the bottleneck:
greedy and Hungarian frame-to-frame assignment scored identically (0.236 vs
0.236) on the tow data, so improving the assignment algorithm buys nothing
while the detector is the limit.

## Measuring the geometry a phantom needs

- **Fibre diameter** — from a resolved scan of the same material if you have
  one, else nominal from the fibre datasheet. Verify against the equivalent
  circular diameter of detected cross-sections once the detector works.
- **Centre spacing** — the first off-centre peak of the slice autocorrelation.
  On the 1K tow this gave 8.4–9.2 µm for ~7 µm fibres.
- **Packing fraction** — thresholded fibre area over masked specimen area.
- **PSF width** — fit blur to the measured contrast-transfer curve
  (`scripts/calibrate.py`), not by eye.
- **CNR** — from adjacent-slice differences, never the conventional
  within-phase variance (`scripts/quality.py`).

## Sanity checks on real data, once tracking runs

These catch tracks that are artefacts rather than fibres:

- **Reverse-z tracking.** Track the volume forward and backward and compare
  centres. Good data reproduced 97.5–98.6% of centres to under 1.5 px; the
  worst sub-scan managed 62–68%.
- **Stationarity.** A track whose transverse position never changes is a ring
  artefact, not a fibre.
- **Detections per slice against expectation.** Compare the count to what the
  packing fraction predicts. Every sub-scan of the tow sat at 0.87–1.07 of
  expectation, which is what made the low tracked fraction a *linking* story
  rather than a detection story.
- **A physical model.** On a twisted tow, helix angle against radius from the
  tow axis fitted a twist model at r = +0.95 over 811 fibres. An artefact
  population does not do that.
