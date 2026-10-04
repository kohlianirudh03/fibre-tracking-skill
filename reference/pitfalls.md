# Pitfalls

Every entry here changed a number or withdrew a claim. None was caught by
looking at output and judging it plausible; each was caught by a test with a
known answer.

## 1. A negative result from one detector is not a negative result

A carbon tow was reported as "cannot be fibre-tracked", supported by a 30-point
parameter sweep, an axial-coherence measurement and a structure-tensor test.
All consistent, all wrong. Difference-of-Gaussians scored 0.55 recall on a
phantom where Laplacian-of-Gaussian scored 1.00; with LoG, 58% of fibres
tracked in the best sub-scan.

**Cost:** a false "impossible" verdict on a usable dataset.

**Rule:** before reporting that a measurement is impossible, build a phantom
matching the real geometry, score the method on it, and confirm the method
reaches a sensible ceiling on easy cases. Report *method-limited* and
*data-limited* as separate conclusions.

## 2. A 3D watershed fragments fibres along their length

Small fluctuations in the distance transform create spurious basins. On a
0.25 µm sub-volume containing roughly 51 fibres, a 3D watershed returned 4975
objects.

**Rule:** segment each transverse slice in 2D and link centroids afterwards.
One label per fibre per slice, one track per fibre.

## 3. Otsu over the whole reconstruction separates air from specimen

That is the dominant intensity step in the volume, so fibre and matrix both
fall on the same side of the threshold and nothing is segmented.

**Rule:** build a specimen mask first, then compute Otsu on masked specimen
voxels only.

## 4. Build the specimen mask by smoothing *before* thresholding

The mask describes geometry, not texture. Thresholding raw data gives a
speckled mask that the subsequent erosion destroys — which silently emptied the
2.5 µm dataset until it was found.

**Rule:** Gaussian σ ≈ 3 voxels, then threshold, then fill holes, then erode by
the analysis footprint with `border_value=0`.

## 5. Conventional CNR counts the material's texture as noise

It read a flat 2.6–2.9 across ten reconstructions spanning a 20× range of voxel
size and two acquisition protocols. That flatness should itself have been the
warning.

Estimating noise from adjacent-slice differences — where the microstructure is
nearly identical but the noise is independent — gave true CNR of 10.4 down to
2.9, strongly resolution-dependent. It also revealed a 29.8% improvement from
thermal stabilisation at one voxel size that had read as exactly zero.

**Cost:** a conclusion about acquisition protocol, reversed.

**Rule:** `noise_sd = SD(I(z+1) - I(z)) / sqrt(2)`. `scripts/quality.py` does it.

## 6. Downsampling a fine volume does not simulate a coarse scan

Binning averages noise down as the cube root of the binning factor, so a
downsampled volume is far cleaner than a real coarse scan of the same material.

**Rule:** generate or acquire every resolution independently at its own voxel
size. Keep the voxel count fixed so the field of view scales the way a real
instrument's does.

## 7. Reproducibility is not accuracy

Two acquisitions of the same specimen agreed to 0.09°. Error against known
truth at that voxel size was 0.26°. Both scans share the same bias, so their
agreement measures precision only.

**Rule:** never offer repeat-scan agreement as evidence of accuracy. Only a
phantom separates the two.

## 8. A scale-free generator makes a scale test vacuous

A test meant to show consistent behaviour across fibre diameters compared
phantoms at D = 3 and D = 10 voxels. The two generated volumes were
bit-identical arrays: the generator is scale-free by construction, so the test
could not have failed.

**Rule:** before running a test, ask what result would falsify it. If nothing
would, it is not a test. Give the test an absolute scale the generator cannot
rescale away — a point-spread floor, a fixed detector pitch.

## 9. Calibrating a phantom is harder than it looks, in both directions

Two plausible calibrations, both wrong:

- **Supplying the real absolute noise.** Phantom contrast is fixed by
  construction while real contrast falls with coarsening, so coarse phantoms
  came out up to 13× too clean.
- **Matching CNR using the real contrast.** The measured contrast at coarse
  voxels is what survived partial-volume averaging, not a material property.
  Dividing by it double-counts the loss; coarse phantoms came out 5.3× too
  noisy.

**Rule:** fit blur per resolution to the measured contrast-transfer curve
first, then set noise from the *phantom's own* observed contrast and the real
true CNR. `scripts/calibrate.py` implements this.

## 10. Never calibrate against the quantity you are correcting

Blur is fitted to contrast transfer; noise to contrast-to-noise ratio. Neither
target is the orientation angle. Tuning a phantom against the misalignment it
is meant to correct makes the bias estimate circular and the correction
worthless.

## 11. Tune every tool before comparing, in both directions

A cross-tool comparison first ran the competitor at defaults and scored it F1
0.306; a 36-point sweep raised it to 0.644. The correction then overshot in the
mirror image — our own pipeline was run on the competitor's phantom with the
wrong fibre diameter, understating its recall as 0.324 against a properly
configured 0.631.

**Rule:** tune both tools on both generators before quoting any comparison, and
say in the write-up that you did.

## 12. Intensity registration fails on stationary texture

Two volumes of the same unidirectional specimen correlate at 0.96–0.98 *before*
any alignment, because one region looks much like another. The correlation
surface is nearly flat and has no usable peak.

**Rule:** register by matching the constellation of fibre centroids, which is
locally unique. Check transitivity (A→C should equal A→B→C). Accept that axial
offset is unrecoverable for straight fibres, and say so rather than reporting a
converged number that means nothing.

## 13. Per-fibre tables above r = 0.3 are fiction

The same code will happily emit a per-fibre table at any voxel size. Above
r ≈ 0.3 the objects it measures are fibre clusters, not fibres, and the
diameters drift away from nominal accordingly — by 3.9× at r = 1.0.

**Rule:** gate per-fibre output on `r`. Write a voxel-wise orientation field
instead, and say why.

## 14. Quality varies between sub-scans of one dataset

Across five sub-scans of one tow, tracked fraction ranged from 1% to 58%.

**Rule:** report per sub-scan. Do not generalise the best or the worst to the
dataset.
