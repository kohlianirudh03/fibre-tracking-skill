# Parameters and conventions

Values validated on 5 µm carbon fibre in a UD laminate at r = 0.05–1.00, and on
~7 µm fibre in a 1K tow at r ≈ 0.14. Re-score on a phantom before relying on
them for a different material.

`D_vox` is the fibre **diameter** in voxels; `r_vox` the **radius**.

## Angle conventions

Volumes are indexed `(z, y, x)`, z along the slice stack and nominally along
the fibres. A fibre axis is a unit vector folded to the +z hemisphere, because
an axis has no head or tail.

```
theta     = arccos(|v_z|)                 polar misalignment, 0-90 deg
phi       = atan2(v_y, v_x) mod 360       azimuth from +x in the x-y plane
phi_axial = phi mod 180                   folded, for circular statistics
```

θ is the quantity that enters kink-band strength models. φ is meaningful only
in aggregate.

## Circular statistics for the azimuth

Two corrections, both required. Omitting either gives a meaningless number.

**Fold the angle.** φ and φ+180° describe the same axis, so an arithmetic mean
of φ is undefined. Use the doubled-angle transform.

**Weight by θ.** A nearly axial fibre has a transverse component near zero, so
its φ is numerically noisy and physically uninformative. Weighting by θ
suppresses that without discarding data.

```
c       = sum(theta_i * exp(2j * radians(phi_axial_i))) / sum(theta_i)
phi_bar = (degrees(angle(c)) / 2) mod 180      mean azimuth of the lean plane
R       = abs(c)                               1 = single plane, 0 = isotropic
```

Two variants of this exist in past outputs — θ-weighted and uniform-weighted —
and they differ. State which you used; do not mix them in one comparison.

## Orientation tensor

```
A2 = (1/N) * sum(outer(v_i, v_i))
FA = sqrt(1.5 * sum((lam - lam.mean())**2) / sum(lam**2))
```

`A2[0,0]` is 1 for perfect alignment with z. The angle between A2's principal
eigenvector and +z is the *axis error* — a check that the nominal fibre
direction is the measured one, and a fast way to catch a mis-stacked volume.

## Detection

| Parameter | Value | Why |
|---|---|---|
| Specimen mask smoothing | Gaussian σ = 3 vox | mask is geometry, not texture |
| Mask erosion (2D) | 3×3, `max(4, 2*D_vox)` iterations | keep the surface out |
| Mask erosion (3D) | 3×3×3, halo iterations, `border_value=0` | volume edge is outside |
| Denoise | Gaussian σ = `max(0.6, D_vox/13)` | tied to fibre size, so equivalent across resolutions |
| Threshold | Otsu on **masked specimen voxels only** | else it finds the air/specimen step |
| Opening | 3×3 | removes single-voxel speckle |
| Seed separation (watershed) | `max(2, round(0.5*D_vox))` vox | half a fibre |
| LoG scale | σ = `0.5 * r_vox` | blob-matched |
| LoG peak separation | 3 vox | at tow packing |
| LoG threshold | `0.5 * response.std()` | |
| Minimum region area | `max(9, 0.15*pi*(D_vox/2)**2)` vox | rejects watershed slivers |
| Border regions | excluded for morphometry | a truncated section under-reports its diameter |

## Linking and fitting

| Parameter | Value | Why |
|---|---|---|
| Slice spacing | 1.0 µm, resampled to whole voxels | same physical length at every resolution |
| Axial span | 180 µm (220 µm for export) | |
| Link tolerance | `0.6 * D` = 3.0 µm at D = 5 µm | |
| Maximum gap | 4 sampled slices | without it one bad slice ends the track — this was the dominant false-negative source |
| Minimum track length | `max(8, min_frac * n_slices)` | `min_frac` 0.25–0.5 |
| Cross-scan match tolerance | `0.45 * D` = 2.25 µm | below half a fibre spacing, so a wrong pairing is geometrically impossible |

Axis fitting is total least squares on the centroid cloud:

```
s_y = sum((z - z.mean()) * (y - y.mean())) / sum((z - z.mean())**2)
s_x = sum((z - z.mean()) * (x - x.mean())) / sum((z - z.mean())**2)
v   = normalise((1, s_y, s_x))
```

The straight-axis model suits fibres that are straight over the traced span.
For wavy or woven architecture, fit per-segment and report curvature; the
phantom generator can produce both to test the limit.

Reported fibre diameter is the **median** of per-slice equivalent diameters,
which is robust to the occasional merged cross-section.

## Structure tensor (voxel-wise, independent of tracking)

| Parameter | Value |
|---|---|
| Derivative scale σ | 1 voxel |
| Integration scale ρ | `max(D/2, 2 voxels)` |
| Halo | `ceil(4*rho/voxel) + 4` voxels |
| Block core | 128 voxels |
| Minimum masked voxels per block | 64 |
| θ histogram | 900 bins over 0–90° |
| Pole figure | 0.5° × 5° |

Tying ρ to the fibre radius keeps the operator measuring the same physical
structure at every resolution; the two-voxel floor stops it collapsing below
the sampling limit at coarse scans.

Erode the mask by the same footprint as the halo. The specimen surface produces
the strongest gradients in the volume and will otherwise bias orientations
toward the surface normal.

Block, and size workers from the **double-precision internal array**, not from
the input dtype — some implementations allocate float64 regardless, which is
~1.9 GB per worker at 224³ and ~14.5 GB unblocked at 600³.

Use this as a second, independent measurement. It fails differently from
tracking, so agreement between the two is informative; a tracking result with
no independent corroboration is weaker than it looks.

## Phantom defaults (UD laminate)

| Parameter | Value |
|---|---|
| Volume | 280³ voxels, same count at every resolution |
| Fibre diameter | 5.0 µm, SD 0.35 µm |
| Volume fraction | 0.52 |
| Packing | disordered hexagonal, pitch `(D/2)*sqrt(2*pi/(sqrt(3)*Vf))` |
| Misalignment plane | φ = 100°, anisotropy 0.72 |
| Noise correlation length | 1.15 voxels |

Random sequential packing saturates near 54% area fraction and cannot reach
what these materials reach; a real laminate is a perturbed hexagonal array, so
solve the lattice pitch from the target fraction and jitter the centres.

Noise must be **spatially correlated**, not white: projection noise propagates
through reconstruction and emerges correlated, and the structure tensor
responds to the noise power spectrum, not only its variance. White noise
understates the orientation error.

## Image quality

```
contrast  = mean(fibre voxels) - mean(matrix voxels)      within the specimen mask
noise_sd  = SD(I(z+1) - I(z)) / sqrt(2)                   adjacent-slice difference
CNR_true  = contrast / noise_sd
```

At coarse voxels `contrast` is not a material property — it is what survived
partial-volume averaging. That distinction is what the two failed calibrations
in `pitfalls.md` got wrong.
