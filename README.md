# fibre-tracking

A Claude skill for tracking individual fibres through X-ray computed tomography
volumes and measuring their orientation, with the validation machinery that
decides whether the numbers mean anything.

Built from a ten-scan synchrotron resolution series on an aligned
carbon-fibre/PEEK specimen (0.25-5 um voxels, 5 um nominal fibre diameter) and
from subsequent work on twisted and hybrid yarns.

## What is here

- `SKILL.md` - when to use this, and the workflow.
- `reference/detectors.md` - watershed-on-distance-transform versus
  Laplacian-of-Gaussian, with the measured comparison.
- `reference/parameters.md` - parameter choices and what each one costs.
- `reference/pitfalls.md` - measurement traps found the hard way.
- `scripts/phantom.py` - synthetic fibre volumes with known ground truth.
- `scripts/fibretrack.py` - detection and axial linking.
- `scripts/quality.py` - scoring against a phantom; true CNR from
  adjacent-slice differences.
- `scripts/calibrate.py` - resolution-dependent bias correction.

## The one rule

**Never report a fibre measurement that has not been scored against a phantom
degraded to match the real scan.** Raw orientation error on this series spanned
23x across resolutions; after phantom-derived correction it spanned 3x. A number
quoted without that step is not a measurement.

## Status

Research code. Paths in the documentation are local to the machine it was
developed on. No warranty of any kind.

## How to cite

If this skill contributes to published work, please cite it. Machine-readable
metadata is in [`CITATION.cff`](CITATION.cff); GitHub renders a "Cite this
repository" button from it.

```
Kohli A. fibre-tracking: per-fibre tracking and orientation measurement in X-ray CT of composites. Version 1.0.0, 2026. https://github.com/kohlianirudh03/fibre-tracking-skill
```

## Method references

The methods implemented here build on:

- Emerson MJ, Jespersen KM, Dahl AB, Conradsen K, Mikkelsen LP. Individual fibre segmentation from 3D X-ray computed tomography for characterising the fibre orientation in unidirectional composite materials. *Composites Part A* 2017;97:83-92. [doi:10.1016/j.compositesa.2016.12.028](https://doi.org/10.1016/j.compositesa.2016.12.028)
- Jeppesen N, Mikkelsen LP, Dahl AB, Christensen AN, Dahl VA. Quantifying effects of manufacturing methods on fiber orientation in unidirectional composites using structure tensor analysis. *Composites Part A* 2021;149:106541. [doi:10.1016/j.compositesa.2021.106541](https://doi.org/10.1016/j.compositesa.2021.106541)
- Advani SG, Tucker III CL. The use of tensors to describe and predict fiber orientation in short fiber composites. *Journal of Rheology* 1987;31(8):751-784. [doi:10.1122/1.549945](https://doi.org/10.1122/1.549945)
- Yuksel O, Broggi G, Hartley R, et al. Benchmarking and sensitivity analysis of segmentation methods for image-based fiber detection in composites. *Composites Part A* 2026;211:110150. [doi:10.1016/j.compositesa.2026.110150](https://doi.org/10.1016/j.compositesa.2026.110150)
- Lindeberg T. Feature detection with automatic scale selection. *International Journal of Computer Vision* 1998;30:79-116.
