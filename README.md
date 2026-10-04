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
