# Adaptive-HRF-GLM

## About
The code within this repository is the reference implementation of adaptive hemodynamic response function modeling within the General Linear Model framework for activation mapping in functional near-infrared spectroscopy analysis.

## Usage
Matlab:
```Matlab
signals_corrected = TDDR(signals, sample_rate);
```

Python:
```Python
from TDDR import TDDR
signals_corrected = TDDR(signals, sample_rate);
```

### Inputs
**signals**: A [sample x channel] matrix of uncorrected optical density data

**sample_rate**: A scalar reflecting the rate of acquisition in Hz

### Outputs
   **signals_corrected**: A [sample x channel] matrix of corrected optical density data
#### Homer2
While Homer2 does not yet contain the TDDR method, a Homer2-compatible script is available in this repository at `toolboxes/Homer2/hmrMotionCorrectTDDR.m`. Usage is similar to other motion correction scripts shipped by Homer2.
