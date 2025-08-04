# Adaptive-HRF-GLM

## About
The code within this repository is the reference implementation of adaptive hemodynamic response function modeling within the General Linear Model framework for activation mapping in functional near-infrared spectroscopy analysis.

## Usage
Matlab:
```Matlab
signals_corrected = TDDR(signals, sample_rate);
```


### Inputs
**signals**: A [sample x channel] matrix of uncorrected optical density data

**sample_rate**: A scalar reflecting the rate of acquisition in Hz

### Outputs
   **signals_corrected**: A [sample x channel] matrix of corrected optical density data
