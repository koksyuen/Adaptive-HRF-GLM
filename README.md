# Adaptive-HRF-GLM

## About
The code within this repository is the reference implementation of adaptive hemodynamic response function (HRF) modeling within the General Linear Model (GLM) framework for activation mapping in functional near-infrared spectroscopy analysis described in:

A Physiological-Constrained Particle Swarm Optimization Approach for fNIRS Hemodynamic Response Function Modeling

---

## Usage

### Matlab

For Oxyhemoglobin (**HbO**) analysis:
```Matlab
[mdl, hb_param] = adaptive_hrf_glm(freq, hbo_signal, stimulus, lambda, P_lb_hbo, P_ub_hbo);
```

For Deoxyhemoglobin (**HbR**) analysis:
```Matlab
[mdl, hb_param] = adaptive_hrf_glm(freq, - hbr_signal, stimulus, lambda, P_lb_hbr, P_ub_hbr);
```

Extract GLM statistics
```Matlab
beta = mdl.Coefficients.Estimate(2);
t_value = mdl.Coefficients.tStat(2);
p_value = mdl.Coefficients.pValue(2);
```

### Inputs
- 'freq': A scalar reflecting the rate of acquisition in Hz

- 'signal': A [1 x time point] vector of hemoglobin (Hb) time-series data for a single channel

- 'stimulus': A [1 x time point] binary vector (boxcar function) indicating task (1) and rest (0) periods, used to model block-design paradigms

- 'lambda': A scalar representing the regularization coefficient governing the weight of physiological polarity–constrained regularization in the objective function

- 'P_lb': A [1 x 6] vector specifying the lower bounds for six HRF parameters *(m1, m2, m3, m4, c1, c2)*

- 'P_ub': A [1 x 6] vector specifying the upper bounds for six HRF parameters *(m1, m2, m3, m4, c1, c2)*

- 'options' *(optional)*: A structure specifying options for particle swarm optimization (see: [PSO documentation](https://uk.mathworks.com/help/gads/particleswarm.html#budidgf-options)

### Recommended Inputs values
'lambda': 0.15

For Oxyhemoglobin (**HbO**) analysis:
```Matlab
P_lb_hbo = [0 4 2 2 0 0];
P_ub_hbo = [3 8 10 8 0.1 0.5];
```

For Deoxyhemoglobin (**HbR**) analysis:
```Matlab
P_lb_hbr = [0 4 2 2 0 0];
P_ub_hbr = [4.5 14 10 12 0.25 0.5];
```

### Outputs
- 'mdl': A `LinearModel` object representing the least-squares fit of regressors to the data [fitlm documentation](https://uk.mathworks.com/help/stats/fitlm.html#bt0ck7o-mdl)

- 'hb_param': A [1 x 6] vector of the optimized HRF parameters *(m1, m2, m3, m4, c1, c2)*


