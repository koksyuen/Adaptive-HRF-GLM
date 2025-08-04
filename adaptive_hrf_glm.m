function [mdl, hb_param] = adaptive_hrf_glm(freq, signal, stimulus, lambda, P_lb, P_ub, options)
% adaptive hemodynamic response function (HRF) modeling within the General Linear Model framework for activation mapping 
% in functional near-infrared spectroscopy analysis 
%
%   INPUTS:
%       freq - a scalar reflecting frequency of time series data being recorded (units: Hz)
%       signal - hemogloblin (Hb) time series data of a channel [size: 1 x number of time points]
%       stimulus - boxcar function, which is a binary time-series that equals one during task periods and zero during rest periods, used to model block desige paradigms [size: 1 x number of time points]
%       lambda - a scalar reflecting regularization coefficient that governs the weight of physiological polarity–constrained regularization within the objective function.
%       P_lb - lower bounds for six HRF parameters (m1, m2, m3, m4, c1, c2) [size: 1 x 6]
%       P_ub - upper bounds for six HRF parameters (m1, m2, m3, m4, c1, c2) [size: 1 x 6]
%       options - options for particle swarm optimization (refer to https://uk.mathworks.com/help/gads/particleswarm.html#budidgf-options)
%
%   OUTPUTS:
%       mdl - LinearModel object representing a least-squares fit of the regressors to the data (refer to https://uk.mathworks.com/help/stats/fitlm.html#bt0ck7o-mdl) 
%       hb_param - optimal HRF parameters (m1, m2, m3, m4, c1, c2) [size: 1 x 6]
%
% Author: Kok Siong Yuen
% Date: 4/8/2025
% Version: 1.0

if nargin < 7 || isempty(options)
    options = optimoptions('particleswarm','SwarmSize',700,'MaxStallIterations', 7,'UseParallel',false,'UseVectorized',true,'HybridFcn',@fmincon,'Display','iter');
end

%% particle swarm optimization
optFun = @(x)objfunc(x,freq,signal,stimulus,lambda);  % objective function definition
hb_param = particleswarm(optFun,length(P_lb),P_lb,P_ub,options);  % estimates optimal HRF parameters

%% general linear model
hrf = half_cosine_hrf(hb_param,freq); % creates HRF based on the optimal parameters
X = conv2(stimulus,hrf');  % convolves hrf with boxcar function to create a regressor
X = X(1:length(stimulus),:); 
mdl = fitlm(X, signal, 'RobustOpts', 'on');  % estimates beta-weight with robust fitting

end

