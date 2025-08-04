function f = objfunc(param, freq, hb, stimulus, lambda)
% Objective function to estimate optimal HRF parameters by minimizing residual error
% while imposing physiological polarity–constrained regularization
%
%   INPUTS:
%       param - six parameters (m1, m2, m3, m4, c1, c2) [size: number of particle x 6]
%       freq - frequency of time series data being recorded (units: Hz)
%       hb - hemogloblin (Hb) time series data of a channel [size: 1 x number of time points]
%       stimulus - boxcar function, which is a binary time-series that equals one during task periods and zero during rest periods, used to model block desige paradigms [size: 1 x number of time points]
%       lambda - regularization coefficient that governs the weight of physiological polarity–constrained regularization within the objective function.
%
%   OUTPUTS:
%       f - objective function [size: number of particle x 1]
%
% Author: Kok Siong Yuen
% Date: 4/8/2025
% Version: 1.0

num_particle = size(param,1); % each row represents one particle
hrf = half_cosine_hrf(param,freq); % creates hemodynamic response function based on the parameters

%% General Linear Model
glm_hrf = conv2(stimulus,hrf');  % convolves hrf with boxcar function to create regressor
glm_hrf = glm_hrf(1:length(stimulus),:)';
offset = ones(num_particle,length(stimulus)); % constant baseline regressor
X = cat(3,glm_hrf,offset);  % forms design matrix

y_hbo = repmat(hb',1,1,num_particle);
design_matrix = permute(X,[2 3 1]);
beta = squeeze(pagelsqminnorm(design_matrix,y_hbo));  % estimates beta-weight via ordinary least square

X_reshaped = permute(X, [3, 2, 1]);
beta_reshaped = permute(beta,[1 3 2]);
y_pho = sum(X_reshaped .* beta_reshaped,1); 
est_y = permute(y_pho,[3 2 1]); % estimation of best-fit line
sse = sum((est_y - hb).^2, 2);  % computes residual

%% Physiological Polarity-constrained Regularization
beta_negative = beta(1,:)' < 0;
penalty = lambda * min(sse); 

%% objective function
f = sse + beta_negative .* penalty;
end