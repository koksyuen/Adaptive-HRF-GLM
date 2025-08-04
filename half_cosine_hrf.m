function y = half_cosine_hrf(P, freq)
% Modeling of hemodynamic response function (HRF) by four half-period cosine, 
% which are parameterized by a total of six parameters.
%
%   INPUTS:
%       P - six parameters (m1, m2, m3, m4, c1, c2) [size: number of particle x 6]
%       freq - frequency of time series data being recorded (units: Hz)
%
%   OUTPUTS:
%       y - four half-period cosine HRF [size: number of particle x length of hrf]
%
% Author: Kok Siong Yuen
% Date: 4/8/2025
% Version: 1.0

num_particles = size(P,1);
period = 1/freq;
t_len = max(sum(P(:,1:4),2));
t = linspace(0, t_len, t_len*freq+1);  % Time vector
hrf_len = length(t);
t_all = repmat(t,num_particles,1);

d1 = P(:,1);
d2 = d1 + P(:,2);
d3 = d2 + P(:,3);
d4 = d3 + P(:,4);

% cosine function
% y = A cos ((2pi/B)*(t - D)) + C
A = zeros(num_particles,hrf_len);
B = zeros(num_particles,hrf_len);
C = zeros(num_particles,hrf_len);
D = zeros(num_particles,hrf_len);

% amplitude, A
A1 = 0.5 * repmat(P(:,5),1,hrf_len);
A2 = -0.5 * repmat((P(:,1)>=period).* P(:,5) + 1,1,hrf_len);
A3 = 0.5 * repmat((P(:,4)>=period).* P(:,6) + 1,1,hrf_len);
A4 = -0.5 * repmat(P(:,6),1,hrf_len);

A(t_all<=d1) = A1(t_all<=d1);
A(t_all>=d1 & t_all<=d2) = A2(t_all>=d1 & t_all<=d2);
A(t_all>=d2 & t_all<=d3) = A3(t_all>=d2 & t_all<=d3);
A(t_all>=d3 & t_all<=d4) = A4(t_all>=d3 & t_all<=d4);

% period, B
B1 = 2 * repmat(P(:,1),1,hrf_len);
B2 = 2 * repmat(P(:,2),1,hrf_len);
B3 = 2 * repmat(P(:,3),1,hrf_len);
B4 = 2 * repmat(P(:,4),1,hrf_len);

B(t_all<=d1) = B1(t_all<=d1);
B(t_all>=d1 & t_all<=d2) = B2(t_all>=d1 & t_all<=d2);
B(t_all>=d2 & t_all<=d3) = B3(t_all>=d2 & t_all<=d3);
B(t_all>=d3 & t_all<=d4) = B4(t_all>=d3 & t_all<=d4);

% offset, C
C2 = 0.5 * (1 - repmat((P(:,1)>=period).* P(:,5),1,hrf_len));  % C2 = (1-f1)/2
C3 = 0.5 * (1 - repmat((P(:,4)>=period).* P(:,6),1,hrf_len));  % C3 = (1-f2)/2

C(t_all<=d1) = -A1(t_all<=d1); % C1 = -f1/2
C(t_all>=d1 & t_all<=d2) = C2(t_all>=d1 & t_all<=d2);
C(t_all>=d2 & t_all<=d3) = C3(t_all>=d2 & t_all<=d3);
C(t_all>=d3 & t_all<=d4) = A4(t_all>=d3 & t_all<=d4); % C4 = -f2/2

% shift, D
D2 = repmat((P(:,1)>=period).* P(:,1),1,hrf_len);
D3 = D2 + repmat((P(:,2)>=period).* P(:,2),1,hrf_len);
D4 = D3 + repmat((P(:,3)>=period).* P(:,3),1,hrf_len);

D(t_all>=d1 & t_all<=d2) = D2(t_all>=d1 & t_all<=d2);
D(t_all>=d2 & t_all<=d3) = D3(t_all>=d2 & t_all<=d3);
D(t_all>=d3 & t_all<=d4) = D4(t_all>=d3 & t_all<=d4);

y = A .* cos((2*pi./B) .* (t_all - D)) + C;
y(isnan(y)) = 0;

end

