function [sysLTI, a] = NormalizeDisturbance(sysLTI,varargin)
%TRANSFORMGAUSSIAN transforms an affine model with Gaussian disturbance
%with mean mu and variance sigma, to a system with Gaussian disturbance
%with mean 0 and variance identity
%
% Inputs:
% -------
%
% Outputs:
% --------
%
% Copyright 2022 Birgit van Huijgevoort, b.c.v.huijgevoort@tue.nl

Bw = sysLTI.Bw;
mu = sysLTI.mu;
sigma = sysLTI.sigma;

dim_dist = size(Bw,2);

a = zeros(sysLTI.dim,1);
if ~isempty(varargin)
    a = varargin{:};   
end

%% Transform system
Bw = Bw*sigma^(1/2);
a = a+Bw*mu;
mu = zeros(dim_dist,1);
sigma = eye(dim_dist);

%% Save transformed system
sysLTI.Bw = Bw;
sysLTI.mu = mu;
sysLTI.sigma = sigma;

end