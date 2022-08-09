function [simRel, interface] = QuantifySim(sys, sysAbs, epsilon, mu, sigma, Bset, varargin)
% This function quantifies the similarity between a continuous-state and finite-state model  
% using simulation relation R = {(xhat,x) | ||x-xhat||_D \leq epsilon}
%
% Outputs:
% interface = interface function
% simRel = simulation relation
%
% Inputs:
% sys = original system  
% sysAbs = abstract (finite-state) system
% deviation bound = 'epsilon' or 'delta'
% value = value of deviation bound
% Bset = disturbance as an ellipsoid
%
% Options:
% % Options:
% 'interface' - specify the interface function by following 'interface'
% with 0, 1 or 2 and with an upperbound on the part of the input for K(x-xhat)
% 0. (default) u=uhat, 1. u=uhat+K(x-xhat), 2. u=R*uhat+Q*xhat+K*(x-P*xhat)
% interface function 2 is only allowed for model-order reduction (MOR).
% 'Weighting' - specify the weigthing matrix D for simulation relation R
% 
% NOT IMPLEMENTED YET
% 'distr' - 'Gaussian' if a Gaussian distribution (default) is considered 
%               'Uniform' if a uniform distribution is considered.
% Type of distribution 'Gaussian' or 'Uniform'
%
% Methods:
% NOT IMPLEMENTED YET
% 'MOR' = false/true
% If MOR is true, sysLTIr and P matrix should follow next...
% 
% NOT IMPLEMENTED YET
% 'ML' = false/true 
% If ML (multi-layered) is true, strategy should follow next...
%
% Written by: Birgit van Huijgevoort

%% Check mu and sigma

dim = sys.dim;

if sum(sum(sigma - eye(dim))) ~= 0 
   error('only sigma equal to identity is allowed')
end
if sum(sum(mu-zeros(dim))) ~= 0 
   error('only mu equal to zero allowed')
end

%% Set options

% default values if unspecified
interfaceK = 0;
givenD = 0;
uuf = zeros(size(sys.B,2),1);
distrUni = 0;

for i = 1:length(varargin)
    % try to find 'interface'
    if strcmp(varargin{i},'interface')
        interfaceK = varargin{i+1};
        if interfaceK
            uuf = varargin{i+2};
        else
            uuf = zeros(size(sys.B,2),1);
        end
    end
    % try to find weighting
    if strcmp(varargin{i},'weighting')
        D = varargin{i+1};
        givenD = 1;
    end
    % try to find distr
    if strcmp(varargin{i},'distr')
        distribution = varargin{i+1};
        if distribution == 'Uniform'
            distrUni = 1;
        else
            distrUni = 0;
        end
    end
end



%% no MOR

% if linear system
if sys.type == 'LTI'
    if givenD % if D is given
        [delta, K] = ComputeDelta2(epsilon,sys,mu,sigma,Bset,D,varargin);
    else % if D is not given
        [delta, Dmin, Kmin] = ComputeDelta(epsilon,sys,mu,sigma,Bset,varargin);
        simRel = SimRel(epsilon,delta,Dmin);
        simRel.NonDetLabels  = NonDeterministicLabelling(sysAbs.outputs, sys.regions, simRel);
        interface = Kmin;
    end
else
% if nonlinear system
disp(['A direct method for quantifying the similarity of a nonlinear system is not implemented yet. \n ...' ...
    'Please perform a piecewise-approximation first.'])
end


%% MOR
if 0
    [delta, Dmin, Kmin, Fmin] = ComputeDelta_intPQRK(epsilon,sys,sysLTIr,mu,sigma,uuf,Z,P)
end

% interface
% simRel
end

