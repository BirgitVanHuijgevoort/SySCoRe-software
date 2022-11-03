function [simRel, interface, varargout] = QuantifySim(sys, sysAbs, epsilon, varargin)
% QUANTIFYSIM Quantify the similarity between a continuous-state (sys) and abstract (finite-state or reduced-order) model (sysAbs)  
% using simulation relation R = {(xhat,x) | ||x-xhat||_D \leq epsilon}
%
% Inputs:
% -------
% sys = original system  
% sysAbs = abstract system
% epsilon = output deviation ||y-yhat|| \leq \epsilon
%
% Outputs:
% -------
% simRel = simulation relation
% interface = interface function. For PWA systems interface equals sys with
% added field Kf for the interface function u=uhat+Kf(x-xhat)
%
% Options (= varargin)
% --------
% 'interface' - specify the interface function by following 'interface'
% with 0, or 1. Here 0. (default) u=uhat, 1. u=uhat+K(x-xhat)
% Example: 
% QuantifySim(sys, sysAbs, 0.3, 'interface', 1)
% see Tutorials/VanderPol for full example
%
% 'weighting' - specify the states used to compute weigthing matrix D for simulation relation R
% Example: 
% D = [1 0; 0 1];
% QuantifySim(sys, sysAbs, 0.3, 'weighting', D)
% see Tutorials/VanderPol for full example
%
% 'MOR' - Quantify the similarity between a full-order and reduced-order model.
% Input sysAbs should be the reduced-order model, and the finite-state model should follow after 'MOR'. 
% Interface function is automatically set to u=uhat+Q*xhat+K*(x-P*xhat). 
% Output interface is K-matrix, Output varargout is F-matrix, 
% with F used to compute the reduced disturbance as w_r = w+F(x-Px_r).
% Example: QuantifySim(sys, sysLTIr, 0.3, 'MOR', sysAbs)
% see Tutorials/BAS for full example
%
% 'fast' - influence the trade-off between computation time and accuracy.
% 'fast' is followed by a number between 0 and 1, here 0 is a fast but
% inaccurate computation, 1 is a slow but accurate computation (of delta)
% Example: QuantifySim(sysLTI, sysAbs, epsilon, 'interface', int_f, 'fast', 0.5);
%
% 'distr' - 'Gaussian' if a Gaussian distribution (default) is considered 
%               'Uniform' if a uniform distribution is considered.
% Type of distribution 'Gaussian' or 'Uniform'. Warning, not fully tested!
%
% This function is based on the method described in:
% van Huijgevoort, B. C., & Haesaert, S. (2022). Similarity quantification 
% for linear stochastic systems: A coupling compensator approach. Automatica. 
% 
% Copyright 2022 Birgit van Huijgevoort b.c.v.huijgevoort@tue.nl

%% Check mu and sigma
disp('<---- Start similarity quantification')

dim = sys.dim;
mu = sys.mu;
sigma = sys.sigma;

if all(diag(sys.sigma)==1) ~= 1
   error('only sigma equal to identity is allowed')
end
if sum(mu) ~= 0 
   error('only mu equal to zero allowed')
end

%% Set options

% default values if unspecified
interfaceK = 0;
givenD = 0;
MOR = 0;
if sys.type == 'PWA'
    uuf = zeros(size(sys.Partition(1).Dynamics.B,2),1);
else
    uuf = zeros(size(sys.B,2),1);
end
distrUni = 0;

for i = 1:length(varargin)
    % try to find 'interface'
    if strcmp(varargin{i},'interface')
        interfaceK = varargin{i+1};
        if interfaceK
            uuf = max(sys.U{3}.V,[],'all');
        else
            uuf = zeros(size(sys.B,2),1);
        end
    end
end
for i = 1:length(varargin)
    % try to find weighting
    if strcmp(varargin{i},'weighting')
        D = varargin{i+1};
        givenD = 1;
    end
end
for i = 1:length(varargin)
    % try to find MOR
    if strcmp(varargin{i},'MOR')
        MOR = 'true';
        sysLTIr = sysAbs; 
        sysAbs = varargin{i+1};
        P = sysLTIr.P;

        uuf = max(sysLTIr.U{3}.V);
    end
end
for i = 1:length(varargin)
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
if ~MOR

    % if linear system
    if sys.type == 'LTI'
        Bset = sysAbs.beta;
        if givenD % if D is given
            [delta, K] = ComputeDelta2(epsilon,sys,mu,sigma,Bset,D,varargin);
        else % if D is not given
            [delta, Dmin, Kmin] = ComputeDelta(epsilon,sys,mu,sigma,Bset,varargin);
            simRel = SimRel(epsilon,delta,Dmin);
            if sys.dim == 1 || size(sys.C,1) == 1
                simRel.NonDetLabels  = NonDeterministicLabelling(sysAbs.outputs, sys.regions, simRel);
            else
                simRel.NonDetLabels  = NonDeterministicLabelling(sysAbs.outputs, sys.regions, simRel, 'Efficient', sysAbs);
            end
            interface = Kmin;
        end
    elseif sys.type == 'PWA'
        if givenD
            
            % Rewrite (struct --> matrices) to be able to use parallel computing
            for i = 1:sys.Np
                Kset(i) = sys.Partition(i).K;
                
                Bset = plus(sys.Partition(i).K,sysAbs.beta);
                Bset = minVRep(Bset);
                sys.Partition(i).Bset = Bset;
                Bset_all(i) = sys.Partition(i).Bset;
                
                Dynamics(i) = sys.Partition(i).Dynamics;
            end
        
            % Compute delta for each partition
            f_delta = ones(1,sys.Np);
            KfK = [];
            parfor i = 1:sys.Np
                Bset = plus(Kset(i),sysAbs.beta);
                Bset = minVRep(Bset);
                
                [delta, Kf] = ComputeDelta2(epsilon,Dynamics(i),mu,sigma,Bset,D,'interface',1,uuf);
            
                f_delta(1,i) = delta;
                KfK = [KfK; Kf];
            end
            for i = 1:sys.Np
                sys.Partition(i).delta = f_delta(1,i);
                sys.Partition(i).Kf = KfK(i,:);
                sys.Partition(i).rel = SimRel(epsilon,sys.Partition(i).delta,eye(2));
            end
            
            %disp([', epsilon = ', num2str(epsilon), 'delta = ', num2str(f_delta)])
        
            % Define simulation relation
            simRel = SimRel(epsilon,f_delta,eye(2));
            
            % Determine labelling 
            simRel.NonDetLabels  = NonDeterministicLabelling(sysAbs.outputs, sys.regions, simRel);

            % Get delta in correct shape for SynthesizeRobustController
            simRel.delta = simRel.delta(sysAbs.Partition)';
    
            interface = sys;
        else
            error(['The similarity quantification for nonlinear systems requires a weighting matrix to be given. \n ...' ...
                'Please compute a weighting matrix before quantifying the similarity.'])
        end
    elseif sys.type == 'NonLin'
    % if nonlinear system
        error(['A direct method for quantifying the similarity of a nonlinear system is not implemented yet. \n ...' ...
            'Please perform a piecewise-approximation first.'])
    end
end


%% MOR
if MOR
    R = 1;
    % Compute polytope 
    beta = sysAbs.beta;
    Uhat = Polyhedron(sysAbs.inputs');
    
    Wlb = sysLTIr.mu-3*sum(sysLTIr.sigma,2);
    Wub = sysLTIr.mu+3*sum(sysLTIr.sigma,2);
    Wset = Polyhedron('lb', Wlb, 'ub', Wub);
    
    % Compute additional error, by truncating the disturbance
    onemindel = mvncdf(Wlb,Wub,mu,sigma);
    del_trunc = 1-onemindel;
    
    Z = (sys.B*R-P*sysLTIr.B)*Uhat+(sys.Bw-P*sysLTIr.Bw)*Wset;
    Zred = Z;
    Zred = Z.minVRep(); 
    
    % Compute MOR simulation relation
    [delta, D, K, F] = ComputeDelta_intPQRK(epsilon,sys,sysLTIr,mu,sigma,uuf,Zred,P);
    delta = delta+del_trunc;

    interface = K;
    varargout{1} = F;
    simRel = SimRel(epsilon,delta,D);
end

% interface
% simRel

disp('----> Finish similarity quantification')
end

