function [sysAbs] = FSabstraction(varargin)
%FSABSTRACTION computes a finite-state abstraction 
%
% sysAbs = FSabstraction(sys,uhat,l,tol,DFA) computes a finite-state 
% abstraction by gridding the state space of dynamic systems sys
%
% Including efficient tensor computation for 2D systems
% sysAbs = FSabstraction(sys,uhat,l,tol,DFA,'TensorComputation',true);
% 
% Inputs
% ------
% sys   = systems of class LinModel, PWAmodel or NonlinModel (see folder Models)
% uhat  = Finite set of inputs,  example uhat = [-1, 0, 1];
% l     = Number of finite states in each dimension  [l1 l2 l3 ...]
% tol   =  tolerance for truncating to 0
% DFA   = deterministic finite-state automation, (see function TranslateSpec on what is contained in this struct) 
%  
% Output struct sysAbs (= finite-state abstraction of sys) consisting of:
% -------
% states
% outputs corresponding to states in sysAbs.states
% inputs = finite number of inputs
% dim = dimension
% orig = original continuous-state system
% P = probability matrix, often as a TensorTransitionProbability (see
% folder Models)
% labels = labels corresponding to states in sysAbs.states
% beta = disturbance caused by finite-state abstraction (vector that pushes
% states to centers)
% hx = outputs corresponding to states in sysAbs.states (equal to sysAbs.outpus if Bw = I)
% zstates = states after transformation when we do not have Bw = I. 
% l = [l1, l2, l3 ...] = number of finite-states in each dimension
% Partition = (for nonlinear/PWA systems) partition corresponding to states
% in sysAbs.states
% outputmap, for LTI systems with output y = Cx, outputmap = C-matrix (when Bw = I)
%
% Options
% -------
% 'TensorComputation' = false/true 
% 'TensorToolbox' = {'2d'(default), 'tensorlab', 'tensortoolbox' }
%
% Copyright 2021 Sofie Haesaert s.haesaert@tue.nl
% 2022 Birgit van Huijgevoort b.c.v.huijgevoort@tue.nl
 
sys = varargin{1};
uhat = varargin{2};
l = varargin{3};
tol = varargin{4};
DFA = varargin{5};

% Default if TensorComputation is not specified
TensorComputation = 0;

for i = 6:length(varargin)
    % try to find 'TensorComputation'
    if strcmp(varargin{i},'TensorComputation')
        TensorComputation = logical(varargin{i+1});
    end
end

if sys.dim == 1
    if TensorComputation
        warning('Tensor computation not implemented for 1D systems, abstract model computed without tensors instead.')
        TensorComputation = 0;
    end
end
    
% 'TensorComputation' is False
if ~TensorComputation 
    [sysAbs] = GridSpace_nd(sys, uhat, l, tol);
    sysAbs.l = l;
end

% 'TensorComputation' is True
if TensorComputation
    TensorToolbox = '2d';
    for i = 6:length(varargin)
        % try to find 'TensorToolbox'
        if strcmp(varargin{i},'TensorToolbox')
            TensorToolbox = varargin{i+1};
        end

    end
    
    if sys.dim > 2 && strcmp(TensorToolbox,'2d')
        error("System has more than 2 dimensions. Specify which tensor toolbox to use. (TensorToolbox = {2d, tensorlab, tensortoolbox})");
    end

    [sysAbs] = GridSpace_nonlin_tensor(sys, uhat,l, tol, TensorToolbox);
end

% Save some extra system parameters into struct
sysAbs.orig = sys;

% Label output space
sysAbs = DeterministicLabelling(sysAbs, DFA, sys.regions, sys.AP);

if sys.type == 'PWA'
    % Determine for each abstract state in which partition it lies
    %%% Assumes rectangular grid cells
    temp = ceil((sysAbs.states-(min(sys.X.V)'))./((max(sys.X.V)-min(sys.X.V))').*(sys.N(1:sys.dim)'-1)); %Coordinate of the partition ([row;column] in 2D)
    sz = sys.N(1:sys.dim)-1;
    M = [];
    for i = 1:size(temp,1)
        M = [M {temp(i,:)}];
    end
    sysAbs.Partition = sub2ind(sz, M{:});
end

disp('----> Finish finite-state abstraction')
end

