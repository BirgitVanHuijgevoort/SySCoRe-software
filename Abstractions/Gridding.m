function [sysAbs] = Gridding(varargin)
%GRIDDING This function grids the state space and input space of dynamic
%systems
% sysAbs = Gridding(sys, uhat, l, tol)
% 
% Inputs: 
% sys   = LTI systems with fields A and B 
% uhat  = Finite set of inputs,  example uhat = combvec(linspace(ul(1),uu(1),3),linspace(ul(2),uu(2),3));
% l     = Number of finite states in each dimension  [l1 l2 l3 ...]
% tol  =  tolerance for truncating to 0
%  
% Options:
% 'TensorComputation' = false/true 
% 'TensorToolbox' = {'2d'(default), 'tensorlab', 'tensortoolbox' }
 
sys = varargin{1};
uhat = varargin{2};
l = varargin{3};
tol = varargin{4};

% Default if TensorComputation is not specified
TensorComputation = 0;

for i = 5:length(varargin)
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
end

% 'TensorComputation' is True
if TensorComputation
    TensorToolbox = '2d';
    for i = 5:length(varargin)
        % try to find 'TensorComputation'
        if strcmp(varargin{i},'TensorToolbox')
            TensorToolbox = varargin{i+1};
        end

    end

    [sysAbs] = GridSpace_nonlin_tensor(sys, uhat,l, tol,TensorToolbox);
end
end

