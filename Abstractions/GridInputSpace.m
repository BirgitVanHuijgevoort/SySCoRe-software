function [uhat,InputSpace] = GridInputSpace(lu,varargin)
%GRIDINPUTSPACE grids the input space 
%
% uhat = GridInputSpace(lu,Uspace) grids the input space described by a 
% rectangular polyhedron Uspace by selecting a finite number of points (lu^dim) inside USpace.
%
% uhat = GridInputSpace(lu,Uspace,'interface',0) is equivalent to the above
%
% uhat = GridInputSpace(lu,lb,ub) grids the input space described by a 
% rectangular polyhedron with lowerbound lb and upperbound ub 
% by selecting a finite number of points inside USpace.
%
% [uhat, InputSpace] = GridInputSpace(lu,Uspace,'interface',1,Act,Fb) considers interface
% function u = uhat +  K(x-xhat), where Act \in [0,1] specifies the part of
% the input space Uspace used for the actuation and 
% Fb \in [0,1] specifies the part of the input space Uspace used for
% feedback. That is, uhat \in Act*Uspace, K(x-xhat) \in Fb*Uspace.
% The abstract input space is constructed based on the smaller input space
% based on the part for Actuation. In this case, InputSpace is supplied as
% an output of the function
%
% Inputs
% ------
% lu = number of abstract inputs in each direction
% InputSpace = input space, either as a polyhedron or by specifying the
% lower- and upperbound of the rectangular polyhedron.
%
% Outputs
% -------
% uhat = finite number of inputs
% InputSpace = continuous input space. If the interface is set to option 1
% u=uhat+K(x-xhat), then InputSpace{1} is the original space,
% InputSpace{2} is the part used for actuation and InputSpace{3} is the
% part used for feedback. 
%
% Options (varargin)
% ------
% 'interface' - Set the interface function by adding 'interface' followed by 0 or 1
% 0: u=uhat, 1: u=uhat+K(x-xhat). This outputs the new InputSpaces as a
% cell. 
% 'order' - order the abstract input space uhat by increasing absolute
% value
% 
%
% Example 
% ------
% Simple 2D input space:
% Uspace = Polyhedron(combvec([-1,1],[-1,1])');
% lu = 3;  % number of abstract inputs in each direction
% uhat = GridInputSpace(lu,Uspace); 
% 
% output equals: uhat = [-1  0  1  -1  0  1  -1  0  1; 
%                        -1 -1 -1   0  0  0   1  1  1]; 
%
% Copyright 2022 Birgit van Huijgevoort b.c.v.huijgevoort@tue.nl

%% Setting up this function based on the input
% check if input is a cell of polyhedra
disp('<---- Start finite-state abstraction')

if isa(varargin{1},'cell') 
    InputSpace = varargin{1};
    InputSpace = InputSpace{2}; % select part for actuation
    Poly = true;
elseif isa(varargin{1},'Polyhedron') % check if input is a polyhedron
    % it is a polyhedron, so compute lower and upperbound
    InputSpace = varargin{1};
    Poly = true;
else
    Poly = false;
end

for i = 1:length(varargin)
    % try to find 'interface'
    if strcmp(varargin{i},'interface')
        int_f = varargin{i+1};
        if int_f > 0 && nargin > i+1
            Act_part = varargin{i+2}; % input part for actuation
            if nargin > i+2
                Fb_part = varargin{i+3};  % Input part for feedback
            else
                Fb_part = 1-Act_part;   % default if unspecified
            end
        else % default values for actuation and feedback part (for more info see DivideInputSpace)
            Act_part = 0.75;
            Fb_part = 0.25;
        end
        break;
    else
            int_f = 0;
    end
end

order = false;
for i = 1:length(varargin)
    % try to find 'order'
    if strcmp(varargin{i},'order')
        order = true;
        break;
    end
end

%% Adjust InputSpace based on inputs of this function

% Divide the input based on these factors if the interface function is not default
if int_f > 0
    InputSpace_temp = InputSpace;
    InputSpace = cell(1,3);
    InputSpace{1} = InputSpace_temp;
    InputSpace = DivideInputSpace(InputSpace,Act_part,Fb_part); 
end

% Compute lower-bound and upperbound if it is not given.
if Poly
    if isa(InputSpace,'cell')
        InputSpace_act = InputSpace{2};
    else
        InputSpace_act = InputSpace;
    end
    dim = size(InputSpace_act.V,2);
    lb = min(InputSpace_act.V);
    ub = max(InputSpace_act.V);

else % else the lower and upperbound are given as an input
    lb = varargin{1};
    ub = varargin{2};
    dim = length(lb);
    
    if length(lb) ~= length(ub)
    error('the lower- and upperbound of the input are not the same size')
    end

end 

%% Construct abstract input space
 M = cell(dim,1);
for i = 1:dim
    M{i} = linspace(lb(i),ub(i),lu);
end

uhat = combvec(M{:});

if order
    % sort inputs based on absolute value
    [~,idx]=sort(abs(uhat));
    uhat=uhat(idx);
end

end


