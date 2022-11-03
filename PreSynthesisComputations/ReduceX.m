function [linMod, R] = ReduceX(linMod, inputset,  goal_region, type, ...
    horizon, varargin)
% REDUCEX reduce the polytopic state space set to the set of states of 
% interest for a time-bounded invariance specification defined via a 
% set in the output space.
% 
%
% Inputs
% ------
% linMod = Object of class LinModel
% inputset = interval [u_min, u_max]
% output_region = output region of interest. In case of an invariance
% specification this is the region in which you have to stay during the given time horizon. 
% type = 'invariance'. Note: currently only invariance computations have ben
% implemented.
% horizon = Number of time steps for the specification type 
%
% Outputs 
% -------
% linMod = linMod object with state space equal to R.
% R =  the set of states from which your output can stay in the output_region for horizon
% number of times. 
% 
% Example
% ------
% sysLTI = LinModel(A,B,C,[],Bw,mu,sigma) %% see Models/LinModel for example
% [sysLTI,~] = ReduceX(sysLTI, [-1,1] , P1, 'invariance', 5);
%
% Acknowledgement: This code uses the MPT toolbox to compute the set of
% states of interest for bounded horizon computations. 
%
% Copyright 2022 Sofie Haesaert s.haesaert@tue.nl, Birgit van Huijgevoort b.c.v.huijgevoort@tue.nl

%% Preliminaries

% Check whether conditions under which this code holds are satisfied
if ~strcmp(type, 'invariance')
    error(['The type = ', type, ' is not yet known or implemented'])

end

if isa(inputset, 'Polyhedron')
    % Compute lower- and upperbound of input set if not given
    ula = min(inputset.V,[], 'all');    % lower-bound on input space for actuation
    uua = max(inputset.V,[], 'all');    % upper-bound on input space for actuation

    inputset = [ula,uua];
end

% Compute output_region by enlarging the goal region
[~, output_region] = IncreaseDecreasePolytope(goal_region, 0.1);

% Compute initial set
R_init = Polyhedron(output_region.A*linMod.C, output_region.b);

% load the linMod class into the LTISystem object of the mpt toolbox
system = LTISystem('A',  linMod.A, 'B', linMod.B);
system.x.min = min(linMod.X.V, [], 1); % compute a box around the given state space of the linMod system
system.x.max = max(linMod.X.V, [], 1); % compute a box around the given state space of the linMod system
system.u.min = inputset(:,1); % compute a box around the given input space of the linMod system
system.u.max = inputset(:,2); % compute a box around the given input space of the linMod system

if size(inputset,1) ~=1
    warning('works for 1 inputs only')
end
U = Polyhedron(inputset');

% %% Debugging figures
% figure(1)
% subplot(2,1,1)
% hold off
% R_init.plot()
% title(['backward invariance set i = ', num2str(1)])

%% Init reach set
R= cell(horizon,1);

R{1}= R_init;
V = [ ];
for i=1:horizon
    R_ = system.reachableSet('X', R{i}, 'U', U, 'N', 1, ...
        'direction', 'backward');
    R{i+1} = R_ & R_init;
    V = [V;R{i+1}.V];
end

%% Debugging figures
%debugfigs(i, R)

%% Retrieve reduced model
X = Polyhedron(V);
X.minVRep;
linMod.X = X;
end


function debugfigs(i, R)
% Debugging figures
figure(1)
tl = tiledlayout(2, i+1);

% Print backward invariant sets
for j = 1:i+1
    ax = nexttile;
    R{j}.plot();
    title(ax, ['Backward invariance, i = ', num2str(j)])
    subtitle(ax, 'Invariant set')
    xl{j} = xlim(ax);
    yl{j} = ylim(ax);
end

% Print vertices of backward invariant sets
for j = 1:i+1
    ax = nexttile;
    scatter(R{j}.V(:,1), R{j}.V(:,2), 60, [204 204 255]/256, 'filled')
    subtitle(ax, 'Vertices')
    xlim(ax, xl{j})
    ylim(ax, yl{j})
    grid on
end
end