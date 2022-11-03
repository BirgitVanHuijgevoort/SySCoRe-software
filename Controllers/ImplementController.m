function [XSIM, QSIM] = ImplementController(x0, N, Controller, nTraj, varargin)
%IMPLEMENTCONTROLLER Simulates the controlled system initialized at x0
%
% xsim = ImplementController(x0, N, Controller); gives the state trajectory
% xsim starting at x0 and ending at x(N-1) by simulating the model connected to the Controller
%
% to obtain multiple (nTraj) state trajectories use 
% xsim = ImplementController(x0, N, Controller, nTraj)
%
% to also obtain the trajectory of the DFA state use 
% [xsim, qsim] = ImplementController(x0, N, Controller, nTraj);
%
% for model-order reduction use ImplementController(x0, N, Controller, nTraj,
% 'MOR', sysLTIr, kernel);
%
% Inputs
% ------
% x0 = initial state for simulation
% N = number of time steps (time horizon)
% Controller = struct of the controller constructed using the class RefineController
% nTraj = number of trajectories to simulate (default/unspecified = 1)
% 
% Outputs
% -------
% XSIM = states obtained during simulation. XSIM{i} corresponds to the i-th
% trajectory and is a matrix of size jxk, with j the state dimension and k
% the time horizon N. The states are given in order of time, [x(0),x(1),
% x(2), ... x(N-1)]
% QSIM = states of the DFA obtained furing simulation. Same structure as
% XSIM.
%
% Options (varargin)
% -------
% 'MOR' - let the program know that moder-order reduction has been used.
% 'MOR' should be followed by the reduced-order system and the kernel.
%
% Copyright 2022 Birgit van Huijgevoort b.c.v.huijgevoort@tue.nl

%% Initialization
disp('<---- Start deployment');

MOR = false;
% Check if MOR is given as an input
for i = 1:length(varargin)
    % try to find 'MOR'
    if strcmp(varargin{i}, 'MOR')
        MOR = true;
        if MOR
            Controller.MOR = true;
            sysLTIr = varargin{i+1};
            Controller.P = sysLTIr.P;
            Controller.Q = sysLTIr.Q;
            Controller.F = varargin{i+2};
        end
        break
    end
end

if nargin == 3
    nTraj = 1;
end

%% Start simulation
% Loop to generate nTraj unique trajectories
XSIM = cell(nTraj, 1);
QSIM = cell(nTraj, 1);
for iTraj = 1:nTraj
    if MOR
        xsim = x0; % Initial state
        % Compute reduced order initial state
        Dr = Controller.simRel.R{1}.R; % weigthing matrix for sim rel between original and reduced-order model
        xr0 = (inv(sysLTIr.P'*Dr*sysLTIr.P))*sysLTIr.P'*Dr*x0;  % Compute suitable value for xr0 (least square optimal solution of min ||x0-P*xr0||_Dr \leq \epsilon)
        xrsim = xr0; % Initial state (reduced)
    else
        xsim = x0; % Initial state
    end

    q_0 = Controller.initDFA(xsim);
    qsim = q_0; % Initial DFA state

    % Start simulating
    for i = 1:N
        % Compute next state
        if MOR
            [xrnext, qnext, u, xnext] = Controller.EvolveSys(xrsim(:, end), qsim(end), xsim(:, end));
            xrsim = [xrsim, xrnext];
        else
            [xnext, qnext, u] = Controller.EvolveSys(xsim(:, end), qsim(end));
        end
        xsim = [xsim, xnext];
        qsim = [qsim, qnext];

        % Check whether specification satisfied or violated
        if qsim(end) == Controller.DFA.F
            fprintf("Satisfaction after %d time steps.\n", i)
            break
        elseif qsim(end) == Controller.DFA.sink
            fprintf("Specification violated after %d time steps.\n", i)
            break
        end
    end

    XSIM{iTraj} = xsim;
    QSIM{iTraj} = qsim;

    % Check whether specification satisfied
    if qsim(end) ~= Controller.DFA.F
        fprintf("Simulation finished after %d time steps without satisfying.\n", i)
    end
end

% Return matrix if only one trajectory is requested
if nTraj == 1
    XSIM = XSIM{:};
    QSIM = QSIM{:};
end

disp('----> Finish deployment');
end