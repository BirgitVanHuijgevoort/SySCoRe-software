function [satProb_lb, pol, varargout] = SynthesizeRobustController( ...
    sysAbs, DFA, rel, thold, varargin)
%SynthesizeRobustController synthesizes a control policy and computes
% the corresponding robust satisfaction probability
%
% Use [satProb_lb, pol] = SynthesizeRobustController(sysAbs, DFA, rel) to
% obtain a satisfaction probablity (lowerbound) and policy. 
%
% To change the threshold for convergence use for example 
% thold = 1e-6;
% [satProb_lb, pol] = SynthesizeRobustController(sysAbs, DFA, rel, thold)
%
% To only compute the satisfaction probability for the initial DFA state
% use [satProb_lb, pol] = SynthesizeRobustController(sysAbs, DFA, rel, [], true)
%
% Inputs
% ------
% sysAbs = abstract finite-state system
% DFA = deterministic finite-state automaton
% rel = simulation relation
% N = maximum number of iterations for value iteration
% 
% Outputs
% --------
% satProp = robust satisfaction probability = lowerbound on actual
% satisfaction probability. 
% satProp(i,j) is satisfaction probability for DFA state i and initial
% abstract state sysAbs.states(:,j).
% pol = abstract control policy.
% pol(:,i,j) is optimal abstract input for state (xhat,q) with
% xhat = sysAbs.states(:,i) and q = DFA.S(j)
%
% Options (varargin)
% -------
% 'initialonly' = true/false. Default is false. If it is set to true, this
% function only computes the value function and satisfaction probability for the initial DFA state
% 'antagonist_pol' = true/false. Default is false.  This option will explicitly compute the
% policy with which the worst-case labels are chosen [Warning not fully
% tested]. This will lead to an extra output:
% [satProp,pol,a_pol] = SynthesizeRobustController(sysAbs,DFA,rel, N, 0, 1)
% 'upperBound' = true/false. Default is false.  This option will compute the
% upper bound on the satisfaction probability in addition to the lower bound 
% [Warning not fully tested]. This will lead to an extra output:
% [satProp,pol,a_pol,satProb_ub] = SynthesizeRobustController(sysAbs,DFA,rel, N, 0, 1, 1) or
% [satProp,pol,satProb_ub] = SynthesizeRobustController(sysAbs,DFA,rel, N, 0, 0, 1)
%
% Copyright 2022 Sofie Haesaert s.haesaert@tue.nl, Birgit van Huijgevoort
% b.c.v.huijgevoort@tue.nl, Oliver Schoen o.schoen2@newcastle.ac.uk
%
% References:
% Haesaert, Sofie, and Sadegh Soudjani. "Robust dynamic programming for temporal
% logic control of stochastic systems." IEEE Transactions on Automatic Control (2020).

%% Init
checkDFAact(DFA) % Check validity of DFA.act

disp('<---- Start synthesizing a robust controller');

% Set default value for thold if it is not given
if nargin<4
    thold = 1e-12;
end
if isempty(thold)
    thold = 1e-12;
end

if nargin>=5
    if islogical(varargin{1})
        initialonly = varargin{1};
    else
        error(['5th input argument should be true/false for ' ...
            'the output of the value function and policy ' ...
            'for initial (discrete) state only'])
    end
else
    initialonly = false;
end

if nargin>=6
    if islogical(varargin{2})
        antagonist_pol = varargin{2};
    else
        error(['6th input argument should be true/false for ' ...
            'the computation of the counter strategy of the labeling'])
    end
else
    antagonist_pol = false;
end

if nargin>=7
    if islogical(varargin{3})
        upperBound = varargin{3};
    else
        error(['7th input argument should be true/false for ' ...
            'the computation of the upper bound on the satisfaction probability'])
    end
else
    upperBound = false;
end

%% Perform value iteration
N = 500; %maximum number of iterations performed.  
uhat = sysAbs.inputs;
delta = rel.delta;
nS = length(DFA.S);
nX = length(sysAbs.states);

if ~isempty(rel.NonDetLabels)
    outputs2act = rel.NonDetLabels;
else
    error("The simulation relation does not have non-deterministic labels. Run the NonDeterministicLabelling function to create one.");
end

% Initialise value function
% V(i,j) is the probability of reaching F from DFA state i and abstract state
% sysAbs.states(:,j)
V_lb = zeros(nS,nX);
if upperBound
    V_ub = zeros(nS,nX);
end
DFA_Active = setdiff(setdiff(DFA.S,DFA.F), DFA.sink);
Converged = ones(1, nS); % create a vector with
Converged(DFA_Active)=deal(0);

V_lb(DFA.F,:) = ones(nX,1); % Set V to 1 for DFA state q = F
if upperBound
    V_ub(DFA.F,:) = ones(nX,1); % Set V to 1 for DFA state q = F
end
tic

% Prepare DFA transitions for all states states
trans = ones(nS,nS,nX);
for i = DFA_Active
    for l = 1: size(outputs2act,1)
        trans(i, DFA.trans(i, l), :) = min( ...
            shiftdim(trans(i, DFA.trans(i, l), :), 1), ...
            1 - outputs2act(l, :));
    end
end
trans = 10 * trans;


for k = 1:N
    % Stop iterating when all values are converged
    if min(Converged) == 1
        disp('Convergence reached!')
        disp(['Number of iteration required, k=', num2str(k)])
        break;
    end

    for i = DFA_Active(end:-1:1) % for each discrete mode
        % check if mode has converged
        if min([Converged(DFA.trans(i, :)), Converged(i)]) == 1
            continue
        else
            Converged(i)= 0;

            % Choose the correct states out of V based on the DFA
            V_sort = min(max(squeeze(trans(i, :, :)), V_lb));

            % Compute probability of stepping to sink state for upper bound
            P_sink = 1 - ones(1, nX) * sysAbs.P;

            % Compute value function for each action uhat
            VVlb = V_sort * sysAbs.P; % P is an object of the TensorTransitionProbability class
            if upperBound
                VVub = VVlb + P_sink;
            end

            % Optimize over uhat and subtract delta
            V_n_lb = max(VVlb, [], 2) - delta;
            if upperBound
                V_n_ub = max(VVub, [], 2) + delta;
            end

            % Make sure that value function is between 0 and 1
            V_n_lb = min(1, max(0, V_n_lb));
            if upperBound
                V_n_ub = min(1, max(0, V_n_ub));
            end

            if min([Converged(DFA.trans(i, :))]) == 1
                Converged(i) = 1;
            elseif max(max(abs(V_n_lb' - V_lb(i, :)))) < thold
                Converged(i) = 1;
            end
            V_lb(i, :) = V_n_lb;
            if upperBound
                V_ub(i, :) = V_n_ub;
            end

        end
    end
end


%% Compute satisfaction probability
if initialonly == true
    % Determine correct q_0 for abstract states
    satProb_lb = V_lb((0:length(V_lb)-1)*nS + DFA.trans(DFA.S0, sysAbs.labels));
    if upperBound
        satProb_ub = V_ub((0:length(V_ub)-1)*nS + DFA.trans(DFA.S0, sysAbs.labels));
    end
else
    satProb_lb = V_lb;
    if upperBound
        satProb_ub = V_ub;
    end
end

%% Compute optimal policy
pol = zeros(size(uhat, 1), nX, nS);
a_pol = zeros(size(uhat, 1), nX, nS);

for i = DFA_Active
    %q_old = DFA.S(i);

    % Choose the correct states out of V based on the DFA
    if antagonist_pol
        [V_sort, a_pol] = min(max(squeeze(trans(i, :,:)),V_lb));

    else
        V_sort = min(max(squeeze(trans(i, :,:)),V_lb));
    end

    % Initialise value iteration with including action uhat.
    VVlb = V_sort*sysAbs.P; % P is a object of the transition_probability class.

    [~,index_pol] = max(VVlb,[],2);   % optimize over uhat
    pol(:,:,i) = uhat(:,index_pol);

end

if initialonly==true

    % Determine correct q_0 for abstract states
    policy = zeros(size(uhat,1),length(sysAbs.states));
    for i = 1:length(sysAbs.states)
        policy(:,i) = pol(:,i,DFA.trans(DFA.S0,sysAbs.labels(i)));
    end
    pol = policy;

end

if antagonist_pol
    varargout{end+1} = a_pol;
end
if upperBound
    varargout{end+1} = satProb_ub;
end

disp('----> Finish synthesizing a robust controller');
end

