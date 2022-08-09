% Copyright 2022 Oliver Schoen
function [satProp, pol] = SynthesizeRobustControllerAugmented( ...
    augProdGMDP, rel, N, initialonly)
%SynthesizeRobustControllerAugmented

%% Init
checkDFAact(augProdGMDP.aut) % Check validity of LDBA.act

% Unpack variables
sysAbs = augProdGMDP.sys;
LDBA = augProdGMDP.aut;
if sum(LDBA.trans == 0, 'all') >= 0
    % LDBA transition elements with value 0 are not supported. Substitute
    % by row number.
    for i = 1:size(LDBA.trans, 1)
        LDBA.trans(i, LDBA.trans(i, :) == 0) = deal(i);
    end
end
% assert(sum(LDBA.trans == 0, 'all') == 0, ...
%     "LDBA transition elements with value 0 are not supported. Consider substituting by row number.")
nQ = length(LDBA.S);
nQ_main = nQ - 1; % One location is the sink location (see below)

uhat = sysAbs.inputs;
delta = rel.delta;

try
    dim = sysAbs.orig.dim;
catch
    dim = sysAbs.det.dim;
end

try
    nS = sysAbs.nS;
catch
    nS = length(sysAbs.states);
end

activeLabels = rel.NonDetLabels; % Active output regions
nLabels = size(activeLabels, 1);

% Initialise value function
% V(i, j) is the probability of reaching phi from LDBA location i and 
% abstract state sysAbs.states(:, j)
%V = zeros(nQ, nS);

% Some helpful variables are introduced
LDBA_Active = setdiff(LDBA.S, augProdGMDP.S_sink); % Recognize that the sink location is not active
Converged = ones(1, nQ_main); % Keep track of converged locations
Converged(LDBA_Active) = deal(0); % Non-active states are already converged

% Set V to 1 for product-state sink state phi which is equivalent with the
% LDBA sink location here. Note that this is differing from the original
% definition, where phi is defined as a definite (xhat, q). Since both the 
% LDBA's and the system's transitions don't determine when the system
% jumps to the sink state phi but this is purely determined by the current
% product state (xhat, q) (and, if (xhat, q)!=phi, by chance) we can model
% the sink state phi by a separate LDBA state only.
% To make the computations more efficient we split V into two matrices as
% follows: V = [V_main; V_sink]
V_main = zeros(nQ_main, nS);
V_sink = ones(1, nS);
%V(augProdGMDP.S_sink, :) = deal(1); 

%% Compute robust policy
disp('Start computing robust policy...')
tic 

% Prepare LDBA transitions for all states
trans = ones(nQ_main, nQ_main, nS); % Possible transitions: impossible transitions receive value inf and zero otherwise
acctrans_nonrobust = ones(nQ_main, nS); % Accepting possible transitions: accepting possible transitions receive value of xi and 1 otherwise. Robustly accepting possible transitions receive value 0
acctrans_robust = ones(nQ_main, nQ_main, nS); % Accepting possible transitions: robustly accepting possible transitions receive value of xi and 1 otherwise
for i = LDBA_Active
    for l = 1:nLabels
        % Transitions from the i-th active location to the location reached
        % under action label l are marked with a zero. Impossible
        % transitions are masked with a value of inf. The value function
        % values for these transitions are hence later-on filtered out.
        trans(i, LDBA.trans(i, l), :) = ...
            min( ...
            shiftdim(trans(i, LDBA.trans(i, l), :), 1), ...
            1 - activeLabels(l, :) ...
            );
    end   
        % Accepting transitions
        % A) Determine for q_i for which xhat_i at least one active 
        % transitions is accepting and label those with xi.
        acc_nonrobust = LDBA.Acc(i, :) * activeLabels;
        acctrans_nonrobust(i, :) = ...
            acctrans_nonrobust(i, :) - acc_nonrobust; %acctrans_nonrobust(i, LDBA.trans(i, :), :)
        % The matrix is now 0 for every LDBA transition that is accepting
        % and 1 otherwise (robustly accepting transitions filtered out
        % later)

%         % B) Determine for q_i for which xhat_i all active transitions are
%         % accepting and label those with zero. Then copy this row vector 
%         % for every active label in q_i
%         acc_robust = LDBA.Acc(i, :)' * subplus( ...
%             LDBA.Acc(i, :) * activeLabels + 1 - sum(activeLabels, 1) ...
%             );
%         acctrans_robust(i, LDBA.trans(i, :), :) = ...
%             shiftdim(acctrans_robust(i, LDBA.trans(i, :), :), 1) - acc_robust;
%         % The matrix is now 0 for every LDBA transition that is accepting
%         % even in the worst-case (aka robustly accepting product states)
end
trans = 1e256 * trans; % inf approximated with 1e256
acctrans_nonrobust = subplus(acctrans_nonrobust); % It only matters if there is at least one accepting transition
acctrans_nonrobust(acctrans_nonrobust == 0) = augProdGMDP.xi; % Now we have only 1 and xi
% acctrans_nonrobust(acctrans_robust == 0) = 0; % Filter out robustly-accepting transitions
% acctrans_robust(acctrans_robust > 0) = acctrans_robust(acctrans_robust > 0) - augProdGMDP.xi; % Substract xi from every non-robustly-accepting element and add it in the next step
% acctrans_robust = acctrans_robust + augProdGMDP.xi; % Now we have only 1 and xi
assert(augProdGMDP.xi > 0, "Xi must be greater zero.")
assert(sum(acctrans_nonrobust < 1, 'all') > 0, "No accepting transitions. Value function will be identical zero.")
toc

% Calculate the value of every product state wrt the worst-case LDBA
% transition considering the epsilon-ball around the system output
for k = 1:N % Iterate over time steps
    % Check whether value function has converged
    if min(Converged) == 1
        disp('Value iteration has converged.')
        disp(['Number of iterations required, k = ', num2str(k)])
        break
    end

    for i = LDBA_Active(end:-1:1) % Iterate over active LDBA modes
        % Check whether mode has converged
        if min([Converged(LDBA.trans(i, :)), Converged(i)]) == 1
            % Mode has converged
            continue
        else
            % Mode not converged yet
            Converged(i) = 0;

            % Apply masking to transitions that are impossible and choose
            % worst-case value wrt all possible LDBA transitions
                % In detail: The max operation overwrites all elements of the
                % value function that are masked by a inf in the transition
                % matrix (because they are impossible considering the applied
                % epsilon-ball and the resulting output labels). Then, the min
                % operator chooses the least value possible under all LDBA
                % transitions left. Note that the min operation returns a value
                % of less than inf by definition, because there is always at
                % least one possible LDBA transition.
            V_sort = max(shiftdim(trans(i, :, :), 1), V_main);
            % Jump to phi with probability (1-xi) if q_i is an accepting
            % location and the accepting action is called (recall that the 
            % accepting transitions have to be travelled infinitely often)
%             V_sort_ = (1 - shiftdim(acctrans_robust(i, :, :), 1)) + ...
%                 shiftdim(acctrans_nonrobust(i, :, :), 1) .* V_sort;
%                 shiftdim(acctrans_nonrobust(i, :, :), 1) .* V_sort;
            factor = (1 - acctrans_nonrobust(i, :));
            V_sort_ = factor + (1 - factor) .* ...
                (min(V_sort, [], 1));

            % Apply masking by filtering out inf elements (aka impossible
            % transitions), remaining with the worst-case value possible
            % wrt possible LDBA transitions
%             V_worstCasePossible = min(V_sort_, 1);
            V_worstCasePossible = V_sort_;

            % Compute value function for each action uhat
            VV = V_worstCasePossible * sysAbs.P; % P is a object of the transition_probability class

            % Optimize over uhat and subtract delta
            V_n = max(VV, [], 2) - delta;

            % Make sure that value function is between 0 and 1
            V_n = min(1, max(0, V_n));

            % Check whether mode has converged
            if min(Converged(LDBA.trans(i, :))) == 1
                % Mode has converged
                Converged(i) = 1;
            elseif max(max(abs(V_n' - V_main(i, :)))) < 1e-256
                % Mode has converged (below tolerance)
                Converged(i) = 1;
            end

            % Update value function
            V_main(i, :) = V_n;
        end
    end
end

% Check whether value function has converged
if min(Converged) ~= 1
    fprintf('Value iteration ended without converging after k = %d iterations.\n', k)
end

V = [V_main; V_sink];


%% Compute satisfaction probability
if initialonly == true
    % Determine correct q_0 for abstract states
    satProp = V( ...
        (0:length(V) - 1) * length(LDBA.S) + ...
        LDBA.trans(LDBA.S0, sysAbs.labels) ...
        );
else
    satProp = V;
end

%% Calculate policy
warning("Policy calculation not implemented yet. Returning empty policy.")
pol = [];

end