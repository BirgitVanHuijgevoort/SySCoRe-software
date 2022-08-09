% Copyright 2022 Oliver Schoen
function augProdGMDP = genAugProdGMDP(sysAbs, LDBA, xi)
%genAugProdGMDP

% Check whether automaton is a LDBA (accepting transitions)
% TODO

% Default xi
if nargin < 3
    xi = 0.1;
end

%% For every accepting transition Acc_i, determine the outgoing DFA location
% q_i and the triggering action labels act_i
idx_Acc = find(LDBA.Acc); % Index of accepting transitions in matrix
[q_Acc, act_Acc] = ind2sub(size(LDBA.Acc), idx_Acc); % Determine accepting locations and action labels
augProdGMDP.q_Acc = q_Acc;
augProdGMDP.act_Acc = act_Acc;

%% Determine accepting states of product gMDP
% -> in synthesize controller function

%% Add sink location 
% Sink states (x, q)_i of the actual LDBA are determined when synthesizing
% the controller
assert(isrow(LDBA.S), "LDBA.S must be a row vector!")
nextIdxS = max(LDBA.S) + 1; % Next location index
LDBA.S = [LDBA.S, nextIdxS]; % Add additional location
LDBA.sink{end+1} = nextIdxS; % Label new location as sink location
LDBA.trans(end+1, :) = deal(nextIdxS); % To transitions from sink location (could be removed but kept for understandability)

%% Build augmented product gMDP
augProdGMDP.sys = sysAbs;
augProdGMDP.aut = LDBA;
augProdGMDP.S_sink = nextIdxS; % Label as accepting location
augProdGMDP.xi = xi; % Probability of jumping from F to phi is (1-xi)
end