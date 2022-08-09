% Copyright 2022 Oliver Schoen
function q_next = updateDFA(DFA, regions, APs, q, ynext)
%updateDFA Returns the updated state of the DFA.
% 
% DFA: Struct DFA of predefined form (See TranslateSpec.m).
% 
% regions: Cell array defining the regions corresponding to the AP labels,
% e.g., {p1x, p1y; p2x, p2y; p3x, p3y}. Every row corresponds to one AP
% label.
% 
% APs: Labels of the regions used in the specification, i.e., on the
% vertices of the DFA. E.g., {'p1', 'p2', 'p3'}.
% 
% q: Current state of the DFA. Must be an integer.
% 
% ynext: Observed output of the system which triggers the update of the
% DFA.
% 

assert((rem(q, 1) ==0) && numel(q) == 1, "The current state q of the DFA must be an integer.")

% Check in which region the system output is
assert(size(regions, 1) == length(APs), "Number of regions and atomic proposition should be the same.")
inRegion = zeros(1, size(regions, 1));
for i = 1:size(regions, 1)
    % Check when in region i
    inRegion(i) = inpolygon(ynext(1,:), ynext(2,:), ...
        regions{i, 1}, regions{i, 2});
end

% Get active AP labels
if sum(inRegion) == 0
    % No active AP (activeAPs = ' ')
    activeAct = findActiveAction({' '}, DFA);
else
    % Retrieve active APs
    activeAPs = APs(logical(inRegion));
    activeAct = findActiveAction(activeAPs, DFA);
end

% Progress DFA
q_next = DFA.trans(q, activeAct);

% Catch if q_next == 0
if q_next == 0
    q_next = q;
end
end

function activeAct = findActiveAction(activeAPs, DFA)
% Find DFA active action / determine which path to take
act = DFA.act;
assert(~isempty(act), "Action 's' is empty. Try ' ' instead.", act)
score = zeros(1, numel(DFA.act));
for i = 1:numel(activeAPs)
    act = erase(act, activeAPs{i});
    for j = 1:numel(act)
        if isempty(act{j})
            score(j) = i;
        end
    end
end

[~, activeAct] = max(score);
assert(~isempty(activeAct), "No action '%s' found. DFA might be incomplete.", [activeAPs{:,:}])
end