% Copyright 2022 Oliver Schoen
function q_next = updateDFA(DFA, regions, APs, q, ynext)
%updateDFA Returns the updated state of the DFA.
% 
% DFA: Struct DFA of predefined form (See TranslateSpec.m).
% 
% regions: array of polyhedra defining the regions corresponding to the AP labels,
% Every row corresponds to one AP label.
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
    inRegion(1,i) = inRectangle(ynext,regions(i));
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
% Find DFA active action(s) / determine which path to take
% Feature: If multiple APs active identifying joint AP in DFA.act
act = DFA.act;
score = zeros(1, numel(act));
for i = 1:numel(activeAPs)
    temp = strfind(act, activeAPs{i});
    for j = 1:numel(temp) % Replace empty with 0 and values greater 0 with 1
        if isempty(temp{j})
            temp{j} = 0;
        elseif temp{j} > 0
            temp{j} = 1;
        end
    end
    score = score + [temp{:}]; % AP i found in these places in DFA.act
end

% Choose AP with highest score
[~, activeAct] = max(score);

assert(~isempty(activeAct), "No action '%s' found. DFA might be incomplete.", [activeAPs{:,:}])
end