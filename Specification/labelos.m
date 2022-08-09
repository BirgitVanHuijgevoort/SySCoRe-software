% Copyright 2022 Oliver Schoen
function sysAbs = labelos(sysAbs, DFA, regions, APs)
%labelss Labels the output space of the system abstraction.
% 
% sysAbs: System abstraction of predefined form (see Gridding.m).
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
assert(size(regions, 1) == length(APs), "Number of regions and atomic proposition should be the same.")
checkDFAact(DFA) % Check validity of DFA.act

% Create empty labels
l = sysAbs.l;
label = ones(1, prod(l)); % Corresponds to the empty label

% Check when in regions
assert(size(regions, 1) == length(APs), "Number of regions and atomic proposition should be the same.")
for i = 1:size(regions, 1)
    % Check when in region i
    inRegion = inpolygon(sysAbs.outputs(1,:), sysAbs.outputs(2,:), ...
        regions{i, 1}, regions{i, 2});
    % Find corresponding label in DFA
    AP = APs{i};
    idx = strcmp(DFA.act, AP);
    assert(sum(idx) == 1, "Atomic proposition %s must be unambigously defined in the DFA.", AP)
    label(inRegion) = deal(find(idx));
end

% Return labels
sysAbs.labels = label;
end