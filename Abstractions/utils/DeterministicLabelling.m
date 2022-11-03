function sysAbs = DeterministicLabelling(sysAbs, DFA, regions, APs)
%DETERMINISTICLABELLING Labels the output space of a finite-state system
%based on the determinstic finite-state automaton (DFA)
%
% Inputs
% ------
% sysAbs: finite-state system abstraction of predefined form (see Gridding.m)
% DFA: Struct DFA of predefined form (see TranslateSpec.m)
% regions: array of polehedra
% APs: Labels of the regions used in the specification, i.e., on the
% vertices of the DFA. E.g., {'p1', 'p2', 'p3'}.
%
% Outputs
% -------
% sysAbs with labels stored in sysAbs.labels
%
% Example
% --------
% P1 = Polyhedron([4, -4; 4, 0; 10, 0; 10 -4]); 
% P2 = Polyhedron([4, 0; 4, 4; 10, 4; 10 0]);  
% regions = [P1;P2]; 
% APs = {'p1', 'p2'};
% formula = '(!p2 U p1)';  
% [DFA] = TranslateSpec(formula,APs);
% sysAbs = DeterministicLabelling(sysAbs, DFA, regions, APs);
%
% Copyright 2022 Oliver Schoen, Birgit van Huijgevoort

assert(size(regions, 1) == length(APs), "Number of regions and atomic proposition should be the same.")
checkDFAact(DFA) % Check validity of DFA.act

% Create empty labels
l = sysAbs.l;
assert(~isempty(l), "Please supply sysAbs with the number of grid cells in sysAbs.l.")
label = ones(1, prod(l)); % Corresponds to the empty label

%% Check when in regions
% Switch w.r.t. dimension of system
dim = size(sysAbs.outputs, 1);
nRegions = size(regions, 1);
nPoints = size(sysAbs.outputs, 2);
if dim == 1
    % 1) 1D case
    inRegion = zeros(nRegions, nPoints);
    for i = 1:size(regions, 1)
        % Check which points are contained in polyhedron i
        inRegion(i, :) = (sysAbs.outputs > min(regions(i).V)) & ...
            (sysAbs.outputs < max(regions(i).V));
    end
else
    % 2) Multidimensional case
    inRegion = EfficientContainment(regions, sysAbs);
end
for i = 1:size(regions,1)
    % Find corresponding label in DFA
    AP = APs{i};
    idx = strcmp(DFA.act, AP);
    assert(sum(idx) == 1, "Atomic proposition %s must be unambigously defined in the DFA.", AP)
    [label(logical(inRegion(i,:)))] = deal(find(idx));
end

% Return labels
sysAbs.labels = label;
end
