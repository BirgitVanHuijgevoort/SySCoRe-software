function rel = CombineSimRel(rel_1,rel_2,sysLTIr,sysAbs)
%COMBINESIMREL Combines two simulation relations. Currently only used for
%model-order reduction.
%
% Inputs:
% -------
% rel_1, rel_2 = simulation relations to be combined
% sysLTIr = reduced-order LTI system
% sysAbs = finite-state abstraction of sysLTIr
%
% Outputs:
% --------
% rel = combined simulation relation
%
% Copyright, 2022 Birgit van Huijgevoort, b.c.v.huijgevoort@tue.nl

% Compute total delta and epsilon
delta = rel_1.delta+rel_2.delta;
epsilon = rel_1.epsilon+rel_2.epsilon;

% Combine simulation relation
rel = rel_1.Combine(rel_2,sysLTIr.X);
%disp(['delta = ', num2str(rel.delta), ', epsilon = ', num2str(rel.epsilon) ])

% additional non deterministic labelling step to take combined simulation
% relation into account.
rel.NonDetLabels = NonDeterministicLabelling(sysAbs.outputs, sysLTIr.regions, rel);


end