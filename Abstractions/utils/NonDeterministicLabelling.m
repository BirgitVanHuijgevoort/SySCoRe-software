function [state2act] = NonDeterministicLabelling(outputs,Polytopes,rel, varargin)
%NonDeterministicLabelling Map the states to the labelling number
%
% Outputs: 
% state2act(i,j) is 1 if the label DFA.act{i} holds for state
% sysAbs.states(:,j) and 0 if the label does not hold.
%
% Inputs:
% - outputs(:,i) = is the output of the system corresponding to state sysAbs.state(:,i) 
% - Polytopes = regions used for labeling the output space
% - rel = simulation relation, with rel.epsilon the output deviation
%
% Options:
% 'Efficient' = false/true
%
% Copyright 2022 Sofie Haesaert s.haesaert@tue.nl

efficient = false;
  for i = 1:length(varargin)
        % try to find 'Efficient'
        if strcmp(varargin{i},'Efficient')
            efficient = true;
            sysAbs = varargin{i+1}; % should be sysAbs.zstates  or sysAbs.states 
        end
  end
  % remark sysAbs.outputs = sys.C*XhatSpace;

epsil = rel.epsilon;
Polytope_array_large = [];
Polytope_array_small = [];
Polytope_containment_large = [];
Polytope_containtment_small = [];

n_p = length(Polytopes); % number of polytopes
for p = 1:n_p
    [PolytopeLarge,PolytopeSmall] = IncreaseDecreasePolytope(Polytopes(p), epsil);
    
    Polytope_array_large = [Polytope_array_large;PolytopeLarge.computeVRep];
    Polytope_array_small = [Polytope_array_small;PolytopeSmall.computeVRep];

end

assert(~isempty(PolytopeLarge.V), "PolytopeLarge of polytope %d has no vertices.", p)
assert(~isempty(PolytopeSmall.V), "The chosen epsilon is too big. PolytopeSmall of polytope %d has no vertices.", p)

if ~efficient
    % for all states compute whether they are contained in the polytopes
    Polytope_containment_large = Polytope_array_large.contains(outputs);
    Polytope_containment_small = Polytope_array_small.contains(outputs);
else
    Polytope_containment_large = EfficientContainment(Polytope_array_large, sysAbs);
    Polytope_containment_small = EfficientContainment(Polytope_array_small, sysAbs);
end

state2act = [];
bin_vals = dec2bin(0:2^n_p-1);
for index = 1:2^n_p
    el = bin_vals(index, :);
    x_true = [];
    for p_index = 1:n_p
        if el(p_index)=='0'
            % could polytope p not hold?
            % check this for all states b looking at the smaller version of
            % polytope p
            x_true = [x_true; 1-Polytope_containment_small(p_index,:)];
        else
            x_true = [x_true; Polytope_containment_large(p_index,:)];

        end
    end
 
    if size(x_true,1)>1
        state2act = [state2act; all(x_true)];
    else
        state2act = [state2act; x_true];
    end
 
end
end
