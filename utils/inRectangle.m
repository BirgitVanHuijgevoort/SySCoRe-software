function [in] = inRectangle(state,poly)
%INRECTANGLE computes if the state is inside the RECTANGULAR polyhedron specified by
%poly
%
% Inputs
% ------
% state = a point (or in this case state of a system). Can be a matrix, then each column
% represents one point.
% poly = the polyhedron described by a Polyhedron of the MPT3 toolbox.
% poly.V gives the points of the vertices.
% 
% Outputs 
% -------
% in = 1 if state is in polyhedron (including boundary)
% in = 0 if state is outside of polyhedron
%
% Copyright 2022 Birgit van Huijgevoort, b.c.v.huijgevoort@tue.nl

%%
n = size(state,2);  % number of states
dim = size(state,1);    % state dimension

in = zeros(dim,n);
for i = 1:dim
    in(i,:) = state(i,:) >= min(poly.V(:,i)) & state(i,:) <= max(poly.V(:,i));
end

in = (sum(in,1) == dim);

end