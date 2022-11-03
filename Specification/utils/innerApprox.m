function [lb, ub] =  innerApprox(Polytope)
% This function implements the inner approximation box introduced in 
% Bemporad, Alberto, Carlo Filippi, and Fabio D. Torrisi. 
% "Inner and outer approximations of polytopes using boxes." 
% Computational Geometry 27.2 (2004): 151-178.
% Given a Polyhedron Ax\geq b it solves the optimization problem
%
% \max_x,y \sum_{j in D} ln y_j
% s.t. Ax+A^+y\geq b
%
% for which the optimizers x^* and y^* define the inner bounding box as 
% Polyhedron('lb', x^*, 'ub', x^*+y^*). 

A  = Polytope.A;
b = Polytope.b;

Apos = max(A, zeros(size(A))); % keep only the positive elements of A
n = size(A,2);
x = sdpvar(n,1);
y = sdpvar(n,1);

tic
options = sdpsettings('verbose',0,'quadprog.maxiter',100, 'cachesolvers', 1);

optimize(A*x+Apos*y <= b, -sum(log(y)), options);
lb = value(x);
ub = value(x)+value(y);
