function [d] = R2G(R,A,BetaSet)
%R2G computes the minimal size of Gamma needed to compensate for the
%disturbance Beta so that the dynamical systems x(t+1) = Ax+b with b in
%Beta stays in R. 

% Todo:
% Check whether the inputs are matrices or zonotopes/polytopes. 

ARB = A*R+BetaSet;
Rmin = ARB&R;
d = hausdorffDist(Rmin,ARB);

end

