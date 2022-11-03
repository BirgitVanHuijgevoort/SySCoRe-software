function [K] = ComputeTaylorAccuracy(Partition,R)
% COMPUTETAYLORACCURACY computes the accuracy of the taylor linearization
%
% Inputs:
% Partition = partition of a PWA system, see Models/PWAmodel for contents
% R = symbolic expression for the remainder computed via PreSynthesisComputations/ComputeTaylorAndRemainder 
%
% Outputs:
% K = set for the error introduced by performing linearization using Taylor
% 
% Details can be found in :
% B. C. van Huijgevoort and S. Haesaert. Temporal logic control
% of nonlinear stochastic systems using a piecewise-affine abstraction.
% 2022. https://www.sofiehaesaert.com/assets/Research/PWA_abstractions.pdf
%
% Example: see Tutorials/VanderPol
%
% Copyright 2022 Birgit van Huijgevoort b.c.v.huijgevoort@tue.nl 

%%
    xc = Partition.Polyhedron.chebyCenter; % center of partition
    n = size(R.func,1); % dimension of system

    %Wrong syntax for fmincon, we need @(x)
    %[X,val] = fmincon(R.Rem2, Partition.Polyhedron.V(1,:)', Partition.Polyhedron.A, Partition.Polyhedron.b)
    
    % Option 2
    % Evaluate the function for many numbers and compute kappa_i_max
    % (maximum approximation error)    
    np = 10; % number of points for evaluation in all directions
    Xsp = cell(1,n);
    for i = 1:n
        points = linspace(min(Partition.Polyhedron.V(:,i)),max(Partition.Polyhedron.V(:,i)),np);
        Xsp(1,i) = {points};
    end

    Xpoints = combvec(Xsp{:})';

    % find maximum error
    kappa_i_max = zeros(n,1);
    for i = 1:np^2
        fx = R.Rem(xc.x',Xpoints(i,:));
        [ind,~] = find(fx>kappa_i_max);
        if ~isempty(ind)
            kappa_i_max(ind) = fx(ind);
        end
    end
    
    % Save matrix K as a Polyhedron
    A = [eye(n); -eye(n)];
    b = repmat(kappa_i_max,[n,1]);

    K = Polyhedron('A', A, 'b', b);
end