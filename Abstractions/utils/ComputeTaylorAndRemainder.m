function [T,R,NL] = ComputeTaylorAndRemainder(sysNonLin)
% function that compute the second order Taylor expansion and corresponding
% remainder of the nonlinear system. 
% Note that this function requires the symbolic toolbox!
%
% Outputs:
% --------
% NL(i,j) indicates if state dynamics x_i(t+1) is nonlinear wrt x (j=1) and
% wrt u (j=2)

nx = sysNonLin.dim;
nu = size(sysNonLin.U.V,2);

x = sym('x', [nx 1]);
u = sym('u',[nu,1]);
c = sym('c', [nx+nu 1]); % linearization point

% PD(i,j) contains the polynomial degree of function i wrt j-th variable in
% [x, u]
PD = zeros(nx,nx);
NL = ones(nx,2);

% Check if the function is already linear (wrt x)
for i = 1:nx
    % Compute polynomial degree of function i (if possible)
    for j = 1:nx
        % with respect to variable x(j)
        try 
            degree = polynomialDegree(sysNonLin.fsym(i),x(j));
        catch
            degree = 10; % any number > 1 leads to nonlinear. 
        end
        PD(i,j) = degree;
    end    
    NL(i,1) = 1-all(PD(i,:)<nx);
end

% Check if the function is already linear (wrt u)
clear PD
PD = zeros(nx,nu);
for i = 1:nx
    % Compute polynomial degree of function i
    for j = 1:nu
        % with respect to variable u(j)
        try 
            degree = polynomialDegree(sysNonLin.fsym(i),u(j));
        catch
            degree = 10;
        end
        PD(i,j) = degree;
    end    
    NL(i,2) = 1-all(PD(i,:)<2);
end

if all(all(NL == 0))
    disp('This system is already linear!')
else
    [i, ~] = find(NL == 1);
    i = unique(i);
    fprintf('Pre-computation of piecewise-affine approximation for state update of %s \n', string(x(i)))
end

fsym = sysNonLin.fsym;

%% Linearize using Taylor expansion around the point c
for i = 1:nx
    if all(NL(i,:) == 0)
        T(i,1) = fsym(i);
    elseif any(NL(i,:) == 1)
        T(i,1) = taylor(fsym(i), [x;u], 'ExpansionPoint', c, 'Order', 2);
        T(i,1) = simplify(T(i));
    end
end

% Get correct format for subs
X = [];
for k=1:nx
    X = [X,x(k)];
end

U = [];
for k = 1:nu
    U = [U,u(k)];
end

C = [];
for k=1:nx+nu
    C = [C,c(k)];
end

% construct affine approximation of format
% Ax+Bu+a

% Compute affine term by substituting zero for the variables
Tx = subs(T,X,zeros(1,nx));
a = subs(Tx,U,zeros(1,nu));

% Compute A matrix 
A = sym('A', [nx nx]);
for i=1:nx
    for j = 1:nx
        I = zeros(1,nx);
        I(1,j) = 1;
        Tx(i) = subs(T(i),U,zeros(1,nu));
        A(i,j) = subs(Tx(i),X,I)-a(i);
    end
end

% Compute B matrix
B = sym('B', [nx,nu]);
for i=1:nx
    for j = 1:nu
        I = zeros(1,nu);
        I(1,j) = 1;
        Tu(i) = subs(T(i),X,zeros(1,nx));
        B(i,j) = subs(Tu(i),U,I)-a(i);
    end
end

% Save Taylor expansion as a matlab function
Tay.func = matlabFunction(T);
Tay.A = matlabFunction(A, 'Vars', {C});
Tay.B = matlabFunction(B, 'Vars', {C});
Tay.a = matlabFunction(a, 'Vars', {C});

%% Compute remainder
% Compute error caused by only using 1st order Taylor expansion
R_vec = num2str(zeros([nx 1]));
R_tot = str2sym(R_vec);
for i = 1:nx
    if any(NL(i,:) == 1)
        theta = sym('theta', [nx+nu 1]);
        %% Compute error caused by Taylor expansion 
        % Compute gradient and hessian
        G = gradient(fsym(i),[x;u]);
        H = hessian(fsym(i),[x;u]);

        % Substitute c+theta.*(x-c) into gradient
        D = subs(H, [x;u], c+theta.*([x;u]-c));

        % Compute remainder
        R_tot(i) = (1/2)*([[x;u]-c]'*H*[[x;u]-c]);
    else
        R_tot(i) = 0;
    end
end


% Save remainder as a matlab function
%Rem = matlabFunction(R_tot, 'Vars', {[c(1), c(2)], [x(1),x(2)]});
Rem = matlabFunction(R_tot, 'Vars', {C, X, U}); % add U
R.func = R_tot;
R.Rem = Rem;

T = Tay;
end

