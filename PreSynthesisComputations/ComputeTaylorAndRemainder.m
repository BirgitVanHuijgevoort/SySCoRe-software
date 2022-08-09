function [T,R,NL] = ComputeTaylorAndRemainder(sysNonLin)
% function that compute the second order Taylor expansion and corresponding
% remainder of the nonlinear system. 
% This function only works for f(x) being a polynomial function with
% maximum degree 3.
% Note that this function requires the symbolic toolbox!

n = sysNonLin.dim;
x = sym('x', [n 1]);
c = sym('c', [n 1]); % linearization point

% PD(i,j) contains the polynomial degree of function i wrt variable x(j)
% NL(i,1) contains 1 if function i is nonlinear
PD = zeros(n,n);
NL = ones(n,1);

% Check if the function is already linear
for i = 1:n
    % Compute polynomial degree of function i
    for j = 1:n
        % with respect to variable x(j)
        degree = polynomialDegree(sysNonLin.fsym(i),x(j));
        PD(i,j) = degree;
    end
    
    NL(i,1) = 1-all(PD(i,:)<2);
end


if all(all(PD == 1))
    disp('This system is already linear!')
else
    i = find(NL == 1);
    fprintf('Pre-computation of piecewise-affine approximation for state update of %s \n', string(x(i)))
end

fsym = sysNonLin.fsym;

%% Linearize using Taylor expansion around the point c
for i = 1:n
    if NL(i) == 0
        T(i,1) = fsym(i);
    elseif NL(i) == 1
        T(i,1) = taylor(fsym(i), x, 'ExpansionPoint', c, 'Order', 2);
        T(i,1) = simplify(T(i));
    end
end

% Get correct format for subs
X = [];
C = [];
for k=1:n
    X = [X,x(k)];
    C = [C,c(k)];
end

% construct affine approximation of format
% Ax+a

% Compute affine term by substituting zero for the variables
a = subs(T,X,zeros(1,n));

% Compute A matrix 
A = sym('A', [n n]);
for i=1:n
    for j = 1:n
        I = zeros(1,n);
        I(1,j) = 1;
        A(i,j) = subs(T(i),X,I)-a(i);
    end
end

% Save Taylor expansion as a matlab function
Tay.func = matlabFunction(T);
Tay.A = matlabFunction(A, 'Vars', {C});
Tay.a = matlabFunction(a, 'Vars', {C});

%% Compute remainder
% Compute error caused by only using 1st order Taylor expansion
for i = 1:n
    Tfull(i,1) = taylor(fsym(i), x, 'ExpansionPoint', c,'Order', max(PD(i,:))+1);
end

E = Tfull-T;

R_vec = num2str(zeros([n 1]));
R_tot = str2sym(R_vec);
for i = 1:n
    if E(i)~=0 
    
        theta = sym('theta', [n 1]);

        % Compute error caused by Taylor expansion in general
        if max(PD(i,:))+1 == 1
            % Compute gradient
            G = gradient(fsym(i),x); 
            % Substitute c+theta.*(x-c) into gradient
            D = subs(G, x, c+theta.*(x-c));

            XminC = [(x-c)]';

            % Compute remainder 
            R = (1/factorial(max(PD(i,:))+1))*XminC*D;
        elseif max(PD(i,:))+1 == 2
            % Compute gradient and hessian
            G = gradient(fsym(i),x);
            H = hessian(fsym(i),x);

            % Substitute c+theta.*(x-c) into gradient
            D = subs(H, x, c+theta.*(x-c));

            % Compute remainder
            R = (1/factorial(max(PD(i,:))+1))*([x-c]'*H*[x-c]);

        elseif max(PD(i,:))+1 == 3  
            % Compute gradient and hessian
            G = gradient(fsym(i),x);
            H = hessian(fsym(i),x);
            
            % Compute derivatives of hessian
            Gm = zeros(n,n^2);
            for j = 1:n
                % j=1
                Gm(:,(j-1)*n+1:j*n) = jacobian(H(:,j),x);
            end
            
            % Substitute c+theta.*(x-c) into gradient
            Dm = subs(Gm, x, c+theta.*(x-c));
            
            % construct matrix to compute (a+b+c)^3
            azerob = num2str(zeros(n,n^2));
            azerob = str2sym(azerob);
            for i = 1:n
                block = diag(repmat(x(i)-c(i),[1 n]));
                azerob(:,(i-1)*n+1:i*n) = block;
            end

            % terms in (a+b+c)^3
            terms = (x-c)*(x-c)'*azerob;

            % Compute remainder
            R = (1/factorial(max(PD(i,:))+1))* ...
                sum(terms.*Dm, 'all');
        else
            error('PWA-approximation not implemented for nonlinear functions with a polynomial degree higher than 3')
        end

        % Add errors together for total remainder function
        R_total = abs(simplify(E(i)+R));
        R_tot(i) = R_total;
    end
end

clear R

% Save remainder as a matlab function
Rem = matlabFunction(R_tot, 'Vars', {[c(1), c(2)], [x(1),x(2)]});
R.func = R_tot;
R.Rem = Rem;

T = Tay;
end

