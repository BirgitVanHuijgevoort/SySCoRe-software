% load Van der Pol model
% author: Birgit van Huijgevoort

%% This code loads a class object for a nonlinear model of the Van der Pol model
% the model has 2 dimensions and 1 input
% the input should be designed to
% achieve an objective

% extra parameters
tau = 0.1; % sampling time
theta = 1;

C =  eye(2);
Bw = 0.2*eye(2);
mu = zeros(2,1);
sigma = eye(2);
dim = size(Bw, 2);

param.tau = tau;
param.theta = theta;

sysNonLin = NonlinModel(@vanderpol_f,C,Bw,mu,sigma,param);

% symbolic version
x = sym('x', [sysNonLin.dim 1]);
u = sym('u', [1 1]);
f1sym = x(1)+x(2)*sysNonLin.param.tau;
f2sym = x(2)+(-x(1)+(1-x(1)^2)*x(2))*sysNonLin.param.tau+u;
sysNonLin.fsym = [f1sym; f2sym];


%% The post function
function xp = vanderpol_f(x, u,sys)

    if nargin ~= 3
        error('Invalid input !');
    end

    if size(x,1)~=2
        error('Wrong size of x');
    end
    xp = zeros(size(x));

    xp(1,:) = x(1,:) + sys.param.tau *x(2,:);
    xp(2,:) = x(2,:)+ (-x(1,:) + (1-x(1,:).^2) .* x(2,:)) .* sys.param.tau + u;
end
