% load Van der Pol model



%% This code loads a class object for a nonlinear model of the Van der Pol model
% the model has 2 dimensions and 1 input
% the input should be designed to
% achieve an objective

% extra paramaters
    tau = 0.1; % sampling time
    theta = 1;
    
C =  eye(2);
Bw = eye(2);
mu = zeros(2,1);
sigma = eye(2);
dim = size(Bw, 2);

param.theta = theta;
param.tau = tau;

sysNonLin = NonlinModel(@vanderpol_f,C,Bw,mu,sigma,param);




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
    xp(2,:) = x(2,:)+ (-x(1,:) + sys.param.theta .* (1-x(1,:).^2) .* x(2,:)) .* sys.param.tau + u;
    end
    
