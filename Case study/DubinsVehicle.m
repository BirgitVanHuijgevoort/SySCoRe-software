% load Dubins' car model
% author: Birgit van Huijgevoort

%% This code loads a class object for a nonlinear model of the Dubins' car model
% the model has 3 dimensions and 1 input
% the input should be designed to
% achieve an objective

% extra parameters
tau = 1; % sampling time
V = 0.1; % fixed velocity
    
C =  eye(3);
Bw = 0.06*eye(3);
mu = zeros(3,1);
sigma = eye(3);
dim = size(Bw, 3);

param.tau = tau;

sysNonLin = NonlinModel(@DubinsCar_f,C,Bw,mu,sigma,param);

% symbolic version
x = sym('x', [sysNonLin.dim+1 1]);  % for Taylor expansion we consider u as an additional state
f1sym = x(1)+(1/x(4))*V*sin(x(3)+x(4)*sysNonLin.param.tau)-(1/x(4))*V*sin(x(3));
f2sym = x(2)-(1/x(4))*V*sin(x(3)+x(4)*sysNonLin.param.tau)+(1/x(4))*V*cos(x(3));
f3sym = x(3)+x(4)*sysNonLin.param.tau;
sysNonLin.fsym = [f1sym; f2sym; f3sym];


%% The post function
function xp = DubinsCar_f(x, u, sys)

    if nargin ~= 3
        error('Invalid input !');
    end
    
    if size(x,1)~=2
        error('Wrong size of x');
    end
    xp = zeros(size(x));

    xp(1,:) = x(1,:) + (1/u)*V*sin(x(3,:)+u*sysNonLin.param.tau)-(1/u)*V*sin(x(3,:));
    xp(2,:) = x(2,:)-(1/u)*V*sin(x(3,:)+u*sysNonLin.param.tau)+(1/u)*V*cos(x(3,:));
    xp(3,:) = x(3,:)+u*sysNonLin.param.tau;
end
    
