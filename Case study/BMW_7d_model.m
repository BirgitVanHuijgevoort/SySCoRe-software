
%% This code loads a class object for a nonlinear model of the BMW_7d model
% the model has 7 dimensions and two inputs
% the input is chosen randomly for now, but it should be designed to
% achieve an objective
 
sysNonLin.type = "NonLin";
sysNonLin.f = @sys_post_7d_model;
sysNonLin.C = [1 0 0 0 0 0 0;
               0 1 0 0 0 0 0];
sysNonLin.Bw = eye(7);
sysNonLin.mu = zeros(7,1);
sysNonLin.sigma = eye(7);
sysNonLin.dim = size(sysNonLin.Bw, 2);



%% The post function
function xp = sys_post_7d_model( x, u)

    if nargin ~= 2
        error('Invalid input !');
    end 
    
    if length(x) ~= 7
        error('Invalid size of x !');
    end     
    if length(u) ~= 2
        error('Invalid size of x !');
    end          
    Ts =0.1;

    % decide a good solving step for this syste based on Ts
    % this is to be chosen imperically for each system and we use
    % matlab's ODE45 solver as a reference for error
    divider = 200;
    
    % now we have new Ts
    Ts = Ts / divider;
    
    f = zeros(7,1);
    
    for k=1:divider
        % kinematic model for small velocities ?
        if abs(x(4)) < 0.1
            %wheelbase
            lwb = 2.578912800000000; 

            % rhs of $\dot x = f(x,u)$
            f(1) = x(4)*cos(x(5));
            f(2) = x(4)*sin(x(5));
            f(3) = u(1);
            f(4) = u(2);
            f(5) = x(4)/lwb*tan(x(3));
            f(6) = u(2)*lwb*tan(x(3)) + x(4)/(lwb*cos(x(3))^2)*u(1);
            f(7) = 0;
        else

            % params
            mu = 1.048900000000000;
            C_Sf = 20.898083706740398; 
            C_Sr = 20.898083706740398; 
            lf = 1.156195706400000;
            lr = 1.422717093600000;
            h = 0.613730040000000;
            m = 1.093295233467405e+03;
            I = 1.791599530012286e+03;
            g = 9.81;

            % rhs of $\dot x = f(x,u)$
            f(1) = x(4)*cos(x(7) + x(5));
            f(2) = x(4)*sin(x(7) + x(5));
            f(3) = u(1);
            f(4) = u(2);
            f(5) = x(6);
            f(6) = -mu*m/(x(4)*I*(lr+lf))*(lf^2*C_Sf*(g*lr-u(2)*h) + lr^2*C_Sr*(g*lf + u(2)*h))*x(6) ...
                +mu*m/(I*(lr+lf))*(lr*C_Sr*(g*lf + u(2)*h) - lf*C_Sf*(g*lr - u(2)*h))*x(7) ...
                +mu*m/(I*(lr+lf))*lf*C_Sf*(g*lr - u(2)*h)*x(3);
            f(7) = (mu/(x(4)^2*(lr+lf))*(C_Sr*(g*lf + u(2)*h)*lr - C_Sf*(g*lr - u(2)*h)*lf)-1)*x(6) ...
                -mu/(x(4)*(lr+lf))*(C_Sr*(g*lf + u(2)*h) + C_Sf*(g*lr-u(2)*h))*x(7) ...
                +mu/(x(4)*(lr+lf))*(C_Sf*(g*lr-u(2)*h))*x(3); 
        end    
        % post state for the given Ts
        xp = x + Ts.*f;
        x = xp;
    end
    
end

function Simulate_bmw(Sys_d)
N = 100000; % number of time steps

% I am not sure  what input is the best to get a meaningful somulation,
% here is some random input to get the simulation going.
us = zeros(2,N);
us(1,:) = 0.01*randn(1,N);
us(2,:) = 0.01+0.01*randn(1,N);
xs = zeros(7,N);

% initial state
k=1;
xs(1,k) = 10; % Initial position coordinate
xs(2,k) = 4; % Initial position coordinate
xs(3,k) = 0.1; % Initial steering angle
xs(4,k) = 15; % Initial heading velocity
xs(5,k) = 0.2; % Initial yaw angle
xs(6,k) = 0.7; % Initial yaw rate
xs(7,k) = 0.01; % Initial slip angle


% simulating the model for N time steps
for k=1:N-1
    xs(:,k+1) = Sys_d.f( xs(:,k), us(:,k));
end


plot(xs(1,:),xs(2,:))

end
