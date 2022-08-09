% This code is created by Sadegh Soudjani
% version 1, date: 22 June 2021
% prepared as part of the ARCH workshop stochastic category
% The purpose is to provide a setting of the 3d model of a robot/car

% the input of model is fixed for now but it  should be synthesized with
% the tool

clc
clear
close all

% Time horizon (for the continous version of the system)
T = 20;
% Sampling  time for discretisation
tau = 0.1;
% Number of discrete time steps
N = ceil(T/tau);
% Initial condition of the system
X0 = [0;0;0];

% velocity is fixed
V = 1; % V = 1 m/s

% support of the noise distribution is [-w1,w1]*[-w2,w2]*[-w3,w3]
% these values should be proportional with sqrt(tau)
w1 = 0.02*sqrt(tau);
w2  = 0.02*sqrt(tau);
w3 = 0.02*sqrt(tau);

% range of the state variables x1\in[x1l,x1u], x2\in[x2l,x2u], x3\in[x3l,x3u]
x1l = -0.6;
x1u = 0.6;
x2l = -1.2;
x2u = 1.98;
x3l = -pi;
x3u = pi;

% range of the inout u\in[ul,uu],
ul = -20;
uu = 20;

%% specification
% We consider three polytopes AX<=B1, AX<=B2, and AX<=B3
% AX<=B1, AX<=B2 specifies polytopes around the limit cycle of the
% system. The specification is to compute the probability that the
% trajectory goes outside of this area around the limit cycle in the time
% internal [0,T] after entering the polytope AX<=B3


%% computing the limit cycle using the deterministic system

% input
u = 0.02*ones(1,N+1);

% initialise the state vectors for Matlab
x = zeros(1,N+1);
y = zeros(1,N+1);
theta = zeros(1,N+1);

% the initial state
x(1) = X0(1);
y(1) = X0(2);
theta(1) = X0(3);

% deterministic dynamics
for k=1:N
    if(u(k) ~= 0)
        x(k+1) = x(k) + (V/u(k))*(sin(theta(k)+u(k)*tau)-sin(theta(k)));
        y(k+1) = y(k) - (V/u(k))*(cos(theta(k)+u(k)*tau)-cos(theta(k)));
        theta(k+1) = theta(k)+u(k)*tau;
    end
    if(u(k) == 0)
        x(k+1) = x(k) + V*tau*cos(theta(k));
        y(k+1) = x(k) + V*tau*sin(theta(k));
        theta(k+1) = theta(k)+u(k)*tau;
    end
end


%% Monte Carlo simulation

% Number of trajectories for Monte Carlo simulation
Ns = 100;


% Initialise the simulations
xn = zeros(Ns,N+1);
yn = zeros(Ns,N+1);
thetan = zeros(Ns,N+1);

xn(1:Ns,1) = X0(1);
yn(1:Ns,1) = X0(2);
thetan(1:Ns,1)= X0(3);
    
for i=1:Ns
    for k=1:N 
        if(u(k) ~= 0)
        xn(i,k+1) = xn(i,k) + (V/u(k))*(sin(thetan(i,k)+u(k)*tau)-sin(thetan(i,k)))+w1*(2*rand-1);
        yn(i,k+1) = yn(i,k) - (V/u(k))*(cos(thetan(i,k)+u(k)*tau)-cos(thetan(i,k)))+w2*(2*rand-1);
        thetan(i,k+1) = thetan(i,k)+u(k)*tau+w3*(2*rand-1);
        end
        if(u(k) == 0)
            xn(i,k+1) = xn(i,k) + V*tau*cos(thetan(i,k))+w1*(2*rand-1);
            yn(i,k+1) = x(i,k) + V*tau*sin(thetan(i,k))+w2*(2*rand-1);
            thetan(i,k+1) = thetan(i,k)+u(k)*tau+w3*(2*rand-1);
        end
    end
end


disp(['The number of trajectories for Monte Carlo Simulation is: ', num2str(Ns)]);
disp(' ');

%% plots
% plot sample trajectories
figure,
h1 = plot(xn,yn,'-.','color','r','linewidth',1);

% plot the deterministic trajectory
hold on,
h2 = plot(x,y,'linewidth',3,'color','b');

% end of the code
