classdef PWAModel
% PWAMODEL Class of piecewise-affine systems with noise on the transitions. 
% Objects of this class define a piecewise-affine model with dynamics
% x(t+1) = A_i x(t) + B_i u(t) + a_i + B_w w(t); y(t)   = C_i x(t), with i\in[1,Np]
% that approximates the nonlinear model from NonLinModel with function f(x(t)).
% Here, w(t) ~ N(mu, sigma), with N a normal distribution with mean mu
% and sigma the noise variance. 
%
% Examples 
% A simple 2 dimensional model with one control input and one
% output:
%   C = [1, 1]; 
%   Bw = [0; 1]; 
%   model = PWAModel(@vanderpol_PWA_f,Partition,C,Bw,[0;0],eye(2));
%
% with function vanderpol_f defined below.
%   function xp = vanderpol_f(x, u,sys)
%     if nargin ~= 3
%         error('Invalid input !');
%     end
% 
%     if size(x,1)~=2
%         error('Wrong size of x');
%     end
%     xp = zeros(size(x));
% 
%     xp(1,:) = x(1,:) + sys.param.tau *x(2,:);
%     xp(2,:) = x(2,:)+ (-x(1,:) + (1-x(1,:).^2) .* x(2,:)) .* sys.param.tau + u;
%     end
% 
% with Partition(i) defined below for i = 1
% Partition(1).Polyhedron = ([combvec([-3 -2.85], [-3, -2.85])]')
% Partition(1).Dynamics.A = [1 0; 0 1];
% Partition(1).Dynamics.a = [1; 1]; 
% Partition(1).Dynamics.B = [0;1];
% Partition(1).Dynamics.C = [1 0];
% Partition(1).Dynamics.D = 0;
% Partition(1).Dynamics.Bw = [0;1];
% Partition(1).Dynamics.dim = sysNonLin.dim;
%  
% Copyright 2022 Birgit van Huijgevoort, b.c.v.huijgevoort@tue.nl
    
    properties
        type = "PWA";
        f % the deterministic function for the state update
        Partition % struct containing the dynamics and polyhedrons of the partitions
        C % C matrix to map from state to to output space
        Bw % Bw matrix for the noise
        Np % number of partitions
        mu % noise mean
        sigma % noise variance
        dim  % the dimension of the state space
        orig % original nonlinear model
        N % number of partitions in each direction used to obtain PWA approximation
        
        X % state space
        U % input space 
        regions % regions labelled with atomic propositions
        AP % atomic propositions
        P % the probability transition matrix. (Class: TransitionProbability)
        
    end
    
    methods
        function obj = PWAModel(Partition,C,Bw,mu,sigma)
            %PWAMODEL Construct an instance of this class
            %   Load all values f, Partition, C, ... as described in the
            %   class file. 
            %   model = PWAModel(f,Partition,C,Bw,mu,sigma) Loads model dynamics.
            % 
            % 
            % Examples
            % A simple 2 dimensional model with one control input and one
            % output:
            %   C = [1, 1]; 
            %   Bw = [0; 1]; 
            %   model = PWAModel(@vanderpol_PWA_f,Partition,C,Bw,[0;0],eye(2));
            %
            % with function vanderpol_f defined below.
            %   function xp = vanderpol_f(x, u,sys)
            %     if nargin ~= 3
            %         error('Invalid input !');
            %     end
            % 
            %     if size(x,1)~=2
            %         error('Wrong size of x');
            %     end
            %     xp = zeros(size(x));
            % 
            %     xp(1,:) = x(1,:) + sys.param.tau *x(2,:);
            %     xp(2,:) = x(2,:)+ (-x(1,:) + (1-x(1,:).^2) .* x(2,:)) .* sys.param.tau + u;
            %     end
            % 
            % with Partition(i) defined below for i = 1
            % Partition(1).Polyhedron = ([combvec([-3 -2.85], [-3, -2.85])]')
            % Partition(1).Dynamics.A = [1 0; 0 1];
            % Partition(1).Dynamics.a = [1; 1]; 
            % Partition(1).Dynamics.B = [0;1];
            % Partition(1).Dynamics.C = [1 0];
            % Partition(1).Dynamics.D = 0;
            % Partition(1).Dynamics.Bw = [0;1];
            % Partition(1).Dynamics.dim = sysNonLin.dim;
            %
            % Copyright 2022 Birgit van Huijgevoort, b.c.v.huijgevoort@tue.nl

            obj.f=@PWA_f;
            obj.Partition = Partition;
            obj.C = C;
            obj.Bw = Bw;
            obj.mu = mu;
            obj.sigma = sigma;
            obj.dim = size(Bw, 1);
            obj.Np = size(Partition,2);
        end

        function x_n = f_det(obj,x,u)
            %F_DET computes the next states undisturbed by noise 
            % The next state is based on the deterministic/nominal
            % dynamics for a nonlinear model and is computed as
            % x(t+1) = f(x(t),u(t))
            %
            % x_next = f_det(obj,x_current,u_current) computes the next
            % state x_next based on x_current, u_current, and based on the
            % model obj. 
            x_n =   obj.f(x,u,obj);
        end

        function [x_n, varargout] = f_stoch(obj, x, u, Pr)
            %F_STOCH computes the next states disturbed by noise 
            % The next state is based on the stochastic
            % dynamics for a PWA model and is computed as
            % x(t+1) = A_i x(t) + B_i u(t) + Bw_i w(t) + a_i
            % where w(t) is sampled at random from it's specified
            % underlying distribution and i is the partition number
            %
            % x_next = f_stoch(obj, x_current, u_current, Partition number) computes the next
            % state x_next based on x_current, u_current, partititon number  and based on the
            % model obj. 
            
            w = mvnrnd(obj.mu, obj.sigma, 1)'; % Sample noise
              
            % Get corresponding matrices
            A_i = obj.Partition(Pr).Dynamics.A;
            B_i = obj.Partition(Pr).Dynamics.B;
            Bw_i = obj.Partition(Pr).Dynamics.Bw;
            a_i = obj.Partition(Pr).Dynamics.a;

            % Determine next state
            x_n = A_i*x + B_i*u + a_i + Bw_i*w;
            varargout{1} = w;
            
        end

        function xp = PWA_f(x, u, sys)
            xp = zeros(size(x));
    
            for l = 1:size(sys.Partition,2)
                ind = find(sys.Partition(l).Polyhedron.contains(x));
                xp(:,ind) = sys.Partition(l).Dynamics.A*x(:,ind)+ sys.Partition(l).Dynamics.B*u(ind);
            end
        end
        
        
    end
end
