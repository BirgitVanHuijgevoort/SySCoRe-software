classdef NonlinModel
    %NONLINMODEL Summary of this class goes here
    %   Nonlinear model with dynamics
    %   x(t+1) = f(x(t), u(t)) + Bw w(t) with w(t)\sim N(mu, \sigma)
    %         f = the deterministic function for the state update
    %           f:X\times U -> X
    %         C = C matrix to map from state to to output space
    %         Bw = B matrix for noise
    %         mu = noise mean
    %         sigma = noise variance
    %         dim  = the dimension of the state space
    %
    %   
    % Examples 
    % A simple 2 dimensional model with one control input and one
    % output:
    %   C = [1, 1]; 
    %   Bw = [0; 1]; 
    %   model = NonlinModel(@vanderpol_f,C,Bw,[0;0], eye(2),param)
    %
    %     x = sym('x', [sysNonLin.dim 1]);
    %     f1sym = x(1)+x(2)*sysNonLin.param.tau;
    %     f2sym = x(2)+(-x(1)+(1-x(1)^2)*x(2))*sysNonLin.param.tau;
    %     sysNonLin.fsym = [f1sym; f2sym];
    % 
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
    % See also: LinModel
    % 
    % Copyright 2021,   Birgit van Huijgevoort, Sofie Haesaert
    % (s.haesaert@tue.nl).


    properties
        type = "NonLin";
        f % the deterministic function for the state update
        fsym % same as f but symbolic
        C % C matrix to map from state to to output space
        Bw % B matrix for noise
        mu % noise mean
        sigma % noise variance
        dim  % the dimension of the state space

        param % any additional matrices to be shipped
        
        X % bounded state space  
        U % bounded input space
        regions % regions labelled with atomic propositions
        AP % atomic propositions
        P % the probability transition matrix.  
        
        
    end
    
    methods
        function obj = NonlinModel(f,C,Bw,mu,sigma,param)
        %   Nonlinear model with dynamics
        %   x(t+1) = f(x(t), u(t)) + Bw w(t) with w(t)\sim N(mu, \sigma)
        %         f = the deterministic function for the state update
        %           f:X\times U -> X
        %         C = C matrix to map from state to to output space
        %         Bw = B matrix for noise
        %         mu = noise mean
        %         sigma = noise variance
        %         dim  = the dimension of the state space
        %
        %   
        % Examples 
        % A simple 2 dimensional model with one control input and one
        % output:
        %   C = [1, 1]; 
        %   Bw = [0; 1]; 
        %   model = NonlinModel(@vanderpol_f,C,Bw,[0;0], eye(2),param)
        %
        %     x = sym('x', [sysNonLin.dim 1]);
        %     f1sym = x(1)+x(2)*sysNonLin.param.tau;
        %     f2sym = x(2)+(-x(1)+(1-x(1)^2)*x(2))*sysNonLin.param.tau;
        %     sysNonLin.fsym = [f1sym; f2sym];
        % 
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
        % See also: LinModel
        % 
        % Copyright 2021,   Birgit van Huijgevoort, Sofie Haesaert (s.haesaert@tue.nl), 

            obj.f = f;
            obj.C =  C;
            obj.Bw = Bw;
            obj.mu = mu;
            obj.sigma = sigma;
            obj.dim = size(Bw, 1);
            obj.param = param;
        end
        
                
        function x_n = f_det(obj,x,u)
            %F compute the next states undisturbed by noise based on the deterministic/nominal dynamics
            x_n =   obj.f(x,u,obj);
        end
        
    end
end
