classdef LinModel
% LINMODEL Class of LTI systems with noise on the transitions. 
% Objects of this class define a model with dynamics
%  x(t+1) = A x(t) + B u(t) + Bw w(t)
%  y(t)   = C x(t)
% with w(t) ~ N(mu, sigma) where N is a normal distribution with mean mu
% and sigma the noise variance
%
% Examples 
% A simple 2 dimensional model with one control input and one
% output:
%   A = [1 0; 0 1]; 
%   B = [1; 0]; 
%   C = [1, 1]; 
%   D = [0]; 
%   Bw = [0; 1]; 
%   model = LinModel(A, B, C, D, Bw, [0; 0], eye(2))
%
% See also: LinModel (Constructor)
% Copyright 2021 Sofie Haesaert s.haesaert@tue.nl
    
    properties
        type = "LTI"; % Type of model
        A % A matrix of LTI dynamics
        B % B Matrix of LTI dynamics
        C  % C matrix of LTI matrix towards the output dimensions
        D % D matrix of LTI matrix (Should be equal to zero)
        Bw % Bw matrix for noise inputs of LTI dynamics
        mu % noise mean vector
        sigma % noise covariance matrix
        wsupport % noise support
        dim  % the dimension of the state space
        
        X % bounded state space  
        U % bounded input space 
        regions % regions labelled with atomic propositions
        AP % atomic propositions
        P % the probability transition matrix.  
    end
    
    methods
        function obj = LinModel(A, B, C, D, Bw, varargin)
            %LINMODEL Construct an instance of this class
            %   Load all values A, B, C, D, Bw, mu, sigma as described in the
            %   class file. 
            %   model = LinModel(A, B, C, D, Bw, mu, sigma) Loads model dynamics.
            % 
            % 
            % Examples
            % A simple 2 dimensional model with one control input and one output:
            %   A = [1 0; 0 1];
            %   B = [1; 0];
            %   C = [1, 1];
            %   D = [0];
            %   Bw = [0; 1];
            %   model = LinModel(A, B, C, D, Bw, [0; 0], eye(2))

            obj.A = A;
            obj.B = B;
            obj.D = D;
            obj.C = C;
            obj.Bw = Bw;
            obj.dim = size(Bw, 1);

            % Allow for different noise distributions
            if nargin <= 4
                error("LinModel requires at least 5 arguments.")
            elseif (nargin == 6 && numel(varargin{1}) >= 2)
                % Bounded uniform distribution
                obj = obj.setWSupport(varargin{1});
            elseif nargin == 7
                % Unbounded Gaussian distribution
                obj.mu = varargin{1};
                obj.sigma = varargin{2};
                obj = obj.setWSupport([-inf, inf]);
            elseif nargin == 8
                % Bounded Gaussian distribution
                obj.mu = varargin{1};
                obj.sigma = varargin{2};
                obj = obj.setWSupport(varargin{3});
            else
                error("Unsupported arguments for calling LinModel.")
            end
            
        end

        function x_n = f_det(obj, x, u)
            %F_DET computes the next states undisturbed by noise 
            % The next state is based on the deterministic/nominal
            % dynamics for a linear model and is computed as
            % x(t+1) = A x(t) + B u(t)
            %
            % x_next = f_det(obj,x_current,u_current) computes the next
            % state x_next based on x_current, u_current, and based on the
            % model obj. 

            x_n = obj.A*x + obj.B*u;
        end

        function x_n = f_stoch(obj, x, u)
            %F_STOCH computes the next states disturbed by noise 
            % The next state is based on the stochastic
            % dynamics for a linear model and is computed as
            % x(t+1) = A x(t) + B u(t) + Bw w(t)
            % where w(t) is sampled at random from it's specified
            % underlying distribution
            %
            % x_next = f_stoch(obj, x_current, u_current) computes the next
            % state x_next based on x_current, u_current, and based on the
            % model obj. 

            w = mvnrnd(obj.mu, obj.sigma, 1)'; % Sample noise
            x_n = obj.A*x + obj.B*u + obj.Bw*w;
        end
    end

    methods (Hidden)
        function obj = setWSupport(obj, wsupport)
            assert(size(wsupport, 2) == 2, "Invalid noise support.")
            if (size(wsupport, 1) == 1 && size(wsupport, 1) ~= obj.dim)
                % Copy support for all dimensions
                obj.wsupport = repmat(wsupport, obj.dim, 1);
            elseif size(wsupport, 1) == obj.dim
                obj.wsupport = wsupport;
            else
                error("Invalid noise support.")
            end
        end
    end
end
