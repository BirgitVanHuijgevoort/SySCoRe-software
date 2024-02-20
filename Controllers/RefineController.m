classdef RefineController
    %REFINECONTROLLER Constructs a refined controller that can be
    %implemented on a continuous-state model.
    %
    % Objects of this class define a controller that yields a certain satisfaction
    % probablity by applying a control policy to a system
    %
    % Copyright 2024 Birgit van Huijgevoort bhuijgevoort@mpi-sws.org

    properties
        Prob % satisfaction probability
        pol % abstract control policy
        sysAbs % abstract model
        simRel % simulation relation
        int_f % Interface function \in [0,1,2] ...
        % 0. (default) u=uhat 1. u=uhat+K(x-xhat) 2. u = uhat + Qxr + K(x-Pxr)
        FbMatrix % matrix K for interface = 1 or 2
        sys % Can be either of class LinModel or NonlinModel (if MOR is 1, this is a reduced-order model, if KKfilter is 1 this is filtered model)
        sys2  % If both MOR and KKfilter, sys is reduced-order model (MOR) and sys2 is KK filtered model.

        DFA % Deterministic finite automaton
        regions % regions labelled with atomic propositions used for the specification (and DFA)
        dim % (state) dimension of the (finite-state) system

        MOR = 0 % \in [0,1], logical value for model-order reducion (MOR), 0 for no MOR, 1 for MOR
        P % the projection matrix for model-order reduction x = Pxr
        Q % for MOR, Q-matrix in interface function 2
        F % for MOR, this matrix gives the coupling of disturbances wr = w+F*(x-P*xr);

        KKfilter = 0 % \in [0,1], logical value for KK filtering. 0 for no KK filtering, 1 for KK filtering used to obtain the model
        K % kalmann matrix K of filtered model
        Cobs % observation matrix in y=Cx (paper Maico)
    end

    methods

        function obj = RefineController(satProb,policy,sysAbs,simRel,sys,DFA,varargin)
            %REFINECONTROLLER Constructs a refined controller
            %   Load all data as described by the properties.
            %   Controller =
            %   RefineController(satProb,policy,sysAbs,simRel,sys,DFA)
            %   loads the controller. 
            %
            %   use Refine_controller(satProb,pol,sysAbs,rel,sysLTIr,DFA,{int_f1,
            %   int_f2}, {K1, K2}); to specify two simulation relations and interface functions
            %   (when MOR is applied).
            %   Here, int_f1 and K1 refine input uhat from sysAbs to input ur from
            %   sysLTIr. Similarly, int_f2 and K2 refine input ur from sysLTIr to u
            %   from sysLTIr.original (full-order model). Note that they always have
            %   to specified to increasing order wrt behavior, that is sysAbs <= sysLTIr <=
            %   sysLTI. 
            %   Example: see the BAS tutorial, instead of Controller = RefineController(satProb,pol,sysAbs,rel,sysLTIr,DFA,int_f,K);
            %   we can obtain the exact same controller using the
            %   following:
            %   int_f1 = 1; K1 = zeros(1,2);
            %   Controller = RefineController(satProb,pol,sysAbs,rel,sysLTIr,DFA,{int_f1,int_f}, {K1, K});
            %
            % when a KK filter is used, use the following syntax: 
            % Controller =
            %   RefineController(satProb,policy,sysAbs,simRel,sys,DFA, [],[], 'KKfilter', sysLTI_KF)
            %   with the last argument the KK filtered model

            obj.Prob = satProb;
            obj.pol = policy;
            obj.sysAbs = sysAbs;
            obj.simRel = simRel;
            obj.sys = sys; % Can be either of class LinModel or NonlinModel
            obj.DFA = DFA;
            obj.regions = sys.regions;


            % Use the (state) dimension of the abstract system if it is
            % given
            if ~isempty(sysAbs)
                obj.dim = size(obj.sysAbs.states,1);
            else
                obj.dim = 1;
            end

            % Load the feedback matrix for interface function u =
            % uhat+K(x-xhat). If it is not given set it to zero (hence we
            % have interface u = uhat).
            if nargin > 6
                obj.int_f = varargin{1};
                if sys.type == 'PWA'
                    obj.FbMatrix = {sys.Partition.Kf};
                else
                    if nargin > 7
                        obj.FbMatrix = varargin{2};
                    else
                        udim = size(obj.sys.B,2);
                        obj.FbMatrix = zeros(udim,obj.dim);
                    end
                    if nargin > 8
                        if varargin{3} == 'KKfilter'
                            obj.KKfilter = 1;
                            obj.sys2 = varargin{4};
                        end
                    end

                end
            else
                obj.int_f = 0;
                udim = size(obj.sys.B,2);
                obj.FbMatrix = zeros(udim,obj.dim);
            end
            if isempty(obj.FbMatrix)
                udim = size(obj.sys.B,2);
                obj.FbMatrix = zeros(udim,obj.dim);
            end
        end

        % Compute state evolution of actual system, that is
        % determine next concrete state and corresponding optimal concrete
        % input
        % option 1) by looking at the maximizer in R wrt the value function.
        function [xnext, qnext, u, varargout] = EvolveSys(obj, x, q, varargin)
            %EVOLVESYS performs a state evolution x(t+1) of the system
            %starting at (x(t),q(t)) and following the evolution described
            %by the object (see folder Models for deterministic and stochastic state evolution). 
            %
            % Outputs
            % -------
            % xnext is the next state x(t+1), qnext is the next DFA state
            % q(t+1), u is the applied input and varargout is used to give
            % the full-order state when model-order reduction is applied. 
            % 
            % Inputs
            % ------
            % obj is the controller object
            % x is the state at time step t
            % q is the DFA state at time step t
            % 
            % Options (varargin, varargout)
            % -------
            % when model-order reduction is applied (obj.MOR=1), varargin
            % is used to supply the full-order state. In this case
            % varargout is the next state of the full-order model. 
            %
            % when also KKfiltering is applied (obj.KKfilter=1),
            % varargin{1} = full-order state, varargin{2} is state of
            % filtered abstract model. varargout{1} = next state of
            % full-order original model, varargout{2} = next state of
            % filtered model. 

            if obj.MOR
                % (xfull, xfullnext) is full state, (x, xnext) is reduced state
                xfull = varargin{1};
            end

            if obj.MOR && obj.KKfilter
                xbar = varargin{2};
            end

            if ~obj.MOR && obj.KKfilter
                xbar = varargin{1};
            end

            % Find abstract state and input (xrhat, urhat for MOR)
            [xhat, uhat] = obj.determAbsXU(x, q);

            % Get initial control input
            if obj.sys.type == 'PWA'
                % Determine partition number
                temp = ceil((x-(min(obj.sys.X.V)'))./((max(obj.sys.X.V)-min(obj.sys.X.V))').*(obj.sys.N(1:obj.dim)'-1)); %Coordinate of the partition ([row;column] in 2D)
                sz = obj.sys.N(1:obj.sys.dim)-1;
                M = [];
                for i = 1:size(temp,1)
                    M = [M {temp(i,:)}];
                end
                Pr = sub2ind(sz, M{:}); % partition number
                
                % Get initial control input
                u = obj.InterfaceFunction(uhat, x, xhat, Pr);
            else
                u = obj.InterfaceFunction(uhat, x, xhat);
            end

            % Compute next state
            if obj.MOR && ~obj.KKfilter
                % Full-order model
                ufull = uhat + obj.Q*x + obj.FbMatrix * (xfull - obj.P*x);
                [xfullnext, wfull] = obj.sys.original.f_stoch(xfull, ufull);
                varargout{1} = xfullnext;

                % Reduced-order model
                wr = wfull + obj.F*(xfull-obj.P*x);
                xnext = obj.sys.f_stoch(x, u, wr); % xrnext
            elseif obj.MOR && obj.KKfilter
                % Full-order model
                ufull = uhat + obj.Q*x + obj.FbMatrix * (xfull - obj.P*x);
                [xfullnext, wfull] = obj.sys.original.f_stoch(xfull, ufull);
                varargout{1} = xfullnext;

                % Reduced-order model
                wr = wfull + obj.F*(xfull-obj.P*x);
                xnext = obj.sys.f_stoch(x, u, wr); % xrnext

                % KK-filtered model 
                ubar = u;
                Obsnext = obj.Cobs*xfullnext; % y_t = Cx_t in paper Maico
                vnext = Obsnext-obj.Cobs*obj.sys2.A*xbar-obj.Cobs*obj.sys2.B*ubar; 
                %%%% ----- Should be replaced by something like:
                % %%%% ------ xbarnext = obj.sys2.f_stoch(xbar, ubar, vnext); --- %%
                xbarnext = obj.sys2.A*xbar+obj.sys2.B*ubar+obj.K*vnext;

                varargout{2} = xbarnext;
            elseif ~obj.MOR && obj.KKfilter
                warning('Control refinement with only knowledge filtering (and no MOR) has not been verified yet.')
                xnext = obj.sys.f_stoch(x, u);

                % KK-filtered model 
                ubar = u;
                Obsnext = obj.Cobs*xnext; % y_t = Cx_t in paper Maico
                vnext = Obsnext-obj.Cobs*obj.sys2.A*xbar-obj.Cobs*obj.sys2.B*ubar; 
                %%%% ----- Should be replaced by something like:
                % %%%% ------ xbarnext = obj.sys2.f_stoch(xbar, ubar, vnext); --- %%
                xbarnext = obj.sys2.A*xbar+obj.sys2.B*ubar+obj.K*vnext;
                xbarnext = obj.sys2.f_stoch(xbar, ubar, vnext); % check if obj.K = obj.sys2.Bw!?!
                varargout{1} = xbarnext;
            elseif obj.sys.type == 'PWA'
                xnext = obj.sys.f_stoch(x, u, Pr);
            else
                xnext = obj.sys.f_stoch(x, u);
            end


            % Check if xnext in sys.X
            if obj.MOR
                if ~(obj.sys.original.X.contains(xfullnext))
                    %error("xnext is outside of sys.X")
                     warning("x(t+1) is outside of the state space, so I changed it to the closest state inside the state space")
                    % find closest state inside state space
                    xnext = obj.sys.X.distance(xnext).y;
                end
            else
                if ~(obj.sys.X.contains(xnext))
                    warning("x(t+1) is outside of the state space, so I changed it to the closest state inside the state space")
                    % find closest state inside state space
                    xnext = obj.sys.X.distance(xnext).y;
                end
            end

            % Determine output of actual system and keep track of DFA
            if obj.MOR
                ynext = obj.getSysOutput(xfullnext);
                qnext = updateDFA(obj.DFA, obj.sys.original.regions, obj.sys.original.AP, q, ynext);
            else
                % PWA and non-reduced
                ynext = obj.getSysOutput(xnext);
                qnext = updateDFA(obj.DFA, obj.regions, obj.sys.AP, q, ynext);
            end
        end

     
        function q_0 = initDFA(obj, x_0)
            %INITDFA computes the initial state of the DFA q_0 by updating it
            %instantely based on the label of the initial output. The initial output is computed 
            % based on the initial state x_0.

            % Initial state of DFA
            q_init = obj.DFA.S0;

            % Get initial system output (should be independent of u)
            y_init = obj.getSysOutput(x_0);

            % Update DFA state (the DFA is always updated instantly)
            if obj.MOR
                q_0 = updateDFA(obj.DFA, obj.sys.original.regions, obj.sys.original.AP, q_init, y_init);
            else
                q_0 = updateDFA(obj.DFA, obj.regions, obj.sys.AP, q_init, y_init);
            end
        end

        function [xhat, uhat] = determAbsXU(obj, x, q)
            %DETERMABSXU  Determines the abstract state xhat and input uhat
            % corresponding to the concrete state (x,q).

            % Find next abstract state by looking at the maximizer in R wrt the value
            % function 
            indexing = 1:length(obj.sysAbs.states);
            if obj.MOR
                inR = obj.simRel.R{2}.inR(x, obj.sysAbs.states);
            else
                inR = obj.simRel.inR(x, obj.sysAbs.states);
            end
            assert(sum(inR) > 0, "No valid state in relation found.")
            indices_valid = indexing(inR);
            if size(obj.Prob,1) > 1
                [~, index_aux] = max(obj.Prob(q, inR));
            else
                [~, index_aux] = max(obj.Prob(inR));
            end
            j = indices_valid(index_aux); % Find maximizing index of abstract state
            xhat = obj.sysAbs.states(:, j);
    
            % Get abstract input from policy
            uhat = obj.evalPol(j, q);
        end

        function y = getSysOutput(obj, x)
            %GETSYSOUTPUT Determines the output of the actual system based
            %on the state x and the system stored in the object controller.

            if obj.MOR
                if isa(obj.sys.original, 'LinModel')
                    assert(sum(obj.sys.original.D, 'all') == 0, "Systems with control inputs in the output equation aren't supported.")
                    y = obj.sys.original.C * x; % + obj.sys.D * u;
                elseif isa(obj.sys.original, 'NonlinModel')
                    error("To be implemented.")
                else
                    error("Model type %s not supported with MOR.", class(obj.sys.original))
                end
            else
                if isa(obj.sys, 'LinModel')
                    assert(sum(obj.sys.D, 'all') == 0, "Systems with control inputs in the output equation aren't supported.")
                    y = obj.sys.C * x; % + obj.sys.D * u;
                elseif isa(obj.sys, 'NonlinModel')
                    error("To be implemented.")
                elseif isa(obj.sys, 'PWAModel')
                    y = obj.sys.C * x;
                else
                    error("Model type %s not supported.", class(obj.sys))
                end
            end
        end

        function uhat = evalPol(obj, j, varargin)
            %EVALPOL evaluates the policy corresponding to the j-th input.
            % varargin is used to supply the DFA state q. 

            if numel(size(obj.pol)) > 2
                % Switch between pages according to current DFA state
                q = varargin{:};
                uhat = obj.pol(:, j, q);
            else
                % Only initial policy
                uhat = obj.pol(:, j);
            end
        end

        function u = InterfaceFunction(obj, uhat, x, xhat, varargin)
            % INTERFACEFUNCTION determine the next input u based on
            % interface function u = uhat + K(x-xhat). 

            if obj.MOR
                if isa(obj.FbMatrix, 'cell')
                    u = uhat + obj.FbMatrix{1} * (x-xhat);
                else
                    % FbMatrix not given, so no feedback-term
                    u = uhat;
                end
            elseif obj.sys.type == 'PWA'
                Pr = varargin{1}; % partition number
                FbMatrix = obj.FbMatrix{Pr};
                u = uhat + FbMatrix*(x-xhat);
            else
                u = uhat + obj.FbMatrix * (x-xhat);
            end
        end
    end
end
