classdef SimRel
    %SIMREL Class of simulation relation Rel
    %   Objects of this class define a simulation relation Rel = {(xhat,x) | ||x-xhat||_R \leq epsilon}
    %   between a continuous-state model with states x and a finite-state
    %   model states xhat. The output deviation is denoted by epsilon and R
    %   is a positive definite weighting matrix. 
    %
    %   This class can also be used to define a simulation relation 
    %   Rel_r = {(xr,x) | ||x-Pxr||_Rr \leq epsilon_r}
    %   between a continuous-state model with states x and a reduced-order
    %   model with states xr. 
    
    properties
        epsilon % output deviation epsilon 
        delta   % probability deviation delta \in [0,1]
        NonDetLabels    % non-deterministic labels of the finite-state model
        R       % weighting matrix R
        P       % projection matrix P for a simulation relation between a model and its reduced-order approximate
        states % in case the simrel is composed of two simrels these are the states of the model that is ommitted
    end
    
    methods
        function obj = SimRel(epsilon,delta,R, varargin)
            %SIMREL Construct an instance of this class
            %   Load all values epsilon, delta and R as described above. 
            %   Rel = SimRel(epsilon,delta,R) loads the simulation relation 
            %   
            %   Examples
            %   A simple two-dimensional simulation relation
            %   epsilon = 0.1;
            %   delta = 0.0016;
            %   R = eye(2);
            %   Rel = SimRel(epsilon,delta,R)
            %
            %   A two-dimensional simulation relation for model order
            %   reduction
            %   epsilon = 0.1;
            %   delta = 0.0016;
            %   R = eye(2);
            %   P = [0.5, 0; 0; 8];
            %   Rel = SimRel(epsilon,delta,R,P);
                
            obj.epsilon = epsilon; 
            obj.delta = delta;
            obj.R = R;
            if nargin >=4
                obj.P = varargin{1};
            else
                obj.P = eye(size(R));   
            end
            if nargin >=5
                obj.states = varargin{2};
            end
        end
        
        function bool = inR(obj,x,xh,varargin)
            %INR Check whether two states belong to the same simulation
            %relation
            % 
            %   Example
            %   epsilon = 0.1;
            %   delta = 0.0016;
            %   R = eye(2);
            %   Rel = SimRel(epsilon,delta,R);
            %   inR = Rel.inR(0.4,0.5);
            if isa(obj.R, 'double')
                len = length(obj.R);
                if nargin==3
                    % Make more efficient if x is single point
                    if size(x, 2) == 1
                        x_c = x - obj.P * xh;
                    else
    %                   x = x.*ones(size(xh));
                        x_c = [eye(size(x,1)), -  obj.P]*combvec(x,xh);
                    end
                elseif strcmp(varargin{1} , '1-to-1')
                    x_c = x-obj.P*xh;
                end
         
                bool = ((ones(1,len)*((obj.R^.5*x_c).^2)).^.5)<=obj.epsilon+eps; %added machine precision
                bool = reshape(bool,size(x,2),size(xh,2));
            
            elseif iscell(obj.R) && isa(obj.R{1},'SimRel')
                if size(x,1) < size(xh,1)

                    bool_states = inR(obj.R{1},x,obj.states); % check which of the states are related
                    states_rel = obj.states(:,bool_states);
                
                    bool_states = inR(obj.R{2},states_rel,xh); % check which of the states are related
                    bool = any(bool_states,1);
                
                elseif size(x,1) > size(xh,1)
                    bool_states = inR(obj.R{2},obj.states,xh); % check which of the states are related
                    states_rel = obj.states(:,bool_states);
                
                    bool_states = inR(obj.R{1},x,states_rel); % check which of the states are related
                    bool = any(bool_states,2)';                  
                elseif size(x,1) == size(xh,1)
                    % determine whether we should check R{1} or R{2}
                    if size(x,1) == size(obj.R{1}.R,1)
                        bool = inR(obj.R{1},x,xh);
                    elseif size(x,1) == size(obj.R{2}.R,1)
                        bool = inR(obj.R{2},x,xh);
                    else
                        error('Check if dimensions of states and simulation relation are the same')
                    end

                else
                    error('Check dimensions of states supplied to function inR')
                end                    
            end
            
            end

        
        function obj_comb = Combine(obj,obj2, states)
            %COMBINE Combines two simulation relations into one.
            epsilonc = obj.epsilon+obj2.epsilon; 
            deltac = obj.delta+obj2.delta;
            Dc = {obj,obj2};
            
            obj_comb = SimRel(epsilonc,deltac,Dc,[],states);
            
        end
        
        
        
    end
end

