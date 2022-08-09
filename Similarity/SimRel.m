classdef SimRel
    %SIMREL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        epsilon
        delta
        NonDetLabels
        R
        P
        states % in case the simrel is composed of two simrels tese are the states of the model that is ommitted
    end
    
    methods
        function obj = SimRel(epsilon,delta,D_m, varargin)
            %SIMREL Construct an instance of this class
            %   Detailed explanation goes here
            obj.epsilon = epsilon; 
            obj.delta = delta;
            obj.R = D_m;
            if nargin >=4
                obj.P = varargin{1};
            else
                obj.P = eye(size(D_m));

                
            end
            if nargin >=5
                obj.states = varargin{2};
            end
            
%             disp(['Set simulation relation with epsilon = ',...
%                 num2str( obj.epsilon), ', delta = ', num2str(obj.delta)])
            
        end
        
        function bool = inR(obj,x,xh,varargin)
            %INR Check whether two states belong to the same simulation
            %relaiton
            %   Detailed explanation goes here
            if isa(obj.R, 'double')
                len = length(obj.R);
                if nargin==3
%                   x = x.*ones(size(xh));
                    x_c = [eye(size(x,1)), -  obj.P]*combvec(x,xh);
                
                elseif strcmp(varargin{1} , '1-to-1')
                    x_c = x-obj.P*xh;
                end
         
                bool = ((ones(1,len)*((obj.R^.5*x_c).^2)).^.5)<=obj.epsilon+eps; %added machine precision
                bool = reshape(bool,size(x,2),size(xh,2));
            
            elseif iscell(obj.R) && isa(obj.R{1},'SimRel')
                if size(x,2)<=size(xh,2)
                bool_states = inR(obj.R{1},x,obj.states); % check which of the states are related
                states_rel = obj.states(:,bool_states);
                
                bool_states = inR(obj.R{2},states_rel,xh); % check which of the states are related
                bool = any(bool_states,1);
                
                else
                bool_states = inR(obj.R{2},obj.states,xh); % check which of the states are related
                states_rel = obj.states(:,bool_states);
                
                bool_states = inR(obj.R{1},x,states_rel); % check which of the states are related
                bool = any(bool_states,2)';
                    
                    
                end
                    
                    
                
            end
            
                
            
            
            end

        
        function obj_comb = Combine(obj,obj2, states)
            %COMBINE Combine two simulation relations
            epsilonc = obj.epsilon+obj2.epsilon; 
            deltac = obj.delta+obj2.delta;
            Dc = {obj,obj2};
            
            obj_comb = SimRel(epsilonc,deltac,Dc,[],states);
            
        end
        
        
        
    end
end

