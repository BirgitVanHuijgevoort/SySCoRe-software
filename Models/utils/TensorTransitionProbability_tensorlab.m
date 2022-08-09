classdef TensorTransitionProbability_tensorlab
    %TENSORTRANSITIONPROBABILITY This class can be used for a unified way of
    %multiplying the value function with the transition probability. 
    % This 2nd version is created to try out the tensor toolbox 
    % for this it is needed to install the software of tensortoolbox.org by
    % adding it to your path. 
    % TODO: Add input arguments here. 

    properties
        P_det
        dim
        l
        Pi
        lact
    end
    
    methods
        function obj = TensorTransitionProbability_tensorlab(l, P_det, varargin)
            %TENSORTRANSITIONPROBABILITY Construct an instance of this class
            %   
            obj.P_det = P_det;
            obj.dim = length(l);
            obj.l=l;
            obj.lact = size(P_det,2)/prod(obj.l);
            if length(varargin) ~= obj.dim
                error('incorrect number of inputs')
            end
            
            obj.Pi = cell(obj.dim,1);
            for i = 1:obj.dim
                obj.Pi{i} = varargin{i}';  % transpose only needed for tensorlab
            end
            
            
        end
        function V_n = dettimes(V, obj)
            %DETTIMES multiply value function with deterministic transitions 
            % If input is V e R^1xn, than output is V_n e R^axn 
            % with a the number of finite actions 
            V_n = reshape(V*obj.P_det,length(V),[]); 

        end
        
        
        function V_n = mtimes(V, obj)
            %MTIMES multiply value function with probability transition 
            % If input is V e R^1xn, than output is V_n e R^lactxn 
            % with a the number of finite actions 


            % operations with tensorlab
            VI = tmprod(reshape(V,obj.l), obj.Pi,1:obj.dim);  
            V_n = reshape(VI(:)'*obj.P_det,prod(obj.l),[]); 

            

        end
    end
end
