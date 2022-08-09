classdef TensorTransitionProbability
    %TENSORTRANSITIONPROBABILITY This class can be used for a unified way of
    %multiplying the value function with the transition probability for
    %dimensions up to 2. For higher dimensions the tensortoolbox needs to
    %be used. 
    
    properties
        P_det
        dim
        l
        Pi
    end
    
    methods
        function obj = TensorTransitionProbability(l, P_det, varargin)
            %TENSORTRANSITIONPROBABILITY Construct an instance of this class
            %   
            obj.P_det = P_det;
            obj.dim = length(l);
            obj.l=l;
            if length(varargin) ~= obj.dim
                error('incorrect number of inputs')
            end
            
            obj.Pi = cell(obj.dim,1);
            for i = 1:obj.dim
                obj.Pi{i} = varargin{i};
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
            % If input is V e R^1xn, than output is V_n e R^axn 
            % with a the number of finite actions 

            if obj.dim ~= 2
                error("Dimension different from 2 has not yet been implemented set 'TensorToolbox' to   'tensorlab' or 'tensortoolbox'  ")
            end

            VI = vec(obj.Pi{1}'*reshape(V,obj.l)*obj.Pi{2}); 
            
            V_n = reshape(VI'*obj.P_det,length(VI),[]);  

        end
    end
end
