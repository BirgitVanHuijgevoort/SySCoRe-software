classdef TransitionProbability
    %TRANSITIONPROBABILITY This class can be used for a unified way of
    %multiplying the value function with the transition probability. 
    
    properties
        Prob
    end
    
    methods
        function obj = TransitionProbability(Prob_value)
            %TRANSITIONPROBABILITY Construct an instance of this class
            %   Detailed explanation goes here
            obj.Prob =  reshape(permute(Prob_value,[2 1 3]),size(Prob_value, 1),[]);

        end
        
        function V_n = mtimes(V, obj)
            %MTIMES multiply value function with probability transition 
            % If input is V e R^1xn, than output is V_n e R^axn 
            % with a the number of finite actions 
            V_n = reshape(V*obj.Prob,length(V),[]); 

        end
    end
end

