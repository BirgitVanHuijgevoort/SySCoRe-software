classdef Policy
    %POLICY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        pol
        DFA
    end
    
    methods
        function obj = Policy(pol,DFA)
            %POLICY Construct an instance of this class
            %   Detailed explanation goes here
            obj.pol = pol;
            obj.DFA = DFA;
        end
        

    end
end

