%% Test Class Definition
%   TODO: 
%       - Add different test for different types of specifications including
%            next, or, and, ...
classdef tSpecifications < matlab.unittest.TestCase



    %% Test Method Block
    methods (Test)

        %% Test basic functinality to load values
        function testReachAvoid(testCase)
            formula = '(p1 U p2)';  % p1 = safe region, p2 = target region
            [DFA] = TranslateSpec(formula,{'p1','p2'});
            testCase.verifyEqual(DFA.S, [1,2,3]);
            testCase.verifyEqual(DFA.S0, 2);
            testCase.verifyEqual(DFA.F, 1);
            testCase.verifyEqual(DFA.sink, 3);

            trans = [ 0     0     0     0;...
                      3     1     2     1;...
                      3     3     3     3];
            
            testCase.verifyEqual(DFA.trans, trans);


        end

    end

end
