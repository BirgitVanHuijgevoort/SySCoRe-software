%% Test Class Definition
classdef tSimRel < matlab.unittest.TestCase

    %% Test Method Block
    methods (Test)

        %% Test basic functinality to load values
        function testBasicLoad(testCase)
            %% Exercise function under test
            % act = the value from the function under test
            rel_act = SimRel(0.1,.2,eye(2));
            %% Verify using test qualification
            testCase.verifyEqual(rel_act.epsilon,0.1);
            testCase.verifyEqual(rel_act.delta,0.2);
            testCase.verifyEqual(rel_act.R,eye(2));


        end
        %% Test basic functinality to load values
        function testAdvanceLoad(testCase)
            %% Exercise function under test
            % act = the value from the function under test
            rel_act = SimRel(0.1,.2,eye(2), 2*eye(2));
            %% Verify using test qualification
            testCase.verifyEqual(rel_act.P,2*eye(2));

        end


        %% Test basic functinality to load values
        function testInR(testCase)
            %% Exercise function under test
            % act = the value from the function under test
            rel_act = SimRel(0.1,.2,eye(2));
            act = rel_act.inR([1;1], [1;1.02]);

            %% Verify using test qualification

            testCase.verifyEqual(act,true);

        end
         function testInRfalse(testCase)
            %% Exercise function under test
            % act = the value from the function under test
            rel_act = SimRel(0.1,.2,eye(2));
            act = rel_act.inR([1;1], [1;1.12]);

            %% Verify using test qualification

            testCase.verifyEqual(act,false);

        end
    end
end