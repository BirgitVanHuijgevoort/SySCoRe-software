%% Test Class Definition
classdef tMDP < matlab.unittest.TestCase
       
    properties
        P
        MDP
    end
    
    methods  (TestClassSetup)
        function createArguments(testCase)
            GW = createGridWorld(2,2);
            Ps = permute(GW.T, [2,1,3]); 
            testCase.P = reshape(Ps, [4, 4,4]);
            hx1 = [1:2];
            hx2 = [1:2];   
            testCase.MDP = MDP_model(testCase.P,{hx1,hx2});
        end
    end

    %% Test Method Block
    methods (Test)

        %% Test basic functinality to load values
        function testBasicLoad(testCase)
            %% Exercise function under test
            % act = the value from the function under test
            createArguments(testCase)
            testCase.verifyEqual(testCase.MDP.hx{1},[1,2]);
            testCase.verifyEqual(testCase.MDP.hx{2},[1,2]);
        end

    % TODO: add tests for the MDP_model class. 


    end
end