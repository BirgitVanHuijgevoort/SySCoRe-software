%% Test Class Definition
classdef tDependencies < matlab.unittest.TestCase

    properties
        sysLTI
        sysAbs
        DFA
        rel
    end

   
    
    %% Test Method Block
    methods (Test)
        %% Test basic functinality of EReachTime
        function testMosek(testCase)
            assert(exist('mosek.lic')==2,...
                "No Mosek license is found on the Matlab path. Please make sure you have a complete licensed installation of Mosek")
        end

        function testMPT(testCase)
           assert(exist("mpt_init")==2, " MPT  cannot be found on the search path")
           assert(exist("mpt_computeTrajectory",'file')==0, "Old version of MPT is used. Please use MPT3. ")
           assert( exist("mptopt",'file')==2, " MPT doesn't  seem to be fully installed, please run mpt_init.")
        end

        function testTensorToolbox(testCase)
        assert(strcmp(class(tensor), 'tensor'), "Check whether the tensor class is known")
        end
    end
end