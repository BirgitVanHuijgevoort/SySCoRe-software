%% Test Class Definition
classdef tPolytope < matlab.unittest.TestCase

    properties
        sysLTI
        sysAbs
        DFA
        rel
    end

    methods  (TestClassSetup)
        function createArguments(testCase)
            Box2d = Polyhedron(combvec([-1,1], [-1,1])');
            Box3d= Polyhedron(combvec([-1,1], [-1,1], [-1,1])');

        end
    end
    %% Test Method Block
    methods (Test)
        %% Test basic functinality of EReachTime
        function testBox2d(testCase)
            Box2d = Polyhedron(combvec([-1,1], [-1,1])');
            
        end

        function testBox3d(testCase)
            Box3d= Polyhedron(combvec([-1,1], [-1,1], [-1,1])');
        end


        function test_IncreaseDecreasePolytope(testCase)
            Box2d = Polyhedron(combvec([-1,1], [-1,1])');
            Box2d.computeHRep()
    
            [Box2d_big,Box2d_small] = IncreaseDecreasePolytope(Box2d, .7);

            Box2d_big.computeVRep()
            Box2d_small.computeVRep()
            V_big = combvec([-1,1], [-1,1])'*1.7;
            testCase.verifyEqual(sortrows(Box2d_big.V),sortrows(V_big), "AbsTol", 0.01)

            V_small = combvec([-1,1], [-1,1])'*0.3;
            testCase.verifyEqual(sortrows(Box2d_small.V),sortrows(V_small),"AbsTol", 0.01)

        end

    end
end