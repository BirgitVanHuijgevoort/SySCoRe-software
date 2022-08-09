%% Test Class Definition
classdef tTensor < matlab.unittest.TestCase

    properties
        sysLTI
        sysAbs
        DFA
        rel
    end

    methods  (TestClassSetup)
        function createArguments(testCase)
            % Define model

            %% Specify system parameters and regions
            % LTI systems of the form
            % x(t+1) = Ax(t) + Bu(t) + Bw w(t)
            % y(t) = Cx(t) + Du(t)
            % Define all system parameters (incl Bw) into a struct
             A = 0.9*eye(2);
             B = 0.7*eye(2);
             C = eye(2);
             D = zeros(2);
             Bw = eye(2);
             dim = length(A);
            
            % Specify mean and variance of disturbance w(t) --> should be 2x1 and 2x2!
             mu = 0; % mean of disturbance
             sigma = 1; % variance of disturbance
             
            testCase.sysLTI = LinModel(A,B,C,D,Bw,mu,sigma);
             
            % Bounds on state space 
                x1l = -10;   % Lowerbound x1
                x1u = 10;   % Upperbound x1
                x2l = -10;   % Lowerbound x2
                x2u = 10;   % Upperbound x2
            testCase.sysLTI.X = Polyhedron(combvec([x1l,x1u],[x2l,x2u])');
            
            % Bounds on  and input space
                ul = [-1;-1];   % Lowerbound input u
                uu = [1;1];     % Upperbound input u
            testCase.sysLTI.U = Polyhedron(combvec([ul(1),uu(1)],[ul(2),uu(2)])');
            
            % Specify regions for the specification
                p1x = [4 4 10 10 4];    % x1-coordinates
                p1y = [0 -4 -4 0 0];    % x2-coordinates
                p1 = [p1x; p1y];         % parking region
            P1 = Polyhedron(p1');
            
                p2x = [4 4 10 10 4];    % x1-coordinates
                p2y = [0 4 4 0 0];      % x2-coordinates
                p2 = [p2x; p2y];        % avoid region
            P2 = Polyhedron(p2');
            
            testCase.sysLTI.regions = [P1;P2]; % regions that get specific atomic propositions
            testCase.sysLTI.AP = {'p1', 'p2'}; % with the corresponding atomic propositions
                        
            
            
            %% Synthesize scLTL formula (or input DFA yourself)
            %%% use LTL2BA and check if determinstic and accepting state with loop with 1.
            % input: (sc)LTL formula and atomic propositions (see readme in folder
            % LTL2BA)
            % output: struct DFA containing (among other) the transitions
            
            formula = '(!p2 U p1)';  % p1 = parking, p2 = avoid region
            % formula should use atomic propositions in sysLTI.AP. 
            
            % Make sure your current folder is the main SySCoRe folder
            [testCase.DFA] = TranslateSpec(formula,testCase.sysLTI.AP);
            
            
            %% Construct abstract model by gridding it
            disp('start gridding');tic
            % input: sysLTI, sigma, space bounds for input (controller) and state space
            
            % Specify division of input space for actuation and feedback
            ula = 1*ul;   % part of input for actuation (lowerbound)
            uua = 1*uu;
            ulf = ul-ula;   % part of input for feedback (lowerbound)
            uuf = uu-uua;
            
            lu = 3;
            uhat = combvec(linspace(ula(1),uua(1),lu),linspace(ula(2),uua(2),lu));
            
            
            l = [300, 300];  % number of grid cells in x1- and x2-direction
            tol=10^-6;
            % sysAbs = GridSpace_nd(sysLTI,uhat,l,tol);
            testCase.sysAbs = Gridding(testCase.sysLTI,uhat,l,tol, 'TensorComputation', 'true');
            
            % Save some extra system parameters into struct
            testCase.sysAbs.orig = testCase.sysLTI;
            
            label = zeros(1,prod(l));
            [label(1:prod(l))] = deal(1);
            inP1 = inpolygon(testCase.sysAbs.states(1,:),testCase.sysAbs.states(2,:),p1x,p1y);
            inP2 = inpolygon(testCase.sysAbs.states(1,:),testCase.sysAbs.states(2,:),p2x,p2y);
            [label(inP1)] = deal(3);
            [label(inP2)] = deal(2);
            testCase.sysAbs.labels = label;

            epsilon = 1.005;    % should be larger than vector beta!
            delta = 0.0163;  
            D_m=eye(2);
            testCase.rel = SimRel(epsilon,delta,D_m);
            testCase.rel.NonDetLabels  = NonDeterministicLabelling(testCase.sysAbs.outputs, ...
                testCase.sysLTI.regions, testCase.rel);
            

             

        end
    end
    %% Test Method Block
    methods (Test)
        %% Test basic functinality of EReachTime
        function testBasicLoad(testCase)
            createArguments(testCase)

            N = 30;     % time horizon
            
            [satProp,pol] = SynthesizeRobustController(testCase.sysAbs,testCase.DFA, ...
                                                                    testCase.rel, N, true);
            
            % check resulting policy
%             VisualizePolicy(testCase.sysAbs, pol, testCase.DFA,l, testCase.sysLTI)

        end

        function testTensorLab(testCase)
%             createArguments(testCase)
            
            %% with original computations
             N = 30;     % time horizon
            [satProp,pol] = SynthesizeRobustController(testCase.sysAbs,testCase.DFA, ...
                testCase.rel, N, true);         

             %% with tensor toolbox computations 
             sysAbsv2 =  testCase.sysAbs;
             sysAbsv2.P = TensorTransitionProbability_tensorlab(sysAbsv2.P.l,sysAbsv2.P.P_det,sysAbsv2.P.Pi{:});

             % change the probability transition matrix to tensorlab 

            [satProp_toolbox,pol_toolbox] = SynthesizeRobustController(sysAbsv2,testCase.DFA, ...
                testCase.rel, N, true);     
            verifyEqual(testCase,max(abs(satProp-satProp_toolbox)),0,"AbsTol", 0.01 )
            verifyEqual(testCase,max(max(abs(pol-pol_toolbox))),0,"AbsTol", 0.01 )

            end

             function testTensorToolbox(testCase)
%             createArguments(testCase)
            
            %% with original computations
             N = 30;     % time horizon
            [satProp,pol] = SynthesizeRobustController(testCase.sysAbs,testCase.DFA, ...
                testCase.rel, N, true);         

             %% with tensor toolbox computations 
             sysAbsv2 =  testCase.sysAbs;
             sysAbsv2.P = TensorTransitionProbability_tensortoolbox(sysAbsv2.P.l,sysAbsv2.P.P_det,sysAbsv2.P.Pi{:});

             % change the probability transition matrix to tensorlab 

            [satProp_toolbox,pol_toolbox] = SynthesizeRobustController(sysAbsv2,testCase.DFA, ...
                testCase.rel, N, true);     
            verifyEqual(testCase,max(abs(satProp-satProp_toolbox)),0,"AbsTol", 0.01 )
            verifyEqual(testCase,max(max(abs(pol-pol_toolbox))),0,"AbsTol", 0.01 )


        end
    end
end