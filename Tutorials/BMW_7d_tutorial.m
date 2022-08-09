%% BMW_7d tutorial

% Load model into sysNonLin
BMW_7d_model

% Bounds on state space 
x1l = -10;   % Lowerbound x1
x1u = 10;   % Upperbound x1
x2l = -10;   % Lowerbound x2
x2u = 10;   % Upperbound x2
sysNonLin.X = Polyhedron(combvec([x1l,x1u],[x2l,x2u])');

% Bounds on  and input space
ul = [-1;-1];   % Lowerbound input u
uu = [1;1];     % Upperbound input u
sysNonLin.U = Polyhedron(combvec([ul(1),uu(1)],[ul(2),uu(2)])');


% Specify regions for the specification
%his should be specified over the output of the  model not over
% the states. (Even though in this case this is the same)
p1x = [4 4 10 10 4];    % x1-coordinates
p1y = [0 -4 -4 0 0];    % x2-coordinates
p1 = [p1x; p1y];         % parking region
P1 = Polyhedron(p1');

p2x = [4 4 10 10 4];    % x1-coordinates
p2y = [0 4 4 0 0];      % x2-coordinates
p2 = [p2x; p2y];        % avoid region
P2 = Polyhedron(p2');

sysLTI.regions = [P1;P2]; % regions that get specific atomic propositions
sysLTI.AP = {'p1', 'p2'}; % with the corresponding atomic propositions

%% Synthesize scLTL formula (or input DFA yourself)
%%% use LTL2BA and check if determinstic and accepting state with loop with 1.
% input: (sc)LTL formula and atomic propositions (see readme in folder
% LTL2BA)
% output: struct DFA containing (among other) the transitions

formula = '(!p2 U p1)';  % p1 = parking, p2 = avoid region
% formula should use atomic propositions in sysLTI.AP. 
