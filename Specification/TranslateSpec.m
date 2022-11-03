function [DFA] = TranslateSpec(formula,AP,debug)
%TRANSLATESPEC translates an scLTL specification to a deterministic
% finite-state automaton (DFA)
%
% Inputs
% ------
% formula = specification written in scLTL, for syntax see LTL2BA/README
% AP = atomic propositions
%
% Output struct DFA consisting of:
% -----
% states S
% initial states S0
% final states F
% sink states sink
% actions act consists of strings that activate the transitions
% transitions trans with columns = actions as in DFA.act, row = state q, value = next state
% 
% Basic reach-avoid example: 
% AP = {'p1','p2'};
% formula = '(!p2 U p1)';  % p1 = reach region, p2 = avoid region
% [DFA] = TranslateSpec(formula,sysLTI.AP)
%
% This code uses the tool LTL2BA, which can be found at http://www.lsv.fr/~gastin/ltl2ba/
% and is written by Denis Oddoux (v1.0) and modified by Paul Gastin (v1.2 & v1.3)
%
% Copyright 2022 Sofie Haesaert s.haesaert@tue.nl, Birgit van Huijgevoort b.c.v.huijgevoort@tue.nl

%% check if formula is in scLTL (instead of LTL)
disp('<---- Start translate specification')

% check if always (G) is not in formula
if ~isempty(strfind(formula, 'G'))
    error('Specification is not an scLTL specification, please input a different specification, or input a DFA.')
end

% check if negations are only given before APs
negs = strfind(formula, '!'); % find negations
% Get indices of atomic propositions
AP_index = [];
for j = 1:length(AP)
    index = strfind(formula,AP(j));
    AP_index = [AP_index, index];
end
flag = 0;
for i = 1:length(negs)
    % ! is directly followed by atomic proposition
    good = find(AP_index==negs(i)+1);
    % alternative: check if there is no operator after !
    if isempty(good)
        if formula(negs+1) ~= '&' & formula(negs+1) ~= '|' & formula(negs+1) ~= 'F' ...
                & formula(negs+1) ~= 'U' & formula(negs+1) ~= 'R' & formula(negs+1) ~= 'X'
            good = 1;
        end
    end
    if isempty(good)
        flag = 1;
    end
end
if flag
    error('Specification is not an scLTL specification, please input a different specification, or input a DFA.')
end


%% Construct Buchi automaton
[B,alphabet] = spec2buchi(formula, AP);

%% Clean up Buchi automaton
% 1. Clean up Buchi Automaton. Remove states from Buchi that do not have a
% possible loop. 
selfloop =[];
for s = B.F
    if isempty(B.trans{s,s})
        selfloop =[selfloop, 0];
    
    else
        selfloop = [selfloop, 1];
        
    end
    
end
B.F = B.F(selfloop==1);

if nargin>2 && debug
    figure;
    h=plot(B.aut,'EdgeLabel',B.aut.Edges.Prop);
end


% Check whether we accepting state help (1) as self loop. 
if length(B.F)~=1
    error('Wrong number of accepting states this software is not able to write the specification to a DFA' )
end

if length(B.S0) ~=1
    error('length(B.S0) ~=1: This feature has not been implemented')
end

%% Translate to DFA
% You now have a NFA. This still needs to be translated to a DFA

% step 1: create a table using cell
% from state | action to set of states | ...
Table_transitions = cell(length(B.S), length(B.S));
for s=B.S
    for act =1:length(alphabet)
        next_states = [];
        for s_n = B.S
           if ismember(act,B.trans{s,s_n})
               next_states = [next_states,s_n];
           end
        end
        Table_transitions{s,act} = next_states;
        
    end
end


% find finite states van DFA and transitions
dfa_states = {[0], B.S0};        % Add initial state of the Buchi to the DFA
s_index = 1;
trans = zeros(1, length(alphabet));
while length(dfa_states) > s_index
    s_index = s_index+1;
    s = dfa_states{s_index};
    trans(s_index, :)=NaN;
   
    
    for act =1:length(alphabet)
        s_found = unique([Table_transitions{s,act}]); % s can be a set of values, this action also sorts if needed       
        if ismember(B.F, s_found)  % transition to target is possible
            trans(s_index, act) = 1;
        else
        
        poss_trans = [cellfun(@(x) isempty(setxor(x, s_found)), dfa_states)]; % check existance of this set of states in the given states
        [m,i]=max(poss_trans);
            
            
            if max(poss_trans)==0 
                %display('add state')
                trans(s_index, act) = length(dfa_states)+1;

                dfa_states = {dfa_states{1:end}, s_found};
                
            else  
                trans(s_index, act) = i;
                
            end
            
        end
        
    end
    
    
end
DFA.S = 1:s_index;  % DFA states
DFA.S0 = 2;         % Initial states
DFA.F = 1;          % Accepting state is the first state. 
DFA.act = alphabet; % alphabet used as inputs for actions 
DFA.trans = trans;  % transitions
DFA.sink = DFA.S(all(DFA.trans == DFA.S',2)); % find modes with only self loops

disp('----> Finish translate specification')
end

