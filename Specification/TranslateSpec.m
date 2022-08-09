function [DFA] = TranslateSpec(formula,AP,debug)
% Written by Birgit van Huijgevoort, Sofie Haesaert
% input: LTL formula and atomic propositions
% output: struct DFA consisting of:
% states S
% initial states S0
% final states F
% actions act consists of strings that activate the transitions
% transitions trans with columns = actions as in DFA.act, row = state q, value = next state

%%
[B,alphabet] = spec2buchi(formula, AP);


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
                display('add state')
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
end

