%% Translate LTL specification to Buchi automaton
% *formula* must be given between apostrophes (') ; propositions defined by A and b from before must be denoted by p1, p2, ..., px, where x=N_p; all N_p propositions should appear in formula (otherwise why did you define them?)
% -boolean operators in formula: & (for AND), | (for OR), ! (for NOT)
% -temporal operators in formula: G (for ALWAYS), F (for EVENTUALLY), U (for UNTIL), R (for RELEASES)
% 
% -it is better to use parantheses (to not mistake because the operators priority) and spaces (after ! you can omit the space)
% (example of formula: '(F p1) & (G !p2)')
% 
% *alphabeth* must be given as 1xn cell format, e.g., 
%    alphabet = {'p1'}    {'p2'}    {'p3'}    {'p4'}    {'p5'}
%
function [Buchi,Alph_s] = spec2buchi(formula,alphabet)
% 

Alph_s = alphabet_set(alphabet);  %set of possible combination of propositions that can appear (trans_system has observables in this set - indices in Alph_s)

Buchi = create_buchi(formula, Alph_s);

% create digraph
accepting =  cell(length(Buchi.S), 1);
init = cell(length(Buchi.S), 1);
for i =1:length(Buchi.S)
    accepting{i} = ismember(i,Buchi.F);
    init{i} = ismember(i,Buchi.S0);
end

NodeTable = table(Buchi.S', init, accepting,  ...
    'VariableNames',{'State' 'init', 'accepting'});  
sources= [];
targets = [];
sigma = {};
prop = {};
% Define edges

for s = Buchi.S
    for s_n =Buchi.S

        if ~isempty(Buchi.trans{s,s_n})
            sources = [sources;s];
            targets = [targets;s_n];
            sigma = {sigma{:}, Buchi.trans{s,s_n}};
            prop = {prop{:}, Buchi.prop{s,s_n}};
        end
    end


end
    EdgeTable = table([sources targets], sigma', prop',  'VariableNames',{'EndNodes' 'Sigma', 'Prop'});


Buchi.aut = digraph(EdgeTable,NodeTable);
end
