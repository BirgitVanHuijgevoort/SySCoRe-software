%% Example of Reach avoid specification
%  Example of working with some LTL spec
%  atomic proposition p1, p2, ... (other formats are nto well translated)
% -boolean operators in formula: & (for AND), | (for OR), ! (for NOT)
% -temporal operators in formula: G (for ALWAYS), F (for EVENTUALLY), U (for UNTIL), R (for RELEASES)

%% 

formula = '( G F p1 ) & G p2';
AP = {'p1', 'p2'};
[B,alphabet] = spec2buchi(formula, AP);

% alph = the powerset with the letters of the alphabet 
% B = the Buchi automation, with
%       - B.S the set of states
%       - B.S0 the initial state
%       - B.F the set of accepting states
%       - B.prop the propositions on the transitions of the Buchi
%       - B.trans the transitions labelled with the indicators of the
%       allowed alphabet letters. (indices of alp)
%       - B.aut digraph representation of  Buchi automaton (Good for
%       plotting)