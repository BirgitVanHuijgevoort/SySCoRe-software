function B = create_buchi(formula, Alph_s)
%Convert LTL formula to a Buchi automaton ltl2ba.exe (Oddoux&Gastin,
%http://www.liafa.jussieu.fr/ltl2ba/,
%http://www.liafa.jussieu.fr/~oddoux/ltl2ba/)(compiled using Cygwin) is
%used to convert an LTL formula to a Buchi automaton; automaton is in text
%form;

%this function will transforms it to a structure (tuple) (B=(S,S0,trans,F),
%where: S-states (identified by integers 1,2,...; S0-initial state(s);
%F-accepting (final) state(s); trans-structure containing transitions
%(trans{i,k} is a vector showing which element(s) of alphabet_set enable
%transition s_i->s_k) elements of alphabet_set are labeled with integers
%(1,2,...) (see function alphabet_set for construction of this set) (atomic
%props give observable space in which we are interested)

%inputs: alphabet_set=cell array of strings (ordered in a specific order -
%see function alphabet_set (2^N_p elements, N_p - no. of at props (without
%left-over space)) formula: string containing LTL formula (all its atomic
%proposition are in alphabet from which alphabet_set was created) (use
%spaces in formula) atomic propositions are of form p0...0x, where x is an
%integer; true, false Boolean operators: ! - negation; & - and; | - or; ->
%- implication; <-> - equivalence; Temporal operators: U - until; R -
%release; F - eventually; G - always; X - next example of formula: (F p1) &
%G !(p2 | p3)

%output: cell aray containing Buchi automaton accepting exactly infinite
%strings satisfying LTL formula

formula=regexprep(formula, '&', '&&');  %make changes in formula string so it matches the sintax required by ltl2ba.exe
formula=regexprep(formula, '\|', '||'); %| has special meaning, so it's preceded by \
formula=regexprep(formula, 'R', 'V');
formula=regexprep(formula, 'F', '<>');
formula=regexprep(formula, 'G', '[]');

currentFile = mfilename( 'fullpath' );
[pathstr,~,~] = fileparts( currentFile );

[s,r]=system([fullfile( './Specification/LTL2BA', '/ltl2ba/ltl2ba' ),' -f "' formula '"']); %sintax for calling ltl2ba.exe (located in subdir .\ltl2ba); use full description result (-d)
                                                        %example of run:
                                                        %[s,r]=dos('ltl2ba
                                                        %-d -f "<> p1"');
                                                        %%LTL 2 Buchi for
                                                        %our formula
%assert(s~=0, string(r)) % Check whether successful
% s = error code, s=0 is a successfull one
% r = output string that includes the Buchi Automaton  

% s = error code, s=0 is a successfull one
% r = output string that includes the Buchi Automaton
if s~=0 %error
    fprintf('\n Error when converting LTL to Buchi\n')
    return
end

%s=0, conversion successfull
sig=1:length(Alph_s);  %numeric labels for observables (integers)

str1=[char(10), 'Buchi automaton after simplification' char(10)];    %first string for match (char(10) - line feed)
%str2=[char(10), char(10), 'never {'];    %second string for match (end of
%interesting "zone" from what ltl2ba returns)
str2=[char(10), 'never {'];    %second string for match (end of interesting "zone" from what ltl2ba returns)

%remove from beginning until line contained by str1 (including it), and
%from beginning of str2 to end
r([1:findstr(r,str2)+length(str2)])=[];
%disp(r)

% we have a string like: "state x
%                         p1 -> state_label . . . state accept_** . . . "
% we want to create a cell B containing the Buchi automaton (see help on
% Characters and Strings, Regular Expressions)

%% check whether Buchi automaton is empty
% states_no=length(findstr(r,'state'));  %number of states in Buchi aut
% B.S=1:states_no;
if ~isempty(regexp(r,'empty automaton'))
    fprintf('\nError - empty Buchi automaton (maybe LTL formula has logical inconsistencies)\n');
end

%% find states, initial states and accepting states: B.S, B.S0, B.F 
[S_names, s_ind, e_ind]=regexp(r, '(?<state>\w*\_\w*)\:', 'match', 'start', 'end'); %S_names is a cell containing sub-cells with state names
                %s_ind, e_ind contain start&end indices for expressions
                %(useful for delimiting parts of r corresponding to a
                %certain state)
%% remove :
for i= 1:length(S_names)
    S_names{i} = S_names{i}(1:end-1);

end
% S_names=[S_names{:}] ;   %get rid of sub-cells (S_names is a cell array of strings)
states_no=length(S_names);  %number of states in Buchi aut
B.S=1:states_no;    %numeric indices for states
B.S0=find( cellfun( 'isempty', regexp(S_names,'init') ) == 0 );    %find indices of state(s) containing word 'init' (initial st.)
B.F=find( cellfun( 'isempty', regexp(S_names,'accept') ) == 0 );    %find indices of state(s) containing word 'accept' (accepting states)

%% find transitions (search succesors for each state and atomic proposition)
B.trans=cell(states_no,states_no);  %trans{i,k} gives the indices of elements from alphabet_set
%                                    (with "or" between them) for which s_i -> s_k
B.prop=cell(states_no,states_no);  %trans{i,k} gives the indices of elements from alphabet_set
%                                    (with "or" between them) for which s_i -> s_k

for i=1:states_no
    % Check whether state is an accept_all state
    if  startsWith('accept_all', S_names{i})  
        B.prop{i,i} = 'True';
        B.trans{i,i} = sig;
        continue
    end
    if i~=states_no % not the final state
        str=r((e_ind(i)+6): (s_ind(i+1)-5));   %select string containing transitions of current state (p1 -> . . .)
    else    %last state in string r
        str=r((e_ind(i)+3): end-2);
    end
    % str now includes the stransitions from state i to the other states
    row=regexp(str,'(\([^\n]*)\n','tokens');  %token: ([^\n] - any character different than new line), (* - any number of such characters)
%     celldisp(row)
    row=[row{:}]; %cell with rows for current state (i) (each row contains one state)
    for j=1:length(row)
        % k=find( cellfun( 'isempty', regexp(row{j}, S_names) ) == 0 );
        % %index of state in which s_i transit (on current row)
        %k should be an integer (could be vector, if there are states like
        %T0_S1, T0_S11, but this happens for unusually large formulas in
        %order to avoid this, use the following line (instead of the above
        %commented one):
        k=strmatch(row{j}((findstr(row{j},' -> ')+9) : end), S_names, 'exact'); %find string with state name and find index for exact matching in possible names
        % = find the goto state (its index in B.S) 

       %ONLY {& (and), ! (not), 1 (any=True)} can appear on each row with
       %respect to propositions (|| (OR) operator results in 2 rows) if 1
       %appears, it is the first element and there is no atomic proposition
       %on current row
        if row{j}(1)==1
            B.trans{i,k}=sig;   %for all possible combinations of propositions there is transition s_i -> s_k
            B.prop{i,k} = 'True';
            continue
        end

       %1 does not appear on current row
        prop=row{j}(1 : (findstr(row{j},' -> ')-1)); %delimitate proposition (expression) involving atomic propositions
                                                     %prop is only of kind
                                                     %"[!]pi & [!]pj &
                                                     %[!]pk" ([!] - !
                                                     %appears or not)
        B.prop{i,k} = prop;

        atom_pr=regexp(prop,'([!p]+\d+)','tokens'); %separate in atomic propositions (possibly preceded by !)
        atom_pr=[atom_pr{:}];
  
        labels=sig;  %will store labels of elements of alphabet_set that enable current transition (respecting expression)
                     %start with labels for whole alphabet_set, because
                     %we'll use intersections when updating vector "labels"
    
        for ap=1:length(atom_pr) %for each atomic prop, modify vector "labels"
            if isempty(findstr(atom_pr{ap},'!'))   %current atomic prop is not negated, so we keep ALL subsets that contain it,
                                                   %not only the subset
                                                   %equal with current
                                                   %atomic proposition
                            %use intersections because atomic propositions
                            %(possibly negated) are linked only by & (AND)
                            %operator
                labels=intersect(labels, find( cellfun( 'isempty', regexp(Alph_s,atom_pr{ap}) ) == 0 ));    %indices of subsets including current atomic prop.
            else    %negated, find all subsets that does NOT contain the current atomic proposition
                labels=intersect(labels, find( cellfun( 'isempty', regexp(Alph_s,atom_pr{ap}(2:end)) ) ~= 0 )); %add indices of all other subsets
            end
        end

        B.trans{i,k}=union(B.trans{i,k},labels);    %add current labels to current transitions (transition s_i -> s_k can be captured by more rows
                                                    %(equivalent with OR
                                                    %operator between
                                                    %propositions)
    end
end

if isempty(B.F) %no final state detected (because it is only one initial state and ltl2ba names it only "init" in description of simplified autmaton)
    B.F=B.S0;
end


end
