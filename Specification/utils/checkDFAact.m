function checkDFAact(DFA)
%checkDFAact Check if DFA.act is defined correctly.

for i = 1:length(DFA.act)
    assert(~isempty(DFA.act{1}), "Atomic proposition '%s' is empty! Try ' ' instead.", DFA.act{1})
end
end