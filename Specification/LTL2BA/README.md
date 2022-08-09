# LTLcontrol

Beta version with basic scripts needed for LTL control and model checking via Matlab.

The routines create_buchi, alphabet_set, spec2buchi are based on code from  http://sites.bu.edu/hyness/ltl-control/
The adapted version now works well on a mac and on windows, it has also been adapted to start from the promela code that the LTL2BA generates. For Linux you can compile a new ltl2ba executable via the ltl2ba website (see http://www.lsv.fr/~gastin/ltl2ba/download.php)


Usage LTL FORMULA:

formula must be given between apostrophes (') ;
-boolean operators in formula: & (for AND), | (for OR), ! (for NOT)
-temporal operators in formula: G (for ALWAYS), F (for EVENTUALLY), U (for UNTIL), R (for RELEASES), X (for next)

-it is better to use parentheses (because of the operators priority) and spaces (after ! you can omit the space)
(example of formula: '(F p1) & (G !p2)')

EXAMPLES:
- Try out ExampleLTL

Authors: Sofie Haesaert, M.Kaleem Salahuddin
