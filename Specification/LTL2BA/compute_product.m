function P = compute_product(TS,B, x2Sigma)
% Compute the product of a transition system and a Buchi automaton. See
% slides lecture 2 module 4 part 1.
varnames = cell(TS.nnodes*B.nnodes,1); % empty variable name list

P.nnodes = TS.nnodes*B.nnodes;
varnames_in = zeros(TS.nnodes,B.nnodes);

P.F = []; % accepting set

%% Build node list for P= B xS 
index= 1;
for x = 1:TS.nnodes
    for s = 1:B.nnodes
        varnames{index} = [ num2str(x) num2str(s)];
        varnames_in(x,s) = index;
        index = index + 1;

    end
end
%% collect edge info for P= B xS 
prop = {};l_n_= [];utr = [];
P.sources = []; % empty source node list
P.targets = []; % empty target node list
index= 1;
for x = 1:TS.nnodes
    for s = 1:B.nnodes
        
        
        for x_n = successors(TS.S, x)'
            % for all next transitions of x,
            % find next transition of Buchi
            l_n = x2Sigma(x_n) ; % next letter
            for s_n = 1:B.nnodes
                if ismember(l_n, B.trans{s,s_n})
                    P.sources = [P.sources; index];
                    P.targets = [P.targets; varnames_in(x_n,s_n)];
                    prop = {prop{:},B.prop{s,s_n} };
                    utr =  [utr; TS.S.Edges.U(all(TS.S.Edges.EndNodes == [x x_n], 2))];
                    l_n_ = [l_n_;l_n];
                end
            end
        end
        
        if ismember(s, B.F)
            P.F = [P.F; index];
        end

        index = index + 1;
    end
end

P.EdgeTable = table([P.sources P.targets], prop', l_n_, utr, 'VariableNames',{'EndNodes' 'Prop' 'Sigma' 'U'});



%% collect node info
% Add initial state to set of initial states
P.S0 = [];
for x =TS.X0
    l_n = x2Sigma(x);
    for s = B.S0
        for s_n = 1:B.nnodes
            if ismember(l_n, B.trans{s,s_n})
                P.S0 = [P.S0; varnames_in(x,s_n)];
            end
        end

    end
end

% build info that can be used via the graph later on
Sstate =cell(P.nnodes, 1);Bstate = cell(P.nnodes, 1);Pstate = cell(P.nnodes, 1);Stateinf =cell(P.nnodes, 1);
index=1;
for x = 1:TS.nnodes
    for s = 1:B.nnodes
        Sstate{index} = num2str(x) ;
        Bstate{index} = num2str(s) ;
        Pstate{index} = num2str(index) ;
        Stateinf{index} = ' ';
        if ismember(index, P.S0)
            Stateinf{index} =  [Stateinf{index} ' init '];
        end
        if ismember(index, P.F)
            Stateinf{index} =  [Stateinf{index} 'accepting'];
        end
        index = index + 1;
    end
end
P.NodeTable = table(Pstate,Stateinf, Sstate, Bstate,...
    'VariableNames',{'PState' 'info' 'SState' 'BState'});    
P.T = digraph(P.EdgeTable,P.NodeTable);