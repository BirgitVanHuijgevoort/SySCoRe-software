function [Polytope_contains] = EfficientContainment(Polytopearray, sysAbs )
containment = []        ;

% The sysAbs variable gives information about the gridding and the output
% map. 
map = sysAbs.outputmap ;
l = sysAbs.l ;
if sysAbs.zstates
    states = sysAbs.zstates; % check which gridding was used: transformed griddingt. 
else
    states = sysAbs.states; % check which gridding was used:  gridding of orig space. 
end
dim = size(states,1);


n_p = length(Polytopearray); % number of polytopes
        % find the domain of x

        X = Polyhedron('lb', min(states,[],2), 'ub', max(states,[],2));


        lb_inner = [];
        ub_inner = [];
        lb_outer = [];
        ub_outer = [];
        for p = 1:n_p
            Polytope_x = intersect(Polyhedron(Polytopearray(p).A*map,Polytopearray(p).b),X); 
            Polytope_x.computeVRep();
            % compute polytope for X and limit to domain of X
            [lb,ub] = innerApprox(Polytope_x) ;
            lb_inner = [lb_inner, lb];
            ub_inner = [ub_inner, ub];


              % TODO: Sofie continue here with the outer approx
              % computations. 
            lb_outer = [lb_outer, min(Polytope_x.V', [], 2) ];
            ub_outer = [ub_outer, max(Polytope_x.V', [], 2) ];

            % add to vector of polytopes
        end
        gridSize = (states(:, end)-states(:,1))./(l-1)';

         
        Polytope_contains = zeros(n_p,size(states,2));
        for p=1:n_p
        z_n = [lb_outer(:,p), lb_inner(:,p), ub_inner(:,p), ub_outer(:,p)];
        z_n_ind = ones(dim,1)+round(diag(gridSize.^-1)*(z_n-states(:,1)));
        
        % indices of states that are always true
        zi_indices= arrayfun(@(i)z_n_ind(i,2)+1:z_n_ind(i,3)-1,1:dim,'UniformOutput',false);
        ind_true = mat2cell(combvec(zi_indices{:}), ones(dim,1));

        Polytope_contains(p,sub2ind(l,ind_true{:})) = deal(1);


        % indices of states that are sometimes true
        zi_indices_outbound= arrayfun(@(i)[z_n_ind(i,1):z_n_ind(i,4) ],1:dim,'UniformOutput',false);
        zi_indices= arrayfun(@(i)[z_n_ind(i,1):z_n_ind(i,2),z_n_ind(i,3):z_n_ind(i,4) ],1:dim,'UniformOutput',false);
        ind_some_pre = [];
        for i = 1:dim
            zi = zi_indices_outbound;
            zi{i} = zi_indices{i};
            ind_some_pre = [ind_some_pre,combvec(zi{:})] ;
        end
        ind_some = mat2cell(ind_some_pre, ones(dim,1));
        indices_tot = unique(sort(sub2ind(l,ind_some{:})));
        critical_states = states(:,indices_tot);
        Polytope_contains(p,indices_tot) = Polytopearray(p).contains(critical_states);


        end


end