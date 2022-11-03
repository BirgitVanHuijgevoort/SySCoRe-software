function EReachTime(sysAbs,DFA,pol,rel,  N)
    %EREACHTIME Compute the expected reach time for the satisfaction of the
    %temporal logic specification of interest. 
    % Written by: Sofie Haesaert
    % input: sysAbs, DFA, P, epsilon, delta, N,dim,l,uhat, initial only=false
    %
    uhat = sysAbs.inputs;

    disp('Warning this is untested code that is still under development')
    
    V = zeros(length(DFA.S),length(sysAbs.states)); 
    DFA_Active = setdiff(setdiff(DFA.S,DFA.F), DFA.sink);
    Converged(DFA_Active)=deal(0);
    outputs2act = rel.NonDetLabels;

    
    % prepare DFA transitions for all states states
    trans = ones(length(DFA.S),length(DFA.S),length(V));
    Pol = cell(1,length(DFA.S));
    for i = DFA_Active
        for l = 1: size(outputs2act,1)
            trans(i,DFA.trans(i,l),:) = min(shiftdim(trans(i,DFA.trans(i,l),:),1), 1-outputs2act(l,:));
        end
                    % compute probability transition matrix in the deterministic
            % case
            pol_d = reshape(pol(:,:,i), [[1],size(pol(:,:,i))]);
            [~,index_el] = min(abs(pol_d - uhat'), [],1);

            Pol{i} = sparse(index_el(:),1:length(index_el), 1, length(uhat),length(index_el) );
        
    end
    trans =10*trans;
    toc


for k = 1:N

    % Stop iterating when all values are converged
    if min(Converged) ==1
        disp('Convergence reached!')
        disp(['Number of iteration required, k=',num2str(k)])
        break;
    end
    
    for i = DFA_Active(end:-1:1) % for each discrete mode
        % check if mode has converged
        if min([Converged(DFA.trans(i,:)), Converged(i)])==1
            continue
        
        else 

            Converged(i)= 0;
        
            % Choose the correct states out of V based on the DFA
             V_sort = min(max(squeeze(trans(i, :,:)),V));
            

            % Compute value function for each action uhat
            Vi = V_sort*sysAbs.P*Pol{i}+1; % P is a object of the transition_probability class. 

           
       if min([Converged(DFA.trans(i,:))])==1
                           Converged(i)=1;

       elseif max(max(abs(V_n'-V(i,:))))<1e-6
                Converged(i)=1;
        
       end
        V(i,:) = Vi;
        
            
            
        end
    end
   
end

end