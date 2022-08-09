function [satProp,pol] = SynthesizeController(sysAbs,DFA,N,delta, initialonly)
% Written by: Birgit van Huijgevoort
% input: sysAbs, DFA, P, epsilon, delta, N,dim,l,uhat, initial only=false
%
% outputs: 
% output 1: robust satisfaction probability satProp. SatProp(1,i) is satisfaction 
% probability when initial abstract state equals sysAbs.states(:,i).

% output 2: robust control policy pol, with pol(:,i) indicating the optimal uhat when
% the current state equals sysAbs.states(:,i).


if nargin<4
  delta = 0;
end
if nargin<5
  initialonly = true;
end

uhat = sysAbs.inputs;
% Initialise value function
% V(i,j) is the probability of reaching F from DFA state i and abstract state
% sysAbs.states(:,j)
V = zeros(length(DFA.S),length(sysAbs.states));
V_new = zeros(length(DFA.S),length(sysAbs.states));
DFA_Active = setdiff(setdiff(DFA.S,DFA.F), DFA.sink);
V_new(DFA.F,:) = ones(length(sysAbs.states),1);
    tic    
for k = 1:N

    % Set V to 1 for DFA state q = F

    % Stop iterating when all values are converged
    if max(max(abs(V_new-V)))<1e-6
        disp('Convergence reached!')
        disp(['Time steps required, k=',num2str(k)])
        break;
    end
    V = V_new;
    
    for i = DFA_Active % for each discrete mode
        %q_old = DFA.S(i);

        % Choose the correct states out of V based on the DFA
        V_sort = V((0:length(V)-1)*length(DFA.S)+DFA.trans(i,sysAbs.labels));

        % Initialise value iteration with including action uhat.  
        % VV = zeros(l^dim+1,N,size(uhat,2));
%         VV = zeros(size(V,2),size(uhat,2));   % without sink state
        
        % Compute value function for each action uhat
        VV = V_sort*sysAbs.P; % P is a object of the transition_probability class. 
        V_new(i,:) = max(max(VV,[],2)-delta,0);   % optimize over uhat and subtract delta

    end
   
end
 toc
%% Compute satisfaction probability

if initialonly==true
    % Determine correct q_0 for abstract states
    satProp = V((0:length(V)-1)*length(DFA.S)+DFA.trans(DFA.S0,sysAbs.labels));
else
    satProp = V;
end

%% Compute optimal policy

pol = zeros(size(uhat,1),length(sysAbs.states),length(DFA.S));
for i = DFA_Active
        %q_old = DFA.S(i);
        
        % Choose the correct states out of V based on the DFA
        V_sort = V((0:length(V)-1)*length(DFA.S)+DFA.trans(i,sysAbs.labels));


        % Initialise value iteration with including action uhat.  
        VV = V_sort*sysAbs.P; % P is a object of the transition_probability class. 

        [~,index_pol] = max(VV,[],2);   % optimize over uhat 
        pol(:,:,i) = uhat(:,index_pol);

end

if initialonly==true
    
% Determine correct q_0 for abstract states 
policy = zeros(size(uhat,1),length(sysAbs.states));
for i = 1:length(sysAbs.states)
    policy(:,i) = pol(:,i,DFA.trans(DFA.S0,sysAbs.labels(i)));
end 
pol = policy;

    
end

end

