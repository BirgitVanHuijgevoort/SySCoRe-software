function [satProp,pol,varargout] = SynthesizeRobustController(sysAbs,DFA,rel, N, varargin)
% Written by: Birgit van Huijgevoort, Sofie Haesaert
% input: sysAbs,DFA,rel, N, initialonly=false
%
% outputs: 
% satProp is robust satisfaction probability. s
% satProp(1,i) is satisfaction probability for DFA state 1 and initial 
% abstract state sysAbs.states(:,i).
%
% 
% Reference:
%  Haesaert, Sofie, and Sadegh Soudjani. 
%   "Robust dynamic programming for temporal logic control of stochastic systems." 
%   IEEE Transactions on Automatic Control (2020).

checkDFAact(DFA) % Check validity of DFA.act

disp(' <---- Start computing robust policy '); tic;

if nargin>=5
    if islogical(varargin{1})
        initialonly = varargin{1};
    else
       error(['5th input argument should be true/false for ' ...
           'the output of the value function and policy ' ...
           'for initial (discrete) state only'])  
    end
else
  initialonly = false;
end

if nargin>=6
   if islogical(varargin{2})
        antagonist_pol = varargin{2};
   else 
       error(['6th input argument should be true/false for ' ...
           'the computation of the counter strategy of the labeling'])

   end
else
  antagonist_pol = false;
end

uhat = sysAbs.inputs;
delta = rel.delta;
outputs2act = rel.NonDetLabels;
    
% Initialise value function
% V(i,j) is the probability of reaching F from DFA state i and abstract state
% sysAbs.states(:,j)
V = zeros(length(DFA.S),length(sysAbs.states));
DFA_Active = setdiff(setdiff(DFA.S,DFA.F), DFA.sink);
Converged = ones(1, length(DFA.S)); % create a vector with 
Converged(DFA_Active)=deal(0);

V(DFA.F,:) = ones(length(sysAbs.states),1); % Set V to 1 for DFA state q = F
tic    

% prepare DFA transitions for all states states
trans = ones(length(DFA.S),length(DFA.S),length(V));
for i = DFA_Active
    for l = 1: size(outputs2act,1)
        trans(i, DFA.trans(i, l), :) = min( ...
            shiftdim(trans(i, DFA.trans(i, l), :), 1), ...
            1 - outputs2act(l, :));
    end
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
            VV = V_sort*sysAbs.P; % P is a object of the transition_probability class

            % Optimize over uhat and subtract delta
            V_n = max(VV,[],2)-delta;  
            % Make sure that value function is between 0 and 1
            V_n= min(1,max(0,V_n));
           
       if min([Converged(DFA.trans(i,:))])==1
                           Converged(i)=1;

       elseif max(max(abs(V_n'-V(i,:))))<1e-6
                Converged(i)=1;
        
       end
        V(i,:) = V_n;
        
            
            
        end
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
a_pol = zeros(size(uhat,1),length(sysAbs.states),length(DFA.S));

for i = DFA_Active
        %q_old = DFA.S(i);
        
        % Choose the correct states out of V based on the DFA
        if antagonist_pol
                    [V_sort, a_pol] = min(max(squeeze(trans(i, :,:)),V));

        else
                    V_sort = min(max(squeeze(trans(i, :,:)),V));
        end

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
disp([' ----> Finished computing robust policy in ', num2str(toc)])

end

