function [] = plotSatProb(satProb,sysAbs,varargin)
%PLOTSATPROB plots the satisfaction probability.
%
% plotSatProb(satProb,sysAbs) plots the satisfaction probability satProb for the
% states of the abstract system sysAbs [for all DFA states]
%
% plotSatProb(satProb,sysAbs, 'initial', DFA) plots the satisfaction probability satProb for the
% states of the abstract system sysAbs, but only for the initial DFA state

%% initialization
Initial = false;
for i = 1:length(varargin)
    % try to find 'initial'
    if strcmp(varargin{i},'initial')
        Initial = true;
        DFA = varargin{i+1};
        nS = length(DFA.S);
        break;
    end
end

StateShift = false;
for i = 1:length(varargin)
    % try to find 'shift'
    if strcmp(varargin{i},'shift')
        StateShift = true;
        xss = varargin{i+1};
        break;
    end
end

MOR = false;
for i = 1:length(varargin)
    % try to find 'MOR'
    if strcmp(varargin{i},'MOR')
        MOR = true;
        break;
    end
end

%% Plot results

if size(sysAbs.states,1) == 2   % 2D abstract model
    if Initial  % only plot satisfaction probability for initial DFA state
        % check if satProb is only computed for initial state
        if ~(size(satProb,1) == 1) % if this is not the case then convert
            % Determine correct q_0 for abstract states
            satProb = satProb((0:length(satProb)-1)*nS + DFA.trans(DFA.S0, sysAbs.labels));
        end
    end
    X1hat = sysAbs.hx{1};
    X2hat = sysAbs.hx{2};
        if StateShift
            X1hat = X1hat + xss(1);
            X2hat = X2hat + xss(2);
        end
    for i = 1:size(satProb,1)
        figure
        %surf(X1hat, X2hat, reshape(satProb(i,:), sysAbs.l(1), sysAbs.l(2)),'EdgeColor','interp')
        imagesc(X1hat, X2hat, transpose(reshape(satProb(i,:), sysAbs.l(1), sysAbs.l(2))))
        set(gca, 'Ydir', 'normal')
        if MOR
            xlabel('x_{r,1}')
            ylabel('x_{r,2}')
        else
            xlabel('x_1')
            ylabel('x_2')
        end
        zlim([0, 1])
        title('Robust satisfaction probability')    % for DFA state i
        colorbar
    end
elseif size(sysAbs.states,1) == 1   % 1D abstract model
    for i = 1:size(satProb,1)
        figure
        plot(sysAbs.states, satProb)
        xlabel('x')
        zlim([0, 1])
        title('Robust satisfaction probability')    % for DFA state i
        colorbar
    end
else
    warning('Plotting the satisfaction probability is only implemented for 1D or 2D abstract models.')
end

end