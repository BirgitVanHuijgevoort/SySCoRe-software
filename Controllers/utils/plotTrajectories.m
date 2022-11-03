function h = plotTrajectories(xsim, xlimits, sysLTI, varargin)
%plotTrajectories Plot ND trajectories supplied by ImplementController
% Can display 2D trajectories in a plane and separate subplots for general
% ND trajectories. Takes the state trajectories and plots the output
% trajectories. Can plot multiple trajectories and add coloring according
% to the DFA state, if this information is provided.
%
% Inputs:
% -------
% xsim = state trajectories (of type double for one trajectory and of type
% cell for several trajectories) as provided by ImplementController
% xlimits = bounds on the state space
% sysLTI = linear model
%
% Outputs:
% -------
% h = figure handle
%
% Options (= varargin)
% --------
% 'DFAinfo' - pass additional information about the evolution of the DFA
% state as provided by ImplementController as qsim. This will lead to the
% trajectories being plotted in colored sub-trajectories indicating the
% evolution of the DFA state. Expects input qsim
%
% 'shift' - add a shift of the state space. Expects input xss
% 
% Example: 
% plotTrajectories(xsim, xBounds, sysLTI, 'DFAinfo', qsim);
% see Tutorials/PackageDelivery for full example
% 
% Copyright Oliver Schoen, 2022

% Check whether DFA information supplied
DFAinfo = false;
StateShift = false;
for i = 1:length(varargin)
    % try to find 'initial'
    if strcmp(varargin{i}, 'DFAinfo')
        DFAinfo = true;
        qsim = varargin{i+1};
        break
    elseif strcmp(varargin{i}, 'shift')
        StateShift = true;
        xss = varargin{i+1};
        break
    end
end

% Plot system
figure
if sysLTI.dim == 2 && ~StateShift
    plot_x = plot(sysLTI.C*sysLTI.X);
    set(plot_x, 'FaceColor', 'None');
    hold on
    for i =1:length(sysLTI.regions)
        plot_x = plot(sysLTI.regions(i));
        set(plot_x, 'FaceColor', 'None');
        hold on
        xc = sysLTI.regions(i).chebyCenter;
        text(xc.x(1), xc.x(2), sysLTI.AP{i})
    end
end

% Add trajectories
if isa(xsim, 'cell')
    % Loop over multiple trajectories
    for i = 1:length(xsim)
        h = cell(length(xsim), 1);

        % Convert state trajectories to output trajectories
        if StateShift
            ysim = sysLTI.C*(xsim{i}+xss);
        else
            ysim = sysLTI.C*xsim{i};
        end

        % Plot i-th trajectory
        hold on
        if DFAinfo
            h{i} = plotTraj(DFAinfo, ysim, qsim{i});
        else
            h{i} = plotTraj(DFAinfo, ysim);
        end

        % Add initial point
        if sysLTI.dim == 2
            hold on
            scatter(ysim(1, 1), ysim(2, 1), 250, 'k+')
        end
    end
elseif isa(xsim, 'double')
    % Only one trajectory supplied
    % Convert state trajectories to output trajectories
    if StateShift
        ysim = sysLTI.C*(xsim+xss);
    else
        ysim = sysLTI.C*xsim;
    end

    % Plot trajectory
    if DFAinfo
        h = plotTraj(DFAinfo, ysim, qsim);
    else
        h = plotTraj(DFAinfo, ysim);
    end

    % Add initial point
    if sysLTI.dim == 2
        hold on
        scatter(ysim(1, 1), ysim(2, 1), 250, 'k+')
    end
else
    error("Type %s of xsim is not supported", class(xsim))
end

% Apply styles
title('')
if sysLTI.dim == 2
    grid on
    ylimits = sysLTI.C*xlimits;
    xlim(ylimits(1, :))
    ylim(ylimits(2, :))
    xlabel('$y_1$', 'Interpreter', 'latex')
    ylabel('$y_2$', 'Interpreter', 'latex')
end
end

function h = plotTraj(DFAinfo, ysim, qsim)
% Check dimension
% assert(size(ysim, 1) == 2, "Only two-dimensional trajectories supported. Dimension received is %d.", size(ysim, 1))

% Plots one trajectory
if size(ysim, 1) == 2
    % 2-dimensional case
    if DFAinfo
        % Plot individual pieces for DFA locations
        colorOrder = get(gca, 'ColorOrder');
        k = 1;
        iColor = 1;
        idcs = 1:length(qsim);
        while true
            % Get active DFA state
            q_active = qsim(k);

            % Get indices for which q stays the same
            log_active = (qsim(k:end) == q_active); % Ones where active
            k_cut = idcs(k:end);
            % Find the end by finding the first zero element
            log_active_stop = find(log_active - 1);
            if isempty(log_active_stop)
                % Only one DFA state in qsim
                log_active_stop = length(log_active);
            end
            log_active = log_active(1:log_active_stop);
            k_active = k_cut(log_active);
            k = k_active(end) + 1;
            % If not last state included make connection to next state
            if k_active(end) ~= length(qsim)
                k_active = [k_active, k_active(end)+1];
            end

            hold on
            if iColor > size(colorOrder, 1)
                % Out of colors
                colorOrder = [colorOrder; colorOrder];
            end
            h{iColor} = plot(ysim(1, k_active), ysim(2, k_active), ...
                'LineWidth', .5, 'Color', colorOrder(iColor, :), ...
                'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', colorOrder(iColor, :));
            hold off
            iColor = iColor + 1;

            % Check for end of trajectory
            if k > length(qsim) - 1
                break
            end
        end
    else
        hold on
        h = plot(ysim(1, :), ysim(2, :), ...
            'LineWidth', .5, 'Color', 'b', ...
            'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
        hold off
    end
else
    % n-dimensional case (n not equal 2)
    nDim = size(ysim, 1);
    for iPlot = 1:nDim
        subplot(nDim, 1, iPlot)
        ax = gca;
        hold on
        h = plot(ax, ysim(iPlot, :), ...
            'LineWidth', .5, 'Color', 'b', ...
            'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
        grid on
        xlabel('$t$', 'Interpreter', 'latex')
        ylabel(sprintf('$y_%d$', iPlot), 'Interpreter', 'latex')
        hold off
    end
end
end