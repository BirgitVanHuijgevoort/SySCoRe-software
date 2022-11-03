function Plot_sysLTI(LinModel)
%PLOT_SYSLTI Plot the important regions of the LTI system.
%
% Plot_sysLTI(LinModel) plots the regions of a 1D or 2D output space 
% based on LinModel.regions


if size(LinModel.C,1)>2
    error('This function does not work for systems with a non-2D space')
end
%1. plot the domain of X
if LinModel.dim==2
    figure;
    plot_x = plot(LinModel.X);
    set(plot_x, 'FaceColor','None');
    title('X domain')
elseif LinModel.dim==1
    figure;
    xpoints = linspace(LinModel.X.V(1), LinModel.X.V(2), 100);
    ypoints = zeros(size(xpoints));
    plot(xpoints,ypoints, 'k', 'LineWidth',2)
    hold on
    xlim([min(LinModel.X.V), max(LinModel.X.V)])
    ylim([-0.1 0.1])
end

%2. plot the regions of the LTI system together with the AP strings
if LinModel.dim == 2
    figure;
    plot_x = plot(LinModel.C*LinModel.X);
    set(plot_x, 'FaceColor','None');
    hold on 
    for i =1:length(LinModel.regions)
        plot_x = plot(LinModel.regions(i));
        set(plot_x, 'FaceColor','None');
        hold on
        xc = LinModel.regions(i).chebyCenter;
        text(xc.x(1),xc.x(2), LinModel.AP{i})
    end
elseif LinModel.dim == 1
    for i = 1:length(LinModel.regions)
        xpoints = linspace(min(LinModel.regions(i).V), max(LinModel.regions(i).V), 100);
        ypoints = zeros(size(xpoints));
        plot(xpoints,ypoints, '+');
        hold on 
        xc = min(LinModel.regions(i).V)+(max(LinModel.regions(i).V)-min(LinModel.regions(i).V))/2;
        text(xc,0.01, LinModel.AP{i})
    end
    xlim([min(LinModel.X.V), max(LinModel.X.V)])
    ylim([-0.1 0.1])
end

title('Output space Y with labelled regions' )

end

