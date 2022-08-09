function [linMod, R] = ReduceX(linMod, inputset,  output_region, type, ...
    horizon, varargin)
% REDUCEX reduce the polytopic statespace set based on an output set and a
% property of either time-bounded invariance or reachability (latter not
% implemented yet)

% Check whether conditions under which this code holds are satisfied
if ~strcmp(type, 'invariance')
    error(['The type = ', type, ' is not yet known or implemented'])

end

%% Preliminaries
R_init = Polyhedron(output_region.A*linMod.C, output_region.b);

system = LTISystem('A',  linMod.A, 'B', linMod.B);
system.x.min = min(linMod.X.V, [], 1);
system.x.max = max(linMod.X.V, [], 1);
system.u.min = inputset(:,1);
system.u.max = inputset(:,2);

if linMod.U.Dim ~=1
    warning('works for 1 inputs only')
end
U = Polyhedron(inputset');

% %% Debugging figures
% figure(1)
% subplot(2,1,1)
% hold off
% R_init.plot()
% title(['backward invariance set i = ', num2str(1)])

%% Init reach set
R= cell(horizon,1);

R{1}= R_init;
V = [ ];
for i=1:horizon
    R_ = system.reachableSet('X', R{i}, 'U', U, 'N', 1, ...
        'direction', 'backward');
    R{i+1} = R_ & R_init;
    V = [V;R{i+1}.V];
end

%% Debugging figures
%debugfigs(i, R)

%% Retrieve reduced model
X = Polyhedron(V);
X.minVRep
linMod.X = X;
end


function debugfigs(i, R)
% Debugging figures
figure(1)
tl = tiledlayout(2, i+1);

% Print backward invariant sets
for j = 1:i+1
    ax = nexttile;
    R{j}.plot();
    title(ax, ['Backward invariance, i = ', num2str(j)])
    subtitle(ax, 'Invariant set')
    xl{j} = xlim(ax);
    yl{j} = ylim(ax);
end

% Print vertices of backward invariant sets
for j = 1:i+1
    ax = nexttile;
    scatter(R{j}.V(:,1), R{j}.V(:,2), 60, [204 204 255]/256, 'filled')
    subtitle(ax, 'Vertices')
    xlim(ax, xl{j})
    ylim(ax, yl{j})
    grid on
end
end