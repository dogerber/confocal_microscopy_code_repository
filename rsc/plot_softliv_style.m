function plot_softliv_style()
% Displaying a graph in the unique SoftLiv style
wh = [1,1,1]; %RGB white
bk = [0,0,0];
dark_grey = 0.7 *wh;

% Setup
%f1 = figure;
%x = linspace(0,2*pi);
%p1 = plot(x,sin(x));
ax1 = gca; % get the axis
f1 = gcf;

% Axes color
ax1.XColor = dark_grey; % x-Axis color
ax1.YColor = dark_grey;
ax1.ZColor = dark_grey;

% Grid
ax1.GridColor = dark_grey;
ax1.GridAlpha = 0.5;

% background color
f1.Color = bk; % of the figure, so outside grid
ax1.Color = 0.5*dark_grey; % change background color of grid box

% Text size
ax1.FontSize = 16;

end