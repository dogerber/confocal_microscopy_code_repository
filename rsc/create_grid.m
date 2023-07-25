function [X, Y, dx] = create_grid(d,max_nr_pts)
% Create the grid to interpolate the particle tracking data onto
% This section should remain unchanged if adapting to 3d TFM.
% [X, Y, dx] = create_grid(d,max_nr_pts)

fracpad = 0;

% Subtract off displacements from reference time
tref = 1;
for i=1:length(d)
    d(i).dr=d(i).dr-d(tref).dr;
end

% Select number of points for interpolated grid
ovr = 1; % Spatial oversampling (ovr=1 gives grid spacing= avg interparticle distance). ovr should be <=1.
nb_beads=length(d(1).r); % Total number of beads
nx=round(ovr*sqrt(nb_beads)); % Number of points on each side of the interpolation grid

if nargin>1
if nx > max_nr_pts
    %fprintf('[Experimental] Reducing grid points in create_grid!\n');
    % to large grids can not be used by futher tfm analysis
    nx = max_nr_pts;
end
end

if mod(nx,2)==0
    nx = nx+1; % Make sure odd number points in grid
end

% We pad the displacement data on each side of the grid. This reduces
% artifacts in the stress calculation. fracpad is the fraction of extra padding on
% each side of the original data. Thus fracpad=0.5 doubles the width and
% height of the original field of view

%fracpad=0.5;
npad = round(fracpad*nx);
% Calculate the boundaries of the data set
xmn = min(d(tref).r(:,1));
xmx = max(d(tref).r(:,1));
ymn = min(d(tref).r(:,2));
ymx = max(d(tref).r(:,2));

dx = max( (xmx-xmn)/nx, (ymx-ymn)/nx); % Distance between the grid points
c=.5*[xmn+xmx,ymn+ymx]; % Centre of data set


% Construct the grid
xi = linspace(-(nx-1)/2-npad,(nx-1)/2+npad,nx+2*npad)*dx+c(1);
yi = linspace(-(nx-1)/2-npad,(nx-1)/2+npad,nx+2*npad)*dx+c(2);
[X,Y]=meshgrid(xi,yi); % Matrix of gridpoints

end
