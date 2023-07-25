function [X, Y, Z] = create_grid_v2(pks,do_equal_spacing,dx)
% [X, Y, Z] = create_grid_v2(pks,do_equal_spacing,dx)
% if do_equal_spacing = false, it's equal number of grid points in all
% directions
%

do_debug = false;
if nargin <2
    do_equal_spacing = true; % grid-points are spaced apart the same distance in all directions
    % if false, there is the same number of grid points in every dimension
    % (n_pts)^(1/3)
end



%% import
n_dim = length(pks(1,:));

% Subtract off displacements from reference time
x = pks(:,1); y = pks(:,2);
xmn = min(x);
xmx = max(x);
ymn = min(y);
ymx = max(y);

if n_dim >2
    z = pks(:,3);
    zmn = min(z);
    zmx = max(z);
end

n_beads=length(x); % Total number of beads


%% Determine grid points
if do_equal_spacing %% determine grid points (equally spaced in all directions, average volume)
    if n_dim ==2
        if ~exist('dx','var'); dx = sqrt( (xmx-xmn) * (ymx-ymn) /n_beads); end
        
        [X,Y] = meshgrid(...
            linspace(xmn,xmx,round((xmx-xmn)/dx)),...
            linspace(ymn,ymx,round((ymx-ymn)/dx)));
        
    elseif n_dim ==3
        if ~exist('dx','var'); dx = ((xmx-xmn) * (ymx-ymn) * (zmx-zmn)/n_beads)^(1/3); end
        
        [X,Y,Z] = meshgrid(...
            linspace(xmn,xmx,round((xmx-xmn)/dx)),...
            linspace(ymn,ymx,round((ymx-ymn)/dx)),...
            linspace(zmn,zmx,round((zmx-zmn)/dx)));
    end
    
    
else %% determine grid points (equal number of pts in all directions)
    
    if n_dim ==2
        [X,Y] = meshgrid(...
            linspace(xmn,xmx,(n_beads)^(1/2)),...
            linspace(ymn,ymx,(n_beads)^(1/2)));
    elseif n_dim ==3
        [X,Y,Z] = meshgrid(...
            linspace(xmn,xmx,(n_beads)^(1/3)),...
            linspace(ymn,ymx,(n_beads)^(1/3)),...
            linspace(zmn,zmx,(n_beads)^(1/3)));
    end
    
    
end

if do_debug
    if n_dim == 2
        figure
        plot(x,y,'r.','DisplayName','input data');
        hold on
        plot(X(:),Y(:),'b.','DisplayName','grid points');
        axis equal; legend;
        shg
        fprintf('Initial points: %d \n Grid points: %d \n',length(x),length(X(:)));
    else
                figure
        plot3(x,y,z,'r.','DisplayName','input data');
        hold on
        plot3(X(:),Y(:),Z(:),'b.','DisplayName','grid points');
        axis equal; legend;
        shg
        fprintf('Initial points: %d \n Grid points: %d \n',length(x),length(X(:)));
    end
end


if size(pks,1)~=length(X(:))
    if do_debug; warning('create_grid_v2 did not return the same number of points as it received.'); end
end

end % function
