function pks_out = outliers3_fast(pks)
%% OUTLIERS3_FAST Detects outliers in a 3d-point cloud that belongs to a surface
% The algorithm compares roughly each point to it's vicinity and the local standarddeviation
% in z-coodinate. It removes points that are to far out.
%
%  Input:
%       pks: matrix with all the points [x,y,z], Nx3
%
%  Output:
%       pks_out: pks without outliers
%
% Options are in the begining of the code.
% 
% ---------------------------------------------- Dominic Gerber, June 2020




%% Parameters
param.do_debug = false; % gives you visual feed back on what's happening

max_derivation = 3; % how many standarddeviation are outliers away from the mean?

%nr_xsec = 10; % into how many tiles should the x-direction be seperated?
%nr_ysec = 10;


%% Preparation
pks = [pks zeros(1,length(pks(:,1)))']; % allocate 4 column
sec_pks_sum = nan(size(pks));

% determine how many tiles to form
nr_pts_total = length(pks(:,1));
estimate_tiles_nr = nr_pts_total/200; % group 200 pts together, roughly
% assuming more or less uniform distribution here
nr_xsec = round(sqrt(estimate_tiles_nr));
nr_ysec = nr_xsec;

% determine xy extents
x0 = min(pks(:,1));
y0 = min(pks(:,2));
xmx = max(pks(:,1));
ymx = max(pks(:,2));

% tile widths
x_sec_width = abs(xmx-x0)/nr_xsec;
y_sec_width = abs(ymx-y0)/nr_ysec;


if param.do_debug
    figure;
    scatter3m(pks(:,1:3));
    hold on
end

%% Actual loop
for x_sec = 0:nr_xsec-1
    xlowend = x0+x_sec*x_sec_width;
    xhighend = x0+(x_sec+1)*x_sec_width;
    
    idx = and(pks(:,1)>(xlowend),...
        (xhighend)> pks(:,1));
    
    sec_pks_x = pks(idx,:);
    
    
    for y_sec = 0:nr_ysec-1
        idx_y = and(sec_pks_x(:,2)> (y0+y_sec*y_sec_width),...
            (y0+(y_sec+1)*y_sec_width)>sec_pks_x(:,2));
        sec_pks = sec_pks_x(idx_y,:);
        
        if length(sec_pks(:,1))<50
            sec_pks(:,4) = mean(sec_pks(:,3)); % set mean so they are not removed
            sec_pks_sum = [sec_pks_sum;sec_pks];
            continue; % not enough pts for good statistics
        end
        
        %actual calculation
        sec_mn = mean(sec_pks(:,3));
        sec_std = std(sec_pks(:,3));
        
        sec_pks(:,4) = abs((sec_pks(:,3)-sec_mn))/sec_std; % distance from mean in realtion to std
        sec_pks_sum = [sec_pks_sum;sec_pks];
        
        if param.do_debug
            scatter3(sec_pks(:,1),sec_pks(:,2),sec_pks(:,3),[],sec_pks(:,4),'x');
            beep;
        end
    end
    
end

%% Export
idx_out = (sec_pks_sum(:,4) < max_derivation); % these are to be cut out
pks_out = sec_pks_sum(idx_out,:);
pks_out = pks_out(~isnan(pks_out(:,1)),:);


if param.do_debug
    figure('Name','Cleaned');
    scatter3(pks_out(:,1),pks_out(:,2),pks_out(:,3));
end


end % main function

