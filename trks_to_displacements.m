function [d,trks_out] = trks_to_displacements(trks,param)
% [d,trks_out] = trks_to_displacements(trks,param)
%
% Alternative to track2disp(), which does not throw away incomplete traces,
% but instead interpolates a flow field at every timepoint.
%
%
% INPUT:
%    - trks as  [x,y, t, id] or [x,y,z, t, id]
%
% OUTPUT:
%    - d
%    - trks
%
% Reasoning: start with grided points and let them go through all
% timepoints. For each timepoint interpolate the expected displacement and
% add them onto the grided points
%
%
% Originally for Bsc thesis Matthias Kuster, improved for Shaohua
%
%
% ------------------------------------------- Dominic Gerber, December 2022

if nargin <2
    param = struct();
end

if ~isfield(param,'do_debug'); param.do_debug = true; end % shows a lot of infos about the internal working
if ~isfield(param,'do_save'); param.do_save = false; end % %
if ~isfield(param,'do_talk'); param.do_talk = true; end %
if ~isfield(param,'do_show'); param.do_show = true; end % %
if ~isfield(param,'do_remove_outliers'); param.do_remove_outliers = false; end %
if ~isfield(param,'do_remove_moved_outside_fov'); param.do_remove_moved_outside_fov = true; end %
if ~isfield(param,'do_equal_spacing'); param.do_equal_spacing = true; end % only when using grid
if ~isfield(param,'do_make_grid'); param.do_make_grid = false; end % if false, will use the initial particle positions



%% prepare trks and starting points

if nargin<1
    [ffname,ffolder] = uigetfile(cd,'Choose the IMAGE you want to analyze','*.mat');
    FULLFILE_IN_MAT = fullfile(ffolder, ffname);
    load(FULLFILE_IN_MAT,'trks');
    save_name = extractBefore(ffname,'_disp');
else
    if param.do_save
        save_name = input('Sample name: ','s');
    end
end

if param.do_remove_outliers
    trks = remove_outliers_RAFT_DG(trks,'dim',3,'abs_diff',100,'angle_diff',40,'show_result',param.do_debug,'min_number_neighbours',20); % remove badly tracked displacements
    if param.do_talk; fprintf('Outliers removed.\n'); end
end

% measure input sizes
trks = double(trks);
n_dim = size(trks,2)-2;
unique_timepoints = unique(trks(:,end-1));
n_timepoints = length(unique_timepoints); % old max(trks(:,end-1));
n_pts = floor(length(trks)/n_timepoints); % average number of points per timepoint

colors_t = cbrewer('div','Spectral',n_timepoints);

xlims=[min(trks(:,1)),max(trks(:,1))];
ylims=[min(trks(:,2)),max(trks(:,2))];
avg_spacing = sqrt((xlims(2)-xlims(1))*(ylims(2)-ylims(1))/n_pts);
if n_dim >2; zlims = [min(trks(:,3)),max(trks(:,3))]; end



% -- make grid
if param.do_make_grid
    if n_dim <3
        if ~param.do_equal_spacing % equal number of points in both dimensions
            if param.do_talk; fprintf('Gridding-Mode: equal number in both dimensions.\n'); end
            n_grid_pts = floor(sqrt(n_pts))^2; %
            [X_grid,Y_grid]=meshgrid(linspace(xlims(:,1),xlims(:,2),sqrt(n_grid_pts)),...
                linspace(ylims(:,1),ylims(:,2),sqrt(n_grid_pts)));
        else % euqal grid step in both dimensions
            if param.do_talk; fprintf('Gridding-Mode: equal spacing in both dimensions.\n'); end
            [X_grid,Y_grid]=meshgrid([xlims(:,1):avg_spacing:xlims(:,2)],...
                [ylims(:,1):avg_spacing:ylims(:,2)]);
            n_grid_pts = length(X_grid(:));
        end
    else % 3d
        [X_grid, Y_grid, Z_grid] = create_grid_v2(trks( trks(:,end-1)==1,1:3 ),param.do_equal_spacing);
        n_grid_pts = length(X_grid(:));
    end
    
    
else % use actual particle positions instead
    
    if false % only use from first timepoint
        if param.do_talk; fprintf('Using initial particle positions instead of a grid.\n'); end
        idx_trks_t0 = (trks(:,end-1)==min(trks(:,end-1)));
        X_grid = trks(idx_trks_t0,1);
        Y_grid = trks(idx_trks_t0,2);
        if n_dim >2; Z_grid = trks(idx_trks_t0,3); end
        
        n_grid_pts = sum(idx_trks_t0);
    else % use from all timepoints
        if param.do_talk; fprintf('Using particle positions instead of a grid, with points from all timepoints.\n'); end
        sensible_dist = 10;
        tps = unique(trks(:,end-1));
        pks = trks( trks(:,end-1)== tps(1), 1:3);
        for i_tp = 2:length(tps)
            % only add the ones that are in an empty space
            pks_itp = trks( trks(:,end-1)== tps(i_tp), 1:3);
            
            idx = rangesearch(pks,pks_itp,sensible_dist);
            idx_no_nearest_neighbours = cellfun(@isempty,idx);
            
            pks = [pks; pks_itp(idx_no_nearest_neighbours,:)];
            
            %figure; scatter3m(pks); title(i_tp);
        end
        
% convert to grid
        X_grid = pks(:,1);
        Y_grid = pks(:,2);
        if n_dim >2; Z_grid = pks(:,3); end
        
        n_grid_pts = length(X_grid);
    end
    
    
end


x_traces = zeros(n_grid_pts,n_timepoints); % all positions (NxM) of N particles over M timepoints
y_traces = zeros(n_grid_pts,n_timepoints);
if n_dim >2; z_traces = zeros(n_grid_pts,n_timepoints); end

% write starting points down
x_traces(:,1) = X_grid(:);
y_traces(:,1) = Y_grid(:);
if n_dim >2; z_traces =  Z_grid(:); end



%% calculate displacements from track

for t=1:n_timepoints-1
    
    % find parts of trks that are in this and the next timepoint
    sub_trks_t = double(trks(trks(:,end-1)==unique_timepoints(t),:));
    sub_trks_t1 = double(trks(trks(:,end-1)==unique_timepoints(t+1),:));
    
    [~,it,it1] = intersect(sub_trks_t(:,end),sub_trks_t1(:,end)); % find trks-id that are in both timepoints
    
    if isempty(it)
        fprintf('Current timepoint: %d\n',t);
        warning('Can not find any connection between these two timepoints!');
        warning('keeping it here the same as last tp to make the code run');
        x_traces(:,t+1) = x_traces(:,t);
        y_traces(:,t+1) = y_traces(:,t);
        if n_dim >2; z_traces(:,t+1) = z_traces(:,t); end
    else
        % reform to only contain trks-pts that are in both timepoints
        pos = sub_trks_t(it,:);
        difference = sub_trks_t1(it1,1:n_dim)-sub_trks_t(it,1:n_dim);
        
        % remove unbelievably large differences
        if false
            difference(isoutlier(difference,'grubbs')) = 0;
        end
        
        % interpolate the movement step ux and uy
        if n_dim <3 % 2d
            F_ux = scatteredInterpolant(pos(:,1),pos(:,2),difference(:,1));
            F_uy = scatteredInterpolant(pos(:,1),pos(:,2),difference(:,2));
            
            % add the movement step to current pos of the traces
            x_traces(:,t+1) = F_ux(x_traces(:,t),y_traces(:,t)) + x_traces(:,t);
            y_traces(:,t+1) = F_uy(x_traces(:,t),y_traces(:,t)) + y_traces(:,t);
        else % 3d
            if true % using scatteredInterpolant
                F_ux = scatteredInterpolant(pos(:,1),pos(:,2),pos(:,3),difference(:,1));
                F_uy = scatteredInterpolant(pos(:,1),pos(:,2),pos(:,3),difference(:,2));
                F_uz = scatteredInterpolant(pos(:,1),pos(:,2),pos(:,3),difference(:,3));
                
                % add the movement step to current pos of the traces
                x_traces(:,t+1) = F_ux(x_traces(:,t),y_traces(:,t),z_traces(:,t)) + x_traces(:,t);
                y_traces(:,t+1) = F_uy(x_traces(:,t),y_traces(:,t),z_traces(:,t)) + y_traces(:,t);
                z_traces(:,t+1) = F_uz(x_traces(:,t),y_traces(:,t),z_traces(:,t)) + z_traces(:,t);
            else % using robustfit
                fprintf('Using robustfit() method, experimental. \n');
                %%%%%%%%%%%%%%%%%%% SKETCH! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % warning('THIS IS A SKETCH');
                % decide on neighbour groups
                idx_near = knnsearch(pos(:,1:3), [x_traces(:,t),y_traces(:,t),z_traces(:,t)],'k',20);
                % maybe rangeserach instead?
                
                % loop over all grid points and get prediction from
                % neighbours
                for idx_pt = 1:length(x_traces(:,t))
                    
                    for idx_dim = 1:n_dim % for each dimension find the coefficients to predict it from the fulldimension input
                        mdllr = fitlm( pos(idx_near(idx_pt,:),1:3),...
                            difference(idx_near(idx_pt,:),idx_dim),...
                            'RobustOpts','on');
                        
                        if idx_dim==1
                            x_traces(idx_pt,t+1) = x_traces(idx_pt,t) + predict(mdllr,...
                                [x_traces(idx_pt,t),y_traces(idx_pt,t),z_traces(idx_pt,t)]);
                        elseif idx_dim==2
                            y_traces(idx_pt,t+1) = y_traces(idx_pt,t) + predict(mdllr,...
                                [x_traces(idx_pt,t),y_traces(idx_pt,t),z_traces(idx_pt,t)]);
                        elseif idx_dim==3
                            z_traces(idx_pt,t+1) = z_traces(idx_pt,t) + predict(mdllr,...
                                [x_traces(idx_pt,t),y_traces(idx_pt,t),z_traces(idx_pt,t)]);
                        end
                        
                        
                        % A(idx_dim, :) = mdllr.Coefficients.Estimate([2:4,1])';
                    end
                    
                end % idx_pt
            end
        end
        
        
        % OPTIONAL: traces moved out of FOV are killed (?)
        if param.do_remove_moved_outside_fov
            if n_dim <3
                idx_out_of_fov = any([x_traces(:,t+1)< xlims(1),...
                    x_traces(:,t+1)> xlims(2),...
                    y_traces(:,t+1)< ylims(1),...
                    y_traces(:,t+1)> ylims(2)],2);
                
                x_traces(idx_out_of_fov,:) = nan;
                y_traces(idx_out_of_fov,:) = nan;
            else
                idx_out_of_fov = any([x_traces(:,t+1)< xlims(1),...
                    x_traces(:,t+1)> xlims(2),...
                    y_traces(:,t+1)< ylims(1),...
                    y_traces(:,t+1)> ylims(2),...
                    z_traces(:,t+1)< zlims(1),...
                    z_traces(:,t+1)> zlims(2)],2);
                
                x_traces(idx_out_of_fov,:) = nan;
                y_traces(idx_out_of_fov,:) = nan;
                z_traces(idx_out_of_fov,:) = nan;
            end
            if t==1; warning('Traces outside of intial FOV are removed'); end
        end
        
        
        % debug: check
        if param.do_debug
            if n_dim ==2
                figure
                quiver(pos(:,1),pos(:,2),difference(:,1),difference(:,2),0,'blue')
                hold on
                quiver(x_traces(:,t),y_traces(:,t),...
                    x_traces(:,t+1)-x_traces(:,t),...
                    y_traces(:,t+1)-y_traces(:,t),...
                    0,'red');
                
                for t_i = 1:t
                    scatter(x_traces(:,t_i),y_traces(:,t_i),[],colors_t(t_i,:),'.')
                end
                
                title(strcat("Timepoint ",num2str(t)));
                subtitle('Blue trks at this tp (input), red traces at this timepoint (interpolated)');
                axis equal
                
            else % 3d debug
                figure
                quiver3(pos(:,1),pos(:,2),pos(:,3),...
                    difference(:,1),difference(:,2), difference(:,3),...
                    0,'blue')
                hold on
                quiver3(x_traces(:,t),y_traces(:,t),z_traces(:,t),...
                    x_traces(:,t+1)-x_traces(:,t),...
                    y_traces(:,t+1)-y_traces(:,t),...
                    z_traces(:,t+1)-z_traces(:,t),...
                    0,'red');
                
                for t_i = max(1,t-1):t
                    scatter3(x_traces(:,t_i),y_traces(:,t_i),z_traces(:,t_i),...
                        [],colors_t(t_i,:),'.')
                end
                
                title(strcat("Timepoint ",num2str(t)));
                subtitle('Blue trks at this tp (input), red traces at this timepoint (interpolated) (can be on top of each other)');
                axis equal
                
            end
            
        end
        
    end
    
    %     % convert to d format as in track2disp() OLD
    %     if n_dim <3
    %         d(t).r = [x_traces(:,t),y_traces(:,t)];
    %         d(t).dr = [x_traces(:,t)-x_traces(:,1),y_traces(:,t)-y_traces(:,1)];
    %     else
    %         d(t).r = [x_traces(:,t),y_traces(:,t), z_traces(:,t)];
    %         d(t).dr = [x_traces(:,t)-x_traces(:,1),y_traces(:,t)-y_traces(:,1), z_traces(:,t)-z_traces(:,1)];
    %     end
    
    
end

if param.do_talk; fprintf('trks_to_displacement calculation done.\n'); end

% remove any trace that has nan in it
idx_nan = any(isnan(x_traces),2);
x_traces = x_traces(~idx_nan,:);
y_traces = y_traces(~idx_nan,:);
if n_dim>2; z_traces = z_traces(~idx_nan,:);end

if param.do_talk; fprintf('Removed %d traces because they were NaN (probably moved out of FOV).\n',sum(idx_nan)); end


% convert to track / trks
trks_out = [];
for tp=1:n_timepoints
    if n_dim>2
    trks_out = [trks_out; x_traces(:,tp), y_traces(:,tp),z_traces(:,tp),...
        unique_timepoints(tp)*ones(length(z_traces(:,tp)),1),... % timepoint
        [1:length(z_traces(:,tp))]'];
    else
            trks_out = [trks_out; x_traces(:,tp), y_traces(:,tp),...
        unique_timepoints(tp)*ones(length(x_traces(:,tp)),1),...
        [1:length(x_traces(:,tp))]'];
    end
end

d = track2disp(trks_out);



if param.do_save
    trks=trks_out;
    save(strcat(save_name,'_t2d.mat'),'x_traces','y_traces','d','trks');
    if param.do_talk; fprintf('Results saved.\n'); end
end


%% show results
if param.do_show
    
    if n_dim ==2
        f_arrows = figure; hold on
        for i=1:length(x_traces(:,1))
            plot(x_traces(i,:),y_traces(i,:),'-');
            
        end
        axis equal
        
        if true % param.do_debug
            figure; plot(x_traces,y_traces,'.'); title('colored by tp'); axis equal
            figure; plot(x_traces',y_traces','-.'); title('colored by particle id'); axis equal
        end
        
        
        % overall displacements
        ux_total = x_traces(:,end)-x_traces(:,1);
        uy_total = y_traces(:,end)-y_traces(:,1);
        %
        %     f_ux = figure;
        %     %imagesc(reshape(ux_total,sqrt(n_grid_pts),sqrt(n_grid_pts)))
        %     imagesc(reshape(ux_total,size(X_grid,1),size(X_grid,2)))
        
        % add overall displacement vector onto this to check the interpolation
        
        d_in = track2disp(trks);
        figure;
        imagesc(X_grid(:),Y_grid(:),reshape(ux_total,size(X_grid,1),size(X_grid,2)))
        hold on
        
        scatter(d_in(1).r(:,1),d_in(1).r(:,2),[],...
            d_in(end).dr(:,1),...
            'blacksquare', 'LineWidth',5,'MarkerFaceColor','black');
        scatter(d_in(1).r(:,1),d_in(1).r(:,2),[],...
            d_in(end).dr(:,1),'diamond','LineWidth',3);
        set(gca,'YDir','normal'); % same orientation as initial image
        colorbar
        title('u_x total');
        axis equal
        
        
        figure;
        imagesc(X_grid(:),Y_grid(:),reshape(uy_total,size(X_grid,1),size(X_grid,2)))
        hold on
        
        scatter(d_in(1).r(:,1),d_in(1).r(:,2),[],...
            d_in(end).dr(:,1),...
            'blacksquare', 'LineWidth',5,'MarkerFaceColor','black');
        scatter(d_in(1).r(:,1),d_in(1).r(:,2),[],...
            d_in(end).dr(:,2),'diamond','LineWidth',3);
        set(gca,'YDir','normal'); % same orientation as initial image
        colorbar
        title('u_y total');
        axis equal
        
        % saving
        if param.do_save
            export_fig(f_arrows,strcat(save_name,'_arrows.png'),...
                '-r300','-transparent')
            export_fig(f_ux, strcat(save_name,'_ux.png'),...
                '-r300','-transparent')
            export_fig(f_uy, strcat(save_name,'_uy.png'),...
                '-r300','-transparent')
            if param.do_talk; fprintf('Saved pictures in same folder. \n'); end
        end
    else % 3d show
        figure; viz_tracks(trks_out);
    end
    
    if param.do_debug; autoArrangeFigures(); end
    
end

