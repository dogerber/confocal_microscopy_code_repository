function res = displacements_to_principal_stretches(PATH_DISPLACEMENTS,param)
%  res = displacements_to_principal_stretches(PATH_DISPLACEMENTS,param)
% Converts from the displacements d().r format to principal stretches using the calcF
% function
%
%
% Branched from displacements_to_strains 10.11.2022 by DG
%
%
% see also: https://en.wikipedia.org/wiki/Finite_strain_theory
% principal_stretch_to_stresses(),
%
%
%
%
% --------------------------------------------- Dominic Gerber, 2022


% Options
if ~exist('param','var');         param = struct();           end
if ~isfield(param,'do_debug');    param.do_debug =  true;     end
if ~isfield(param,'do_talk');     param.do_talk = true;       end
if ~isfield(param,'do_show');     param.do_show = false;       end
if ~isfield(param,'do_save');     param.do_save = false;      end
if ~isfield(param,'do_waitbar');  param.do_waitbar = true;    end

if isunix; param.do_waitbar = false; param.do_show = false;   end

% strain calucaltion parameters
if ~isfield(param,'neighbour_mode'); param.neighbour_mode = 'distance'; end% 'number' or 'distance'
if ~isfield(param,'neighbour_param');param.neighbour_param = 10; end% number of neighbours to use, or neighbours within this distance to use
if ~isfield(param,'target_state');param.target_state = 'undeformed'; end% for which state do you want to calculate the strains? 'deformed' or 'undeformed'
if ~isfield(param,'use_backslash');param.use_backslash = false; end % use fast and sensitive backsash mthod
param.code_version = 1001; %

if ~isfield(param,'dimensions_to_use');  param.dimensions_to_use = 'xz'; end % reduce to 2d or keep 3d data
% options: 'xyz', 'xy', 'xz'

%% Data import

% Quick and dirty way to get d() from Nicos code:
% d(1).r = r_ref_tracked; d(1).dr = r_ref_tracked*0; d(2).r = r_def_tracked; d(2).dr = r_def_tracked-r_ref_tracked;


if ~exist('PATH_DISPLACEMENTS','var')
    [ffname,ffolder] = uigetfile('Choose the folder in which the slices are','*.mat*');
    PATH_DISPLACEMENTS = fullfile(ffolder, ffname);
    param.fullfilename_in = PATH_DISPLACEMENTS;
end

if isstruct(PATH_DISPLACEMENTS) % input is d()
    d = PATH_DISPLACEMENTS;
else
    [ffolder,ffname] = fileparts(PATH_DISPLACEMENTS); % input is path
    load(PATH_DISPLACEMENTS,'d'); % dispalcements structure e.g. from confocal_to_displacments.m
end

if false
    warning('temporary real space correction for Nicos data is ON');
    muperpx = 0.329131;
    pause;
    for tp=1:length(d)
        d(tp).r = d(tp).r .*[muperpx,muperpx,1];
        d(tp).dr = d(tp).dr .*[muperpx,muperpx,1];
    end
end


% fake data for testing
if false
    warning('FAKING DATA'); pause;
    d = d(1);
    Ef = [1 0.15/2 0.1;...
        0.15/2 1 0.3/2;...
        0.1 0.3/2 1.1] ;
    d(2).r = d(1).r(:,:) * Ef.*d(1).r(:,3)/max(d(1).r(:,3));
    d(2).dr = d(2).r-d(1).r;
    
end


% Debug: show input
if param.do_debug
    figure;
    if length(d(1).r(1,:))==3
        plot3(d(1).r(:,1),d(1).r(:,2),d(1).r(:,3),'r.','DisplayName','reference');
        hold on; plot3(d(end).r(:,1),d(end).r(:,2),d(end).r(:,3),'b.','DisplayName','deformed');
    else
        plot(d(1).r(:,1),d(1).r(:,2),'r.','DisplayName','reference');
        hold on; plot(d(end).r(:,1),d(end).r(:,2),'b.','DisplayName','max. deformed');
    end
    title('Data input RAW'); legend;
end


% Convert everything to meters
if false
    for tp=1:length(d)
        for dim = 1:length(d(1).r(1,:))
            d(tp).r(:,dim) = d(tp).r(:,dim)/10^6;
            d(tp).dr(:,dim) = d(tp).dr(:,dim)/10^6;
        end
    end
    warning('conversion to SI units turned ON!');
else
    if param.do_debug;  warning('conversion to SI units turned OFF!'); end
end


%% Convert to stretches with calcF (local)
r_start = d(1).r;
%wrong? r_start(:,3) = r_start(:,3)*0; % z=0 for reference state
r_end = d(end).r;
%wrong? r_end(:,3) = d(end).dr(:,3); % relative to first tp

% if param.do_debug
%     figure;
%     plot3(r_start(:,1),r_start(:,2),r_start(:,3),'r.','DisplayName','reference');
%     hold on; plot3(r_end(:,1),r_end(:,2),r_end(:,3),'b.','DisplayName','deformed');
%     title('Data input'); legend;
% end

switch param.dimensions_to_use
    case 'xyz' % full 3d
        if false
            [X_grid,Y_grid,Z_grid] = create_grid_v2([d(1).r(:,1),d(1).r(:,2),d(1).r(:,3)],true,25e-6);
            warning('I am oversampling when making the grid');
            res = local_calcF_for_principal_stretches(r_start,r_end,[X_grid(:),Y_grid(:),Z_grid(:)],param);
        else
            warning('Not using a grid, taking the r_start points as locations to calculate stretches on');
            res = local_calcF_for_principal_stretches(r_start,r_end,r_start,param);
        end
        n_dim = 3;
        
    case 'xy' % 2d, looking from z
        if false
        [X_grid,Y_grid] = create_grid_v2([d(1).r(:,1),d(1).r(:,2)]);
        res = local_calcF_for_principal_stretches(r_start(:,1:2),r_end(:,1:2),[X_grid(:),Y_grid(:)],param);
        n_dim = 2;
           else
                warning('Not using a grid, taking the r_start points as locations to calculate stretches on');
            res = local_calcF_for_principal_stretches(r_start(:,[1,2]),r_end(:,[1,2]),r_start(:,[1,2]),param);
        end
    case 'xz' % 2d, looking from y-direction
        if false
        [X_grid,Z_grid] = create_grid_v2([d(1).r(:,1),d(1).r(:,3)]);
        res = local_calcF_for_principal_stretches(r_start(:,[1,3]),r_end(:,[1,3]),[X_grid(:),Z_grid(:)],param);
        else
                warning('Not using a grid, taking the r_start points as locations to calculate stretches on');
            res = local_calcF_for_principal_stretches(r_start(:,[1,3]),r_end(:,[1,3]),r_start(:,[1,3]),param);
        end
        n_dim = 2;
    case 'yz'
        [Y_grid,Z_grid] = create_grid_v2([d(1).r(:,1),d(2).r(:,3)]);
        res = local_calcF_for_principal_stretches(r_start(:,[2,3]),r_end(:,[2,3]),[Y_grid(:),Z_grid(:)],param);
        n_dim = 2;
end



%% Visualize and clean up results
% for plotting strains epsilon, referr to displacements_to_strains.m
if param.do_show

if n_dim==2 % 2D
    
    
viz_strains(res);
    
else % 3D data =========================================================================
    
    % ----- show input
    if true
        % show input data and grid
        figure; hold on
        plot3(res.r_start(:,1),res.r_start(:,2),res.r_start(:,3),...
            'bo');
        plot3(res.r_end(:,1),res.r_end(:,2),res.r_end(:,3),...
            'rx');
        plot3(res.X_grid(:,1),res.X_grid(:,2),res.X_grid(:,3),...
            'black.');
        title('Input data and grid points');
        axis equal
        view(45,45);
        
        % as a reference show displacements in 2d
        figure; hold on
        quiver(res.r_start(:,1),res.r_start(:,3),...
            res.r_end(:,1)-res.r_start(:,1),...
            res.r_end(:,3)-res.r_start(:,3),...
            0);
        grid on
        axis equal
        xlabel('x'); ylabel('z'); title('Input data as quivers 2D');
        
        % in 3d
        figure; hold on
        quiver3(res.r_start(:,1),res.r_start(:,2),res.r_start(:,3),...
            res.r_end(:,1)-res.r_start(:,1),...
            res.r_end(:,2)-res.r_start(:,2),...
            res.r_end(:,3)-res.r_start(:,3),...
            0);
        grid on; axis equal
        xlabel('x'); ylabel('y');
        title('Input data as quivers 3D');
    end
    
    
    % ----- show strains (as a check)
    if true
        
        figure('Name','Strain'); hold on
        subplot(2,2,1); hold on
        scatter3(res.X_grid(:,1),res.X_grid(:,2),res.X_grid(:,3),...
            [],res.epsilon_xx(:),'.');
        colorbar; view([45,45]);
        caxis([quantile(res.epsilon_xx,0.01),quantile(res.epsilon_xx,0.99)]);
        subtitle('\epsilon_{xx}'); xlabel('x'); ylabel('y'); grid on;
        
        local_styler_symmetrical()
        
        subplot(2,2,2); hold on
        scatter3(res.X_grid(:,1),res.X_grid(:,2),res.X_grid(:,3),...
            [],res.epsilon_yy(:),'.');
        colorbar; view([45,45]);
        caxis([quantile(res.epsilon_xx,0.01),quantile(res.epsilon_xx,0.99)]);
        subtitle('\epsilon_{yy}');xlabel('x'); ylabel('y'); grid on;
        
        local_styler_symmetrical()
        
        subplot(2,2,3); hold on
        scatter3(res.X_grid(:,1),res.X_grid(:,2),res.X_grid(:,3),...
            [],res.epsilon_zz(:),'.');
        colorbar; view([45,45]);
        caxis([quantile(res.epsilon_xx,0.01),quantile(res.epsilon_xx,0.99)]);
        local_styler_symmetrical()
        subtitle('\epsilon_{zz}');xlabel('x'); ylabel('y'); grid on;
        
        
    end
    
    
    
    % ---- show principal stretches
    if true
        % reform to arrays to speed up quiver3
        eigen_values_stretch =  nan(length(res.principal_stretches),3,3);
        eigen_vectors = nan(length(res.principal_stretches),3,3);
        
        for i=1:length(res.principal_stretches)
            try
                eigen_vectors(i,:,:) = [res.principal_stretches(i).vectors];
                eigen_values_stretch(i,:,:) = [res.principal_stretches(i).values];
            end
        end
        
        % show principal stretches eigenvalues
        if true
            figure
            for i=1:3
                subplot(1,3,i);
                scatter3(res.r_end(:,1),res.r_end(:,2),res.r_end(:,3),[],...
                    squeeze(eigen_values_stretch(:,i,i)),...
                    '.');
                hold on
              %  plot3(res.r_end(:,1),res.r_end(:,2),res.r_end(:,3),'o','Color',0.3*ones(3,1),'MarkerSize',3); % to give outline of particles
                colorbar('Location','south');
                
                if true
                    plot_softliv_style; %set(gca,'Color',0.3*ones(3,1));
                    max_diff_abs = max(abs(1-caxis));
                    caxis(1+[-max_diff_abs,max_diff_abs])
                   % caxis([0.9,1.1])
                    colormap_ice
                   % axis equal
                end
            end
            title('Magnitude of the three eigenvalues of stretch (principal stretches)');
            
            %             figure;
            %             for i=1:3
            %                 subplot(3,1,i)
            %                 histogram(L_store(:,i));
            %                 title(strcat('Principal stretch \lambda_',num2str(i)));
            %             end
            
        end
        
        % quiver3 of principal stretches
        figure; hold on
        subplot(1,2,1); hold on
        for i=1:3
            quiver3(res.r_end(:,1),res.r_end(:,2),res.r_end(:,3),...
                eigen_vectors(:,1,i).*eigen_values_stretch(:,i,i), ...
                eigen_vectors(:,2,i).*eigen_values_stretch(:,i,i), ...
                eigen_vectors(:,3,i).*eigen_values_stretch(:,i,i), ...
                1);
        end
        title('Prinicpal stretches * eigenvalues'); subtitle('autoscaled'); view([-45,45]);
        grid on;
        
        subplot(1,2,2); hold on
        for i=1:3
            quiver3(res.r_end(:,1),res.r_end(:,2),res.r_end(:,3),...
                eigen_vectors(:,1,i), ...
                eigen_vectors(:,2,i), ...
                eigen_vectors(:,3,i));
        end
        title('Only eigenvectors'); subtitle('autoscaled'); view([-45,45]); grid on;
        
    end
    
    % ---- Debugging visualizations
    % show one group of traces at a time (debugging results)
    if false
        figure;
        for idx_grid = randperm(length(res.X_grid))
            inds_n_idx_grid = res.all_inds(idx_grid).inds_n;
            subplot(1,2,1); cla; hold on; % show data selection
            scatter3(res.r_end(:,1),res.r_end(:,2),res.r_end(:,3),[],res.r_end(:,3),'.');
            plot3(res.r_start(inds_n_idx_grid,1),res.r_start(inds_n_idx_grid,2),res.r_start(inds_n_idx_grid,3),'blueo');
            plot3(res.r_end(inds_n_idx_grid,1),res.r_end(inds_n_idx_grid,2),res.r_end(inds_n_idx_grid,3),'redo');
            title('data selection'); view([45,45]); grid on
            
            subplot(1,2,2); cla; hold on;
            if length(res.all_inds(idx_grid).inds_n) > 5
                
                
                % show grid, start and end points for this group
                plot3(res.X_grid(idx_grid,1),res.X_grid(idx_grid,2),res.X_grid(idx_grid,3),'black.'); % (can be same as start_pts)
                plot3(res.r_start(inds_n_idx_grid,1),res.r_start(inds_n_idx_grid,2),res.r_start(inds_n_idx_grid,3),'blue.');
                plot3(res.r_end(inds_n_idx_grid,1),res.r_end(inds_n_idx_grid,2),res.r_end(inds_n_idx_grid,3),'red.');
                
                for i=1:n_dim
                    % principal stretch arrow
                    quiver3(res.X_grid(idx_grid,1),res.X_grid(idx_grid,2),res.X_grid(idx_grid,3),...
                        res.principal_stretches(idx_grid).vectors(1,i)*res.principal_stretches(idx_grid).values(i,i), ...
                        res.principal_stretches(idx_grid).vectors(2,i)*res.principal_stretches(idx_grid).values(i,i),...
                        res.principal_stretches(idx_grid).vectors(3,i)*res.principal_stretches(idx_grid).values(i,i),...
                        1*10^(-5),'red');
                end
                
                % displacement arrows
                quiver3(res.r_start(inds_n_idx_grid,1),res.r_start(inds_n_idx_grid,2),res.r_start(inds_n_idx_grid,3),...
                    res.r_end(inds_n_idx_grid,1)-res.r_start(inds_n_idx_grid,1)-1*mean(res.r_end(inds_n_idx_grid,1)-res.r_start(inds_n_idx_grid,1)),...
                    res.r_end(inds_n_idx_grid,2)-res.r_start(inds_n_idx_grid,2)-1*mean(res.r_end(inds_n_idx_grid,2)-res.r_start(inds_n_idx_grid,2)),...
                    res.r_end(inds_n_idx_grid,3)-res.r_start(inds_n_idx_grid,3)-1*mean(res.r_end(inds_n_idx_grid,3)-res.r_start(inds_n_idx_grid,3)),...
                    0,'Color',ones(3,1)*0.4); % gray color
                
                axis equal; view(45,45); grid on
                title('Principal stretches and the displacements used for it');
                shg
                pause;
            end % if
        end % for idx_grid
        axis equal
    end
    
    
    % for debugging show groups of inds and what happens to them
    if false
        figure; hold on
        for idx_grid = 1:200:length(res.X_grid)
            if length(res.all_inds(idx_grid).inds_n) > 5
                inds_n_idx_grid = res.all_inds(idx_grid).inds_n;
                plot3(res.X_grid(idx_grid,1),res.X_grid(idx_grid,2),res.X_grid(idx_grid,3),'black.');
                
                plot3(res.r_start(inds_n_idx_grid,1),res.r_start(inds_n_idx_grid,2),res.r_start(inds_n_idx_grid,3),'blue.');
                plot3(res.r_end(inds_n_idx_grid,1),res.r_end(inds_n_idx_grid,2),res.r_end(inds_n_idx_grid,3),'red.');
                
                
                for i=1:n_dim
                    try
                        quiver3(res.X_grid(idx_grid,1),res.X_grid(idx_grid,2),res.X_grid(idx_grid,3),...
                            res.principal_stretches(idx_grid).vectors(1,i)*res.principal_stretches(idx_grid).values(i,i), ...
                            res.principal_stretches(idx_grid).vectors(2,i)*res.principal_stretches(idx_grid).values(i,i),...
                            res.principal_stretches(idx_grid).vectors(3,i)*res.principal_stretches(idx_grid).values(i,i),...
                            1*10^(-5));
                    end
                end
                axis equal; view(45,45); grid on
                %             shg
                %             pause;
            end % if
        end % for idx_grid
        axis equal
    end
    

    
end
end % do_show



end % function


        
        
%% Sub-Functions


function res = local_calcF_for_principal_stretches(X,x,X_grid,param)
% Code to take two sets of particle positions, with an index that links the
% first set to the second, and calculate the local deformation gradient,
% the local principle stretches

% see also https://en.wikipedia.org/wiki/Finite_strain_theory


% X is a three column matrix with the X,Y,Z coordinates of the
% pre-deformation data set
% x is a three column matrix with the x,y,z coordinates of the
% post-deformation data set
% X_grid is the two column matrix with the interpolation points from the
% reference frame
% n is the number of nearest neighbours that you want to use to calculate
% the strain

if ~isfield(param,'target_state'); param.target_state = 'deformed'; end

res = struct();
param.date_done = datetime('now');
param.done_by_code = 'local_calcF_for_principal_stretches in displacements_to_principal_stretches';

if ~isfield(param,'neighbour_mode'); param.neighbour_mode = 'number'; end % 'number' or 'distance'
if ~isfield(param,'neighbour_param'); param.neighbour_param = 20; end

% note down options
res.param = param; 

n_dim = size(X,2);

%% allocate space
epsilon_xx = zeros(length(X_grid),1);
epsilon_yy  = epsilon_xx;
epsilon_xy  = epsilon_xx;
principal_stretches = struct();

if n_dim>2 % 3D
    epsilon_zz = epsilon_xx; epsilon_xz =epsilon_xx;
    epsilon_yz = epsilon_xx;  epsilon_kk = epsilon_xx;
end


if param.do_waitbar
    wb = waitbar(0,'calcF @ displacements to principal stretches');
end

%% find neighbours 
min_num_pts_for_deformation_tensor = n_dim+1;

if strcmp(param.neighbour_mode,'number') % old method
    % fprintf('Deciding on neighbours by fixed number of neighbours\n');
    inds_n_all = knnsearch(X,X_grid,'K',param.neighbour_param); % can make this faster by putting outside the loop
elseif strcmp(param.neighbour_mode,'distance')
    % fprintf('Deciding on neighbours by distance\n');
    % Do by distance instead and make NAN if not enough neighbours?
    if false
        dx_guess = abs(X_grid(1,1)-X_grid(1,2)); % this is wrong if data is not gridded
        neighbour_range = 20*dx_guess;
    else
        warning('manual range for neighbours decision');
        neighbour_range = 10;
    end
    inds_n_all = rangesearch(X,X_grid,neighbour_range);
 
    
    if param.do_debug % show num_neighbour distribution
        num_neighbours = cellfun(@length,inds_n_all);
        figure; hold on
        histogram(num_neighbours);
        xline(min_num_pts_for_deformation_tensor);
        xlabel('number of neighbours in range'); ylabel('count');
        title('number of neighbours per X_grid','Interpreter','none');
        subtitle(sprintf('with range %d, line for given minimum amount',neighbour_range));
    end
end

fprintf('Turning iterations limit warnings off!\n')
warning('off','stats:statrobustfit:IterationLimit')

if param.do_debug
    fprintf('Check now all messages and warnings and press ENTER, before the console gets spammed.\n');
    pause;
end


%% to calculate local deformation gradient
for i=1:length(X_grid)

    % uses inds_n from rangesearch above
    try
        inds_n = inds_n_all{i};
    catch
        inds_n = inds_n_all(i,:);
    end
    n = length(inds_n);
    
    if length(inds_n) > min_num_pts_for_deformation_tensor
        
        dX=X(inds_n,:)-repmat(X_grid(i,:),n,1); % The matrix of distances from the gridpoint
        dx=x(inds_n,:)-repmat(X_grid(i,:),n,1);
        
        % for debugging
        if false
            figure('Position',[200,0,900,1200]);
            subplot(2,1,1)
            plot(X_grid(:,1),X_grid(:,2),'black.');
            hold on
            plot(X_grid(i,1),X_grid(i,2),'ro');
            plot(X(:,1),X(:,2),'.','Color',[0.5,0.5,0.5]);
            plot(X(inds_n,1),X(inds_n,2),'bo');
            plot(x(inds_n,1),x(inds_n,2),'go');
            subtitle('Black grid total, red current grid point, grey all X, blue selected X, green deformed selected x');
            axis image
            
            subplot(2,1,2);
            plot(dX(:,1),dX(:,2),'b.'); hold on; plot(dx(:,1),dx(:,2),'g.');
            title('relative coordinates');
            pause
        end
        
        % -- Determine F
        if param.use_backslash % linear fit with matlab backslash (susceptible to outliers)
            % Backslash operation to get gradient Tensor
            %warning('Backslash method active.');
            if true % linear fit
                D=[ones(n,1),dX]; %
            else % quadratic fit
                D = [ones(n,1),dX, dX.^2, dX(:,1).*dX(:,2), dX(:,1).*dX(:,3), dX(:,2).*dX(:,3)];
            end
            A=D\dx; % Matlab-magic fitting with backslash
            
            % Derive all the tensor terms from A
            x0(i,1:n_dim)=A(1,:); % Translation terms
            F = A(2:n_dim+1,:); % "Deformation gradient"
            
        else % more robust fitting with fitlm (slower)
            A = nan(n_dim+1,n_dim+1);
            try
                if size(dX,1)>n_dim+2 % needed for robust fit
                    for idx_dim = 1:n_dim % for each dimension find the coefficients to predict it from the fulldimension input
                        mdllr = fitlm( dX, dx(:,idx_dim),'RobustOpts','on');
                        A(idx_dim, :) = mdllr.Coefficients.Estimate([2:n_dim+1,1])';
                    end
                    
                    % Derive all the tensor terms from A
                    x0(i,1:n_dim)=A(1:n_dim,end); % Translation terms
                    F = A(1:n_dim,1:n_dim); % "Deformation gradient"
                    
                    if false % debugging
                        figure; hold on
                        plot(dX(:,1),dX(:,2),'blue.','DisplayName','origin');
                        plot(dx(:,1),dx(:,2),'red.','DisplayName','target');
                        %                       plot3(dX(:,1),dX(:,2),dX(:,3),'blue.','DisplayName','origin');
                        %                         plot3(dx(:,1),dx(:,2),dx(:,3),'red.','DisplayName','target');
                        est_dx =A*[dX, ones(length(dX),1)]';
                        est_dx = est_dx(1:3,:)';
                        plot(est_dx(:,1),est_dx(:,2),'redo','DisplayName','est_target');
                        plot3(est_dx(:,1),est_dx(:,2),est_dx(:,3),'redo','DisplayName','est_target');
                        shg; grid on; legend; view([45,45]);
                        pause;
                    end
                end
            catch
                % beep;
            end
        end
        
        % convert F to other tensors, determine eigenvalues
        if exist('F','var')
            C = F'*F; % right Cauchy–Green deformation tensor
            B = F*F'; % "Finger tensor" or left Cauchy-Green tensor
            %wrong? E = (F+ F')/2 - eye(n_dim); %old from Rob, "strain tensor" see https://en.wikipedia.org/wiki/Finite_strain_theory#Finite_strain_tensors
            E = (C - eye(n_dim))/2;  % Stefanie 230525, "Lagrangian finite strain tensor, also called the Green-Lagrangian strain", for undeformed configuration
            e = (eye(n_dim)-inv(B))/2; % Eulerian-Alamnsi finite strain tensor, for deformed configuration

            
            if strcmp(param.target_state,'deformed')
                [principal_stretches(i).vectors, principal_stretches(i).values] = eig(B); % B in deformed state (eigenvalues are the same)
            elseif strcmp(param.target_state,'undeformed')
                [principal_stretches(i).vectors, principal_stretches(i).values] = eig(C); % C in reference state (eigenvalues are the same)
            end
             principal_stretches(i).values = sqrt( principal_stretches(i).values); % principal stretch are the square root of the eigenvalues, Stefanie 230525
        else
            E = nan(3,3);
            principal_stretches(i).vectors = nan(n_dim,n_dim);
            principal_stretches(i).values = nan(n_dim,n_dim);
        end
        
    else % did not fin enough neighbours
        E = nan(3,3);
        principal_stretches(i).vectors = nan(n_dim,n_dim);
        principal_stretches(i).values = nan(n_dim,n_dim);
    end
    
    % note down
    res.all_inds(i).inds_n = inds_n;
    
    if strcmp(param.target_state,'undeformed')
        strain_matrix = E; % Green-Lagrange strain, for undeformed configuration
    elseif strcmp(param.target_state,'deformed')
        strain_matrix = e; % Eulerian-Alamnsi finite strain tensor, for deformed configuration
    end
    
    epsilon_xx(i) = strain_matrix(1,1); 
    epsilon_yy(i) = strain_matrix(2,2);
    epsilon_xy(i) = strain_matrix(1,2);
    
    epsilon_kk(i) = principal_stretches(i).values(1,1)*principal_stretches(i).values(2,2);
    
    if n_dim>2 % 3D
        epsilon_zz(i) = strain_matrix(3,3);
        epsilon_xz(i) = strain_matrix(1,3);
        epsilon_yz(i) = strain_matrix(2,3);
        
        % wrong, changed 230608 epsilon_kk(i) = strain_matrix(1,1)+strain_matrix(2,2)+strain_matrix(3,3);
          epsilon_kk(i) = principal_stretches(i).values(1,1)*principal_stretches(i).values(2,2)*principal_stretches(i).values(3,3);
    end
    
    if param.do_waitbar
        waitbar(i/length(X_grid),wb);
    end
end % for i

if param.do_waitbar; close(wb); end

%% note down results for output
res.r_start = X;
res.r_end = x;
res.X_grid = X_grid;
res.translation=x0;
res.epsilon_xx=epsilon_xx;
res.epsilon_yy=epsilon_yy;
res.epsilon_xy=epsilon_xy;
res.epsilon_kk  = epsilon_kk;

res.principal_stretches = principal_stretches;

if n_dim>2 % 3D
    res.epsilon_zz  = epsilon_zz;
    res.epsilon_xz  = epsilon_xz;
    res.epsilon_yz  = epsilon_yz;

end


end % local_calcF...





function local_styler_symmetrical()
     plot_softliv_style; %set(gca,'Color',0.3*ones(3,1));
                    max_diff_abs = max(abs(caxis));
                    caxis([-max_diff_abs,max_diff_abs])
                    colormap_ice
end

