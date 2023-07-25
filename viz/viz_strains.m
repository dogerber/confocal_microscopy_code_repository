function viz_strains(res)
%{



%}

%% Options

param.do_show_input = true;
gl_caxis = [-1,1]*0.8;


%% Look at data
n_dim = size(res.r_start,2);


%% Viz
if n_dim==2 % 2D
    
    
    if param.do_show_input
        % show input data and grid
        figure; hold on
        plot(res.r_start(:,1),res.r_start(:,2),...
            'bo');
        plot(res.r_end(:,1),res.r_end(:,2),...
            'rx');
        plot(res.X_grid(:,1),res.X_grid(:,2),...
            'black.');
        
        figure;
        diff_r = res.r_end-res.r_start;
        quiver(res.r_start(:,1),res.r_start(:,2),...
            diff_r(:,1),diff_r(:,2),0);
    end
    
    % - normal epsilon
    if false % on undeformed configuration
        figure('Name','Strain'); hold on
        subplot(2,2,1); hold on
        scatter(res.X_grid(:,1),res.X_grid(:,2),...
            [],res.epsilon_xx(:),'.');
        colorbar;
        caxis([quantile(res.epsilon_xx,0.01),quantile(res.epsilon_xx,0.99)]);
        subtitle('\epsilon_{xx}'); xlabel('x'); ylabel('y'); grid on;
        
        %  local_styler_symmetrical()
        
        subplot(2,2,2); hold on
        scatter(res.X_grid(:,1),res.X_grid(:,2),...
            [],res.epsilon_yy(:),'.');
        colorbar;
        caxis([quantile(res.epsilon_xx,0.01),quantile(res.epsilon_xx,0.99)]);
        subtitle('\epsilon_{yy}');xlabel('x'); ylabel('y'); grid on;
        
        % local_styler_symmetrical()
        
        subplot(2,2,3); hold on
        scatter(res.X_grid(:,1),res.X_grid(:,2),...
            [],res.epsilon_xy(:),'.');
        colorbar;
        local_styler_symmetrical()
        subtitle('\epsilon_{xy}');xlabel('x'); ylabel('y'); grid on;
    else % on deformed configuration
        figure('Name','Strain'); hold on
        subplot(2,2,1); hold on
        scatter(res.r_end(:,1),res.r_end(:,2),...
            [],res.epsilon_xx(:),'.');
        colorbar;
        if isempty(gl_caxis)
            caxis( [-1,1]*max(abs(caxis)));
        else
        caxis(gl_caxis);
        end
        colormap_ice;
        subtitle('\epsilon_{xx}'); xlabel('x'); ylabel('y'); grid on;

        subplot(2,2,2); hold on
        scatter(res.r_end(:,1),res.r_end(:,2),...
            [],res.epsilon_yy(:),'.');
        colorbar;
            if isempty(gl_caxis)
            caxis( [-1,1]*max(abs(caxis)));
        else
        caxis(gl_caxis);
        end
        colormap_ice;
        subtitle('\epsilon_{yy}');xlabel('x'); ylabel('y'); grid on;
        subplot(2,2,3); hold on
        scatter(res.r_end(:,1),res.r_end(:,2),...
            [],res.epsilon_xy(:),'.');
        colorbar;
            if isempty(gl_caxis)
            caxis( [-1,1]*max(abs(caxis)));
        else
        caxis(gl_caxis);
        end
        colormap_ice;
        subtitle('\epsilon_{xy}');xlabel('x'); ylabel('y'); grid on;
    end
    
    % strain with axis equal
    if true
            figure('Name','Strain, axis equal'); hold on
        subplot(2,2,1); hold on
        scatter(res.r_end(:,1),res.r_end(:,2),...
            [],res.epsilon_xx(:),'.');
        colorbar;
        axis equal;
            if isempty(gl_caxis)
            caxis( [-1,1]*max(abs(caxis)));
        else
        caxis(gl_caxis);
        end
        colormap_ice;
        subtitle('\epsilon_{xx}'); xlabel('x'); ylabel('y'); grid on;

        subplot(2,2,2); hold on
        scatter(res.r_end(:,1),res.r_end(:,2),...
            [],res.epsilon_yy(:),'.');
        colorbar;
        if isempty(gl_caxis)
            caxis( [-1,1]*max(abs(caxis)));
        else
        caxis(gl_caxis);
        end
        colormap_ice;
        axis equal;
        subtitle('\epsilon_{yy}');xlabel('x'); ylabel('y'); grid on;
        subplot(2,2,3); hold on
        scatter(res.r_end(:,1),res.r_end(:,2),...
            [],res.epsilon_xy(:),'.');
        colorbar;
            if isempty(gl_caxis)
            caxis( [-1,1]*max(abs(caxis)));
        else
        caxis(gl_caxis);
        end
        colormap_ice;
        axis equal;
        subtitle('\epsilon_{xy}');xlabel('x'); ylabel('y'); grid on;
    end
    
    
    % -- show principal stretches
    eigen_values_stretch =  nan(length(res.principal_stretches),2,2);
    eigen_vectors = nan(length(res.principal_stretches),2,2);
    
    % reshape to arrays to make quiver faster
    for i=1:length(res.principal_stretches)
        try
            eigen_vectors(i,:,:) = [res.principal_stretches(i).vectors];
            eigen_values_stretch(i,:,:) = [res.principal_stretches(i).values];
        end
    end
    
    % plot quiver
    if false
        figure; hold on
        for i=1:2
            quiver (res.X_grid(:,1),res.X_grid(:,2),...
                eigen_vectors(:,1,i).*eigen_values_stretch(:,i,i), ...
                eigen_vectors(:,2,i).*eigen_values_stretch(:,i,i), ...
                1);
        end
        title('Prinicpal stretches * eigenvalues'); subtitle('on undeformed coordinates, autoscaled');
    end
    
    figure; hold on
    for i=1:2
        quiver (res.r_end(:,1),res.r_end(:,2),...
            eigen_vectors(:,1,i).*eigen_values_stretch(:,i,i), ...
            eigen_vectors(:,2,i).*eigen_values_stretch(:,i,i), ...
            1);
    end
    title('Prinicpal stretches * eigenvalues'); subtitle('on deformed coordinates, autoscaled');
    axis equal
    
    
    % principal stretch components
    figure('Name','Prinicpal stretches'); hold on
    subplot(2,1,1); hold on
    scatter(res.r_end(:,1),res.r_end(:,2),...
        [],eigen_values_stretch(:,1,1),'.');
    colorbar;
    %    caxis([quantile(res.epsilon_xx,0.01),quantile(res.epsilon_xx,0.99)]);
    subtitle('\epsilon_{1}'); xlabel('x'); ylabel('y'); grid on;
    colormap_ice; caxis([0,2])
    
    subplot(2,1,2); hold on
    scatter(res.r_end(:,1),res.r_end(:,2),...
        [],eigen_values_stretch(:,2,2),'.');
    colorbar;
    %    caxis([quantile(res.epsilon_xx,0.01),quantile(res.epsilon_xx,0.99)]);
    subtitle('\epsilon_{2}'); xlabel('x'); ylabel('y'); grid on;
    colormap_ice; caxis([0,2])
    
    
    figure('Name','Prinicpal stretches - volumetric'); hold on
    scatter(res.r_end(:,1),res.r_end(:,2),...
        [],eigen_values_stretch(:,1,1).*eigen_values_stretch(:,2,2),'.');
    colorbar;
    subtitle('\epsilon_{volume}'); xlabel('x'); ylabel('y'); grid on;
    colormap_ice; caxis([0,2])
    
    
    
    % plot principal strain vs position
    if false
        figure; hold on
        plot(res.r_end(:,1),eigen_values_stretch(:,1,1),'.');
        plot(res.r_end(:,1),eigen_values_stretch(:,2,2),'.');
        
        
        figure; hold on
        gradient = 1.0763/10^3;
        x_temp_T = (res.r_end(:,1)*(-1)+200)*gradient;
        hold on
        plot(x_temp_T,eigen_values_stretch(:,1,1),'.');
        plot(x_temp_T,eigen_values_stretch(:,2,2),'.');
        plot(x_temp_T,eigen_values_stretch(:,1,1).*eigen_values_stretch(:,2,2),'.');
        title('volumetric strain vs position');
        grid on
        xlabel('check this');
    end
    
else % 3D data =========================================================================
    
    % ----- show input
    if false
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
        view_angle = [-45,45]; 
        
     figure('Name','Strain'); hold on
        subplot(2,3,1); hold on
        scatter3(res.r_end(:,1),res.r_end(:,2),res.r_end(:,3),...
            [],res.epsilon_xx(:),'.');
        colorbar;
        if isempty(gl_caxis)
            caxis( [-1,1]*max(abs(caxis)));
        else
            caxis(gl_caxis);
        end
        colormap_ice;
        subtitle('\epsilon_{xx}'); xlabel('x'); ylabel('y'); grid on;
             view(view_angle);
        
        subplot(2,3,2); hold on
        scatter3(res.r_end(:,1),res.r_end(:,2),res.r_end(:,3),...
            [],res.epsilon_yy(:),'.');
        colorbar;
        if isempty(gl_caxis)
            caxis( [-1,1]*max(abs(caxis)));
        else
            caxis(gl_caxis);
        end
        colormap_ice;
        subtitle('\epsilon_{yy}');xlabel('x'); ylabel('y'); grid on;
             view(view_angle);
        
        subplot(2,3,3); hold on
        scatter3(res.r_end(:,1),res.r_end(:,2),res.r_end(:,3),...
            [],res.epsilon_zz(:),'.');
        colorbar;
        if isempty(gl_caxis)
            caxis( [-1,1]*max(abs(caxis)));
        else
            caxis(gl_caxis);
        end
        colormap_ice;
        subtitle('\epsilon_{zz}');xlabel('x'); ylabel('y'); grid on;
             view(view_angle);
        
                subplot(2,3,4); hold on
        scatter3(res.r_end(:,1),res.r_end(:,2),res.r_end(:,3),...
            [],res.epsilon_xy(:),'.');
        colorbar;
        if isempty(gl_caxis)
            caxis( [-1,1]*max(abs(caxis)));
        else
            caxis(gl_caxis);
        end
        colormap_ice;
        subtitle('\epsilon_{xy}');xlabel('x'); ylabel('y'); grid on;
             view(view_angle);
        
                        subplot(2,3,5); hold on
        scatter3(res.r_end(:,1),res.r_end(:,2),res.r_end(:,3),...
            [],res.epsilon_xz(:),'.');
        colorbar;
        if isempty(gl_caxis)
            caxis( [-1,1]*max(abs(caxis)));
        else
            caxis(gl_caxis);
        end
        colormap_ice;
        subtitle('\epsilon_{xz}');xlabel('x'); ylabel('y'); grid on;
             view(view_angle);
        
                        subplot(2,3,6); hold on
        scatter3(res.r_end(:,1),res.r_end(:,2),res.r_end(:,3),...
            [],res.epsilon_yz(:),'.');
        colorbar;
        if isempty(gl_caxis)
            caxis( [-1,1]*max(abs(caxis)));
        else
            caxis(gl_caxis);
        end
        colormap_ice;
        subtitle('\epsilon_{yz}');xlabel('x'); ylabel('y'); grid on;
        view(view_angle);
        
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
end % function