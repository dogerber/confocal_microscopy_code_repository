function viz_strains_components(res)
%{



%}

%% Options

param.do_show_input = true;
gl_caxis = ''; %[-1,1]*0.5;

do_axis_equal = false; 


%% Look at data
n_dim = size(res.r_start,2);


%% Viz
if n_dim==2 % 2D
    
  
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
        if do_axis_equal; axis equal; end

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
        if do_axis_equal; axis equal; end
        
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
        if do_axis_equal; axis equal; end

    
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
    
  
    
    
    
subplot(2,2,4);
    scatter(res.r_end(:,1),res.r_end(:,2),...
        [],eigen_values_stretch(:,1,1).*eigen_values_stretch(:,2,2),'.');
    colorbar;
    subtitle('\epsilon_{volume}'); xlabel('x'); ylabel('y'); grid on;
    colormap_ice; caxis([0,2])
    if do_axis_equal; axis equal; end
    
    
  
    
else % 3D data =========================================================================
    error('not done for 3d data');
    
    
    
end
end % function