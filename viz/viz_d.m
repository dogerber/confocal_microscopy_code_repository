function viz_d(d)

% show contents of d = trck2dsp(tks);


if nargin==0
    [ffname,ffolder] = uigetfile;
    load(fullfile(ffolder,ffname),'d');
end

if length(d) <3 % not a timeseries
    
    
    if length(d(1).r(1,:)) >2 %3d
        
        % figure;
        subplot(2,2,1)
        scatter3(d(1).r(:,1),d(1).r(:,2),d(1).r(:,3),[],d(1).r(:,3),'.');
        xlabel('x');ylabel('y');
        title('d(1).r')
        
        subplot(2,2,2)
        scatter3(d(2).r(:,1),d(2).r(:,2),d(2).r(:,3),[],d(2).r(:,3),'.');
        xlabel('x');ylabel('y');
        title('d(2).r')
        
        subplot(2,2,3)
        quiver3(d(1).r(:,1),d(1).r(:,2),d(1).r(:,3),d(2).dr(:,1),d(2).dr(:,2),d(2).dr(:,3),0);
        xlabel('x');ylabel('y');
        title('d(1).r with d(2).dr')
        
        subplot(2,2,4)
        quiver3(d(1).r(:,1),d(1).r(:,2),d(1).r(:,3),d(2).dr(:,1),d(2).dr(:,2),d(2).dr(:,3),0); hold on
        scatter3(d(1).r(:,1),d(1).r(:,2),d(1).r(:,3),'bo');
        scatter3(d(2).r(:,1),d(2).r(:,2),d(2).r(:,3),'r.');
        subtitle('blue to red');
        
        
        
        % quiver3(d(2).r(:,1),d(2).r(:,2),d(2).r(:,3),d(2).dr(:,1),d(2).dr(:,2),d(2).dr(:,3));
        % xlabel('x');ylabel('y');
        % title('d(2).r with d(2).dr')
        
        if false
            figure
            quiver3(d(1).r(:,1),d(1).r(:,2),d(1).r(:,3),d(2).dr(:,1),d(2).dr(:,2),d(2).dr(:,3));
            xlabel('x');ylabel('y');
            title('d(1).r with d(2).dr AUTO SCALED')
        end
        
        % HISTOGRAM OF DX DY DZ
        tp=2;
        colors = cbrewer('div','Spectral',5);
        if tp>1 % skip first timepoint
            figure;
            for dim_i =[1:3]
                subplot(2,2,dim_i); hold on;
                
                if false % normal histogram on top of each other
                    
                    histogram(d(tp).dr(:,dim_i),...
                        'FaceColor',colors(tp,:),'FaceAlpha',0.4);
                else % 3d bars
                    [N,edges] = histcounts(d(tp).dr(:,dim_i));
                    edges = (edges(1:end-1)+edges(2:end))/2;
                    hold on;
                    plot3(edges,tp*ones(length(N),1),N,...
                        'Color',colors(tp,:));
                    view([45,45]);
                end
            end
            subplot(2,2,1); title('dx, tp as color'); xlabel('dx'); ylabel('t-idx'); zlabel('count'); grid on
            subplot(2,2,2); title('dy, tp as color');xlabel('dy'); ylabel('t-idx'); zlabel('count'); grid on
            subplot(2,2,3); title('dz, tp as color');xlabel('dz'); ylabel('t-idx'); zlabel('count'); grid on
            
            subplot(2,2,4); hold on;
            scatter3(d(tp).dr(:,1),d(tp).dr(:,2),d(tp).dr(:,3),...
                [],colors(tp,:),'.');
            xlabel('dx'); ylabel('dy'); zlabel('dz');
            view([45,45]); grid on
        end
        
        if true
            figure
            quiver3(d(1).r(:,1),d(1).r(:,2),d(1).r(:,3),d(2).dr(:,1),d(2).dr(:,2),d(2).dr(:,3),0); hold on
            scatter3(d(1).r(:,1),d(1).r(:,2),d(1).r(:,3),'bo');
            scatter3(d(2).r(:,1),d(2).r(:,2),d(2).r(:,3),'r.');
            subtitle('blue to red');
            axis equal
        end
        
        % show components
        if true
            figure('Name','Displacement components shown on final position')
            marker_scatterplot = 'x';
            for ii=1:3
                subplot(2,2,ii); hold on;
                scatter3(d(end).r(:,1),d(end).r(:,2),d(end).r(:,3),...
                    [], d(end).dr(:,ii),marker_scatterplot);
                title(sprintf('d(end).dr(:,%d)',ii));
                colorbar('Location','southoutside')
            end
            subplot(2,2,4); hold on; 
            scatter3(d(end).r(:,1),d(end).r(:,2),d(end).r(:,3),...
                [], sqrt(d(end).dr(:,1).^2+d(end).dr(:,2).^2+d(end).dr(:,3).^2),marker_scatterplot);
            title(sprintf('|d(end).dr|',ii));
            colorbar('Location','southoutside')
            
            figure('Name','Displacement components shown on reference position')
            for ii=1:3
                subplot(2,2,ii);
                scatter3(d(1).r(:,1),d(1).r(:,2),d(1).r(:,3),...
                    [], d(end).dr(:,ii),'.');
                title(sprintf('d(end).dr(:,%d)',ii));
                colorbar('Location','southoutside')
            end
            subplot(2,2,4);
            scatter3(d(1).r(:,1),d(1).r(:,2),d(1).r(:,3),...
                [], sqrt(d(end).dr(:,1).^2+d(end).dr(:,2).^2+d(end).dr(:,3).^2),'.');
            title(sprintf('|d(end).dr|',ii));
            colorbar('Location','southoutside')
            
        else % 2d
            
            
            figure;
            subplot(2,2,1)
            scatter(d(1).r(:,1),d(1).r(:,2));
            xlabel('x');ylabel('y');
            title('d(1).r')
            
            subplot(2,2,2)
            scatter(d(2).r(:,1),d(2).r(:,2));
            xlabel('x');ylabel('y');
            title('d(2).r')
            
            subplot(2,2,3)
            quiver(d(1).r(:,1),d(1).r(:,2),d(2).dr(:,1),d(2).dr(:,2));
            xlabel('x');ylabel('y');
            title('d(1).r with d(2).dr')
            
            subplot(2,2,4)
            quiver(d(2).r(:,1),d(2).r(:,2),d(2).dr(:,1),d(2).dr(:,2),0);
            xlabel('x');ylabel('y');
            title('d(2).r with d(2).dr')
            
            
        end
        
        
        if true
            figure('Name','3D displacement relative to first tp')
            subplot(2,2,1)
            scatter3(d(2).r(:,1),d(2).r(:,3),d(2).dr(:,2),[],d(2).dr(:,2))
            xlabel('x'); ylabel('z'); zlabel('dy');
            subplot(2,2,2)
            scatter3(d(2).r(:,2),d(2).r(:,3),d(2).dr(:,1),[],d(2).dr(:,1))
            xlabel('y'); ylabel('z'); zlabel('dx');
            
            subplot(2,2,3)
            scatter3(d(2).r(:,1),d(2).r(:,2),d(2).dr(:,3),[],d(2).dr(:,3))
            xlabel('x'); ylabel('y'); zlabel('dz');
        end
        
        
        %% Interpolate displacement field onto regular grid
        if false
            max_nr_pts = inf; % ODD please, type memory into console to see max array size
            [X,Y,dx] = create_grid(d,max_nr_pts);
            param.dx = dx;
            i=2;
            
            
            for j=1:3
                F=scatteredInterpolant(d(1).r(:,1), d(1).r(:,2), d(2).dr(:,j)); % From tp = 1, with displacement of tp = 2
                
                d(2).u(j).u = F(X,Y);
            end
            
            
            figure('Name','Interpoalted displacement components');
            for j=1:3
                subplot(2,2,j)
                imagesc('XData',X(:),'YData',Y(:),'CData',[d(i).u(j).u]);axis image; hold on
                % scatter(d(i).r(:,1), d(i).r(:,2),[], d(i).dr(:,j),'x'); % raw
                % input to double check
            end
        end
        
    else % length(d) >2
        fprintf('d as timeseries detected, using viz_d_timeseries instead... \n');
        viz_d_timeseries(d);
        
    end
    
    
    
end