function viz_d_timeseries(d)
% show contents of d = trck2dsp(tks) with d a (1xTIMEPOINTS) structure



if length(d(1).r(1,:)) >2 %3d
    figure
    subplot(3,1,1);
    quiver3(d(1).r(:,1),d(1).r(:,2),d(1).r(:,3),...
        d(end).dr(:,1),d(end).dr(:,2),d(end).dr(:,3),0);
    axis equal
    title('start to end');
    subtitle('not autoscaled');
    
    subplot(3,1,2);
    quiver3(d(1).r(:,1),d(1).r(:,2),d(1).r(:,3),...
        d(end).dr(:,1),d(end).dr(:,2),d(end).dr(:,3),1);
    subtitle('autoscaled');
    
    subplot(3,1,3); hold on;
    plot3(d(1).r(:,1),d(1).r(:,2),d(1).r(:,3),'r.','DisplayName','T0');
    plot3(d(end).r(:,1),d(end).r(:,2),d(end).r(:,3),'bo','DisplayName','T end');
    title('Start and end positions'); grid on; legend
end


fprintf('Continue with detailed analysis? Press a button.\n'); pause

f1 = figure;
f_histograms = figure;
hold on

%random_colors = rand(length(d(1).r(:,1)),3);
if length(d) >2
colors = cbrewer('div','Spectral',length(d));
else
    colors = cbrewer('qual','Set1',length(d));
end

for tp = 1:length(d)-1
    
    diff_step = d(tp+1).r(:,:)-d(tp).r(:,:);
    
    % if length(d(tp).r(1,:)) >2 %3d
    % quiver3(d(tp).r(:,1),d(tp).r(:,2),d(tp).r(:,3),...
    %     d(tp+1).dr(:,1),d(tp+1).dr(:,2),d(tp+1).dr(:,3));
    % else % 2d
    % quiver(d(tp).r(:,1),d(tp).r(:,2),...
    %     d(tp+1).dr(:,1),d(tp+1).dr(:,2));
    % end
    
    if length(d(tp).r(1,:)) >2 %3d
        figure(f1);
        subplot(3,1,1); hold on
        quiver3(d(tp).r(:,1),d(tp).r(:,2),d(tp).r(:,3),...
            diff_step(:,1),diff_step(:,2),diff_step(:,3),...
            0,... %disable scaling
            'Color',colors(tp,:));
        view([45,45]);
        title('quiver plot of displacements, color for timepoint')
        
        subplot(3,1,2); hold on
        scatter3(d(tp).r(:,1),d(tp).r(:,2),d(tp).r(:,3),[],...
            ones(length(d(tp).r(:,1)),1)*tp,'.')
        view(45,45);
        title('colored by timepoint');
        colorbar
        
        subplot(3,1,3); hold on
        scatter3(d(tp).r(:,1),d(tp).r(:,2),d(tp).r(:,3),[],...
            1:length(d(tp).r(:,1)),'.')
        view(45,45);title('colored by particlenumber');
        
        colormap(colors);
        % styling
        xlabel('x');ylabel('y');
        
        

        % histogram of dx,dy,dz over time
        if tp>1% (length(d)-1) % skip first timepoint
            figure(f_histograms);
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
        axis equal
        
        
        end 

        
    else % 2d
        
%         quiver(d(tp).r(:,1),d(tp).r(:,2),...
%             diff_step(:,1),diff_step(:,2),0,...
%             'Color','b');
        
                subplot(3,1,1); hold on
        quiver(d(tp).r(:,1),d(tp).r(:,2),...
            diff_step(:,1),diff_step(:,2),...
            0,...
            'Color',colors(tp,:));

        
        subplot(3,1,2); hold on
        scatter(d(tp).r(:,1),d(tp).r(:,2),[],...
            ones(length(d(tp).r(:,1)),1)*tp,'.')
        title('colored by timepoint');
        colorbar
        
        subplot(3,1,3); hold on
        scatter(d(tp).r(:,1),d(tp).r(:,2),[],...
            1:length(d(tp).r(:,1)),'.')
       title('colored by particlenumber');
        
        colormap(colors);
        if false
            
            hold on
            scatter(d(tp).r(:,1),d(tp).r(:,2),[],[1:length(d(tp).r(:,1))],'o');
            
        end
    end
    
end % fro tp






%% more extensive plots for debugging
if false
if input_YN('Extensive plots is on. Type [Y/N] to continue. ')
    
    for tp = 1:length(d)
        figure
        subplot(2,1,1)
        plot3(d(tp).r(:,1),d(tp).r(:,2),d(tp).r(:,3),'.');
        title(sprintf('d(%d).r',tp));
        
        subplot(2,1,2)
        plot3(d(1).r(:,1),d(1).r(:,2),d(1).r(:,3),'b.'); hold on
        quiver3(d(1).r(:,1),d(1).r(:,2),d(1).r(:,3),...
            d(tp).dr(:,1),d(tp).dr(:,2),d(tp).dr(:,3),0);
        title(sprintf('d(%d).dr from d(1).r',tp));
        plot3(d(tp).r(:,1),d(tp).r(:,2),d(tp).r(:,3),'r.');
    end
    
    
end



if input_YN('Save d_z images here? [Y/N]')
    for tp = 1:length(d)
        figure('Visible','off');
        scatter(d(1).r(:,1),d(1).r(:,2),[],d(tp).dr(:,3),'.');
        c = colorbar;
        c.Label.String  = 'd_z';
        [~,f_name] = fileparts(d(tp).path_stressed);
        save_name = strcat('dz_',f_name,'.png');
                title(save_name);
        export_fig(gcf,save_name,'-dpng','-r300');
    end
    
    
end

end

end