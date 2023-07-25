function viz_d_components(d)

% show contents of d = trck2dsp(tks);
dark_style = true; 

if nargin==0
    [ffname,ffolder] = uigetfile;
    load(fullfile(ffolder,ffname),'d');
end

if length(d)>2; warning('Only taking first and last timepoint'); end


% show components

figure('Position', [1.4277e+03 317 1.0847e+03 988.6667]);
for ii=1:3
    subplot(2,2,ii);
    scatter3(d(end).r(:,1),d(end).r(:,2),d(end).r(:,3),...
        [], d(end).dr(:,ii),'.');
    title(sprintf('d(end).dr(:,%d)',ii));
    colorbar('Location','southoutside')
    
    if dark_style
   set(gca,'Color',0.3*ones(3,1));
   colormap_ice;
   
 %  caxis(max(abs(caxis))*[-1,1]);
   alpha = 0.01;
   caxis(max(abs([quantile(d(end).dr(:,ii),alpha), quantile(d(end).dr(:,ii),1-alpha), ]))*[-1,1]);
    
    end
    
    view([0,0]);

end
subplot(2,2,4);
scatter3(d(end).r(:,1),d(end).r(:,2),d(end).r(:,3),...
    [], sqrt(d(end).dr(:,1).^2+d(end).dr(:,2).^2+d(end).dr(:,3).^2),'.');
title(sprintf('|d(end).dr|',ii));
colorbar('Location','southoutside')

 


% show histogram
if false
    figure('Position', [1.4277e+03 317 1.0847e+03 988.6667]);
    for ii=1:3
        subplot(2,2,ii);
        histogram( d(end).dr(:,ii));
        title(sprintf('d(end).dr(:,%d)',ii));
        
    end
    subplot(2,2,4);
    histogram( sqrt(d(end).dr(:,1).^2+d(end).dr(:,2).^2+d(end).dr(:,3).^2));
    title(sprintf('|d(end).dr|',ii));
    
    
end

end

