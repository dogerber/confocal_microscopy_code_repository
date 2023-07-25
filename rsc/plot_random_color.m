function color_out = plot_random_color(num_colors)
% color_out = plot_random_color(num_colors)
% returns "random" number of rgb color tripletes
%

if nargin<1
    num_colors = 1;
end

colors = cbrewer('qual','Set1',max(100,num_colors));
random_number = randperm(length(colors));

color_out = colors(random_number(1:num_colors),:);


end