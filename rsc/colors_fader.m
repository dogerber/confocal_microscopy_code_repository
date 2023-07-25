function cmap_out = colors_fader(color_start,color_end,n_colors)
% cmap_out = colors_fader(color_start,color_end,n_colors)
% return colormap that fades between the two input colors (linear, not
% optimized for human perception, see cbrewer for fancier solutions)

cmap_out = nan(n_colors,3);
cmap_out(1,:) = color_start;
cmap_out(end,:) = color_end;

diff_color = (color_end - color_start)/n_colors;


for i=2:n_colors-1
    
   cmap_out(i,:) =  cmap_out(i-1,:) + diff_color;
end




end