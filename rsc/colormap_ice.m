
ctype = 'div';
cname = 'RdBu';
ncol = 100;

colormap_ice_real = cbrewer(ctype, cname, ncol);

colormap_ice_real = flip(colormap_ice_real); % so blue is negative (freezing pushing) red is up (melting up)

colormap(gca,colormap_ice_real);