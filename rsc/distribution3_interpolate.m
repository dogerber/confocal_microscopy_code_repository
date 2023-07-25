function distribution3_interpolate(pks)
cd_start =cd;
do_pts = false;

if nargin == 0
    [file_name,path_file_OR_mat,~]=uigetfile('*.mat');
    cd(path_file_OR_mat);
    load(file_name)
end

x = pks(:,1); y = pks(:,2); z = pks(:,3);

% Prepare Meshgrid
stepsize = 10;
[xq,yq] = meshgrid(min(x):stepsize:max(x), min(y):stepsize:max(y));
range_x = max(x)-min(x);
range_y = max(y)-min(y);

vq = griddata(x,y,z,xq,yq);

s= mesh(xq,yq,vq);
s.FaceAlpha = 0.9;
s.FaceColor = 'interp';
s.EdgeColor = 'interp';
s.FaceLighting = 'none';

% correct aspect ratio
pbaspect([range_x,range_y,range_x]);

% caxis([15 45]); % sets the borders of the c-axis color

if do_pts
    hold on
    plot3(x,y,z,'b.')
    hold off;
end

%shg

% Contour Plot
% figure;
% levels_contour = 20;
% contour(xq,yq,vq, levels_contour);

cd(cd_start);

end