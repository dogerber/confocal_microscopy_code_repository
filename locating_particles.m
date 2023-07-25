function [pks,pks_px,param] = locating_rob2023(FILEPATH_IN,image_stack,param)
% [pks,pks_px,param] = locating_rob2023(FILEPATH_IN)
% Locating spherical blobs in 3d confocal images
% can be batched with files_feeder().m
%
% ** Version history **
param.code_version = '1.0.4';
% 1.0.4 - added proper param integration
% 1.0.3 - added outputing the maximum intensty to better identify particles
%       - added timing capability
% 1.0.2 - added ability to read coordinates from the filename (enables
%         stitching of multiple positions)
%       - added ability to detect if target file already exists
%       - use tiffreadVolume to speed up import
% 1.0.1 - initial commit
%
% ------------------------------------- Rob Style, February 2023

% Options
if ~exist('param','var'); param = struct(); end
if ~isfield(param,'save_folder'); param.save_folder = pwd; end%fullfile('C:\wf_local\SHaohua check localization');
if ~isfield(param,'save_in_relative_folder_structure'); param.save_in_relative_folder_structure = false;end

if ~isfield(param,'do_show'); param.do_show = false; end
if ~isfield(param,'do_debug'); param.do_debug =  false; end
if ~isfield(param,'do_talk');  param.do_talk = true; end
param.do_time = true; 

% Locating parameters
if ~isfield(param,'threshold'); param.threshold = 400; end% Intensity param.threshold for local maxima
if ~isfield(param,'subvolume_size'); param.subvolume_size = [5,5,7]; end% size xyz in voxel around the particle point spread function used for sub-pixel accuracy fitting
%param.min_separation = [11,11,11]; % minimum distance in [x,y,z] between particles, must be odd!

% Microscope parameters
if ~isfield(param,'muperpixel'); param.muperpixel =  0.329131; end%  0.329131 [mum/pxel] given for 20x air on conny
if ~isfield(param,'z_step'); param.z_step = 1; end% [mum]
if ~isfield(param,'width_fov_mum'); param.width_fov_mum = param.muperpixel * 2048; end% width of total FOV in [mum]

%% check input
if any(~mod(param.subvolume_size,2))
    error('param.subvolume_size must contain all odd numbers!')
end

%% Load image data and prepare stuff
if param.do_time; t_start = tic; end

if nargin<1
    FILEPATH_IN = uigetfile_to_fullpath('*.tif*');
end

% get file name
[filepath,fname,ext] = fileparts(FILEPATH_IN);

% get file infos
info = imfinfo(FILEPATH_IN);
num_slices = numel(info);



% get xyz position of objective, if files were exported in 3i SlideBook
% with [%N_T%T_x%x_y%y_z%z_C%c]
x0 = str2double(extractBetween(string(fname),'_x','_y'));
y0 = str2double(extractBetween(string(fname),'_y','_z'));
z0 = str2double(extractBetween(string(fname),'_z','_C'));

if any([isempty(x0), isempty(y0), isempty(z0)])
    fprintf('[!] Recognition of real space parameters failed. This only works if files were exported in 3i SlideBook with [%N_T%T_x%x_y%y_z%z_C%c].\n Showing pixel and slice numbers instead.\n');
    x0 = 0; y0 = 0; z0 = 0;
end

% determine save path
save_name = strcat(fname,'_pos.mat');
if param.save_in_relative_folder_structure
    rel_folder_path = extractAfter(filepath,':\');
    target_folder_path = fullfile(param.save_folder,rel_folder_path);
    if ~exist(target_folder_path,'dir') mkdir(target_folder_path); end
    ff_save = fullfile(target_folder_path, save_name);
else
    ff_save = fullfile(param.save_folder,save_name);
end

% skip if file already exists
if exist(ff_save,'file')
    fprintf('File %s does already exist, skipping locating_rob2023\n',string(ff_save));
    return;
end

%% import image
if nargin<2
try
    image_stack = tiffreadVolume(FILEPATH_IN);
catch
    fprintf('It seems like you dont have the matlab function tiffreadVolume, which was introduced in 2020. The old method i use instead makes the code 15 times slower');
    image_stack = zeros(info(1).Height, info(1).Width, num_slices);
    for k = 1:num_slices
        image_stack(:,:,k) = imread(FILEPATH_IN, k, 'Info', info);
    end
    
end
end

% adaptive thresholding (WIP)
if false
    T = adaptthresh(image_stack,0.4);
    image_stack_bin = imbinarize(image_stack,T);
    J = image_stack_bin(1:100,1:100,:);
    figure
    subplot(2,1,1);
    slice(double(J),size(J,2)/2,size(J,1)/2,size(J,3)/2)
    colormap gray
    shading interp
    subplot(2,1,2);
    slice_midplane(image_stack(1:100,1:100,:))
end

%% Find local maxima

structuring_element = true(3,3,3);
structuring_element(2,2,2) = false; % hole in the middle
blurred_stack = imgaussfilt3(image_stack); % Apply Gaussian blur to remove noise
dilated_stack = imdilate(blurred_stack, structuring_element);

above_threshold = (image_stack > param.threshold);
local_maxima_mask = (blurred_stack > dilated_stack) & above_threshold;
% the first part makes only the middle of the blob to be true, the second
% only passes things which are above the threshold
if sum(local_maxima_mask(:))<1; error('No local maxima found, adjust threshold!'); end
%  figure;
% slice_nr = 100;
% tmp = (blurred_stack > dilated_stack);
% subplot(2,2,1); imagesc(image_stack(:,:,slice_nr)); title('image_stack');
% subplot(2,2,2); imagesc(above_threshold(:,:,slice_nr)); title('above_threshold');
% subplot(2,2,3); imagesc(tmp(:,:,slice_nr)); title('(blurred_stack > dilated_stack)');
% subplot(2,2,4); imagesc(local_maxima_mask(:,:,slice_nr)); title('local_maxima_mask');

if param.do_debug
    maxima_properties = regionprops3(local_maxima_mask, image_stack,'Centroid','MaxIntensity');
else
    maxima_properties = regionprops3(local_maxima_mask, 'Centroid');
end

maxima_points = maxima_properties.Centroid;

% show estimates of locations with threshold value
if param.do_debug
    figure;
    subplot(2,1,1);
    scatter3(maxima_properties.Centroid(:,1),maxima_properties.Centroid(:,2),maxima_properties.Centroid(:,3),...
        [],maxima_properties.MaxIntensity,'.');
    colorbar;
    caxis([quantile(maxima_properties.MaxIntensity,0.05),quantile(maxima_properties.MaxIntensity,0.95)]);
    title('MaxIntensity n the frist step');
    subplot(2,1,2);
    histogram(maxima_properties.MaxIntensity);
    xlabel('Intensity'); ylabel('count');
end

% figure;
% idx = maxima_properties.MaxIntensity>800;
% scatter3(maxima_properties.Centroid(idx,1),maxima_properties.Centroid(idx,2),maxima_properties.Centroid(idx,3),...
%     [],maxima_properties.MaxIntensity(idx),'.');

if and(false,param.do_debug)
    % Plot local maxima on top of each image (optional)
    figure;
    if true
        for i = 1:num_slices
            imagesc(image_stack(:,:,i));
            hold on;
            plot(maxima_points(maxima_points(:,3) == i,1), maxima_points(maxima_points(:,3) == i,2), 'rx');
            hold off;
            colorbar;
            caxis([0, 800]);
            pause;
        end
    else
        volshow(image_stack)
    end
end

if param.do_debug; fprintf('Found %d local maxima in the first step.\n',length(maxima_points)); end

%% remove particles close to the edges
dist_to_edge = ceil(param.subvolume_size/2);

[image_height, image_width, image_depth] = size(image_stack);
filtered_points = maxima_points(maxima_points(:,1) >= dist_to_edge(1) & maxima_points(:,1) <= image_width-dist_to_edge(1)-1 & ... % x
    maxima_points(:,2) >= dist_to_edge(2) & maxima_points(:,2) <= image_height-dist_to_edge(2)-1 & ... % y
    maxima_points(:,3) >= dist_to_edge(3) & maxima_points(:,3) <= image_depth-dist_to_edge(3)-1, :); % z

%% Subpixel fitting
subvolume_halfsizes = floor(param.subvolume_size/2);

subpixel_points = zeros(size(filtered_points,1),size(filtered_points,2)+1);
for i = 1:length(filtered_points(:,1))
    % Extract a subvolume around each local maximum (needs to have
    % odd dimensions)
    subvolume_at_estimated_z = double(blurred_stack(filtered_points(i,2)-subvolume_halfsizes(1):filtered_points(i,2)+subvolume_halfsizes(1), ...
        filtered_points(i,1)-subvolume_halfsizes(2):filtered_points(i,1)+subvolume_halfsizes(2),...
        filtered_points(i,3))); % should this really be double???
    
    % Perform subpixel fitting for x-y coordinates
    [x_offset, y_offset, ~] = radialcenter(subvolume_at_estimated_z);
    
    % Perform parabola fit for z coordinate
    z_data = (filtered_points(i,3)-3:filtered_points(i,3)+3)';
    a_data = double(blurred_stack(filtered_points(i,2), filtered_points(i,1), filtered_points(i,3)-3:filtered_points(i,3)+3));
    quadratic_coeffs = [ones(7,1), z_data, z_data.^2] \ a_data(:);
    
    % Check if the parabola fit is accurate and discard if not
    z_peak = -quadratic_coeffs(2) / (2*quadratic_coeffs(3));
    if abs(z_peak - filtered_points(i,3)) > 1
        subpixel_points(i,:) = [NaN, NaN, NaN, NaN];
    else
       % integrated_intensity = sum(subvolume_at_estimated_z(:)); % at z plane of max intensity, using blurred picture
        maximum_value = double(blurred_stack(filtered_points(i,2),filtered_points(i,1),filtered_points(i,3)));
        subpixel_points(i,:) = [filtered_points(i,1)-3 + x_offset, ...
            filtered_points(i,2)-3 + y_offset, z_peak, maximum_value];
    end
end

% Remove bad particles
subpixel_points(isnan(subpixel_points(:,1)),:)=[];

if param.do_debug
    fprintf('%d [%2.2g %%] particles were removed, because subpixel accuracy fitting for z-coordinate failes.\n',...
        length(maxima_points)-length(subpixel_points),...
         100*(length(maxima_points)-length(subpixel_points))/ length(maxima_points)); 
end

if param.do_show
    figure
    plot3(subpixel_points(:,1),subpixel_points(:,2),subpixel_points(:,3),'.');
    
    figure; scatter3(subpixel_points(:,1),subpixel_points(:,2),subpixel_points(:,3),[],subpixel_points(:,4)); 
     pd= fitdist(subpixel_points(:,4),'Normal');
 caxis([pd.mu-pd.sigma,pd.mu+pd.sigma]);
end

%% export
pks_px = subpixel_points;

% convert to real space
pks = yflip(subpixel_points); % this makes the orientation the same like in the original image when plotting the data. Plot has (0,0) in the bottom left corner, images in the top left.

% convert to real coordinates
pks(:,1) = x0*(-1)+param.muperpixel*pks(:,1) - param.width_fov_mum/2; % axis of conny are reverted, so *(-1)
pks(:,2) = y0*(-1)+ param.muperpixel*pks(:,2) - param.width_fov_mum/2;
pks(:,3) =  z0*(1)+ param.z_step *(pks(:,3)-num_slices/2);
% s.t. x0,y0,z0 is the middle of the whole stack,


% gather information
locating_infos.param = param;
locating_infos.code_used = 'locating_rob2023.m';
locating_infos.code_used_verison = param.code_version;
locating_infos.pks = pks;
locating_infos.pks_px = pks_px;
locating_infos.orignal_path_of_image = filepath;
locating_infos.date_done = datetime('now');

% actual saving
save(ff_save, 'locating_infos','pks','pks_px');

if param.do_talk
    fprintf('File saved as \n \t %s \n',string(ff_save));
end

if param.do_time; fprintf('Finished in %d seconds\n',toc(t_start)); end

end % function