function [sigma, roi_ref] = tfm_3d_DG_v2(PATH_DISPLACEMENTS, param,roi_ref)
%%
%
%  Function to run 3d TFM analysis WITHOUT LN2 regularization
%  for displacements located on a plane 
%
%
% INPUT
%    - xyz positions of beads in stressfree and deformed state given as
%    path to .mat file containing var 'pks' or directly as a matrix
%    - param various control or '' to take default
%    - MM, ROI: from other displacmeents in same series to speed things up
%
%
% OUTPUT
%    - Images of stresses (saved in same folder)
%    - Matrix of stresses (where, which direction)

% Dominic Gerber, 2020; original by Dr. Kathryn A. Rosowski and Dr. Robert Style
% Branched from tfm_3d_DG_LN2.m 2.2.2022

% Notes:
% use [xq,szq] = parallel_slicer_v2([s.x,s.y,s.sz]); figure; plot(xq,szq)
% to check result

% Options
if ~exist('param','var');      param = struct(); end
if isempty(param); param = struct(); end
if ~isfield(param,'do_debug'); param.do_debug =  true; end
if ~isfield(param,'do_talk');  param.do_talk = true; end
if ~isfield(param,'do_show');  param.do_show = true; end
if ~isfield(param,'do_save');  param.do_save = true; end% saving plots
if ~isfield(param,'remove_bad_vectors'); param.remove_bad_vectors = false; end
if ~isfield(param,'fourier_filter'); param.fourier_filter = 2; end % fourier low-pass filter, removing noise
% this is equivalent to the min. feature size expected

% Physical sample substrate Parameters
if ~isfield(param,'Ey'); param.Ey = 265e3; end%150*10^3; %[Pa]
if ~isfield(param,'max_nr_pts'); param.max_nr_pts = inf; end% sqrt(nr_grid_points), ODD please, type memory into console to see max array size
%if ~isfield(param,'substrate_thickness'); param.substrate_thickness = 1e-6; end% [meters]
if ~isfield(param,'nu'); param.nu = 0.4999; end% not exactly 0.5 to avoid problems



% ---
if isunix
    param.do_show = false; param.do_debug = false;
end

if param.do_show
    figure_visible = 'on';
else
    figure_visible = 'off';
end


if param.do_debug; tic; disp(param); end

%% Load and prepare displacements

% if nargin >0
%     if isstruct(PATH_DISPLACEMENTS)
%         param = PATH_DISPLACEMENTS;
%         clear PATH_DISPLACEMENTS
%     end
% end

if ~exist('PATH_DISPLACEMENTS','var')
    [ffname,ffolder] = uigetfile('Choose the folder in which the slices are','*.mat*');
    PATH_DISPLACEMENTS = fullfile(ffolder, ffname);
    param.fullfilename_in = PATH_DISPLACEMENTS;
else
    try
    [ffolder,ffname] = fileparts(PATH_DISPLACEMENTS);
    catch
        
    end
end

if isstruct(PATH_DISPLACEMENTS) % input was a d
    d = PATH_DISPLACEMENTS;
else
    load(PATH_DISPLACEMENTS,'d'); % dispalcements structure e.g. from confocal_to_displacments.m
end


% Convert everything to meters
if true
    for tp=1:length(d)
        for dim = 1:3
            d(tp).r(:,dim) = d(tp).r(:,dim)/10^6;
            d(tp).dr(:,dim) = d(tp).dr(:,dim)/10^6;
        end
        % d(tp).r(:,3) = d(tp).r(:,3) + 100e-6;
    end
end


% look up substrate thickness from measruement
param.substrate_thickness = mean(d(1).r(:,3)); % convert to [meters]

if param.do_talk
    fprintf('Using the substrate thickness from displacement data d(): %d m \n',param.substrate_thickness);
end





%% Remove Bad vectors
if param.remove_bad_vectors
    %bad_vector_area = 30e-6; % area where a vector is comapred to it's neighbours,
    % should be big enough to include many neighbours, but smaller than the
    % area where changes of the vetctors are expected.
    bad_vector_area = 3* sqrt((max(d(1).r(:,1))-min(d(1).r(:,1)))*(max(d(1).r(:,2))-min(d(1).r(:,2)))/length(d(1).r(:,1)));
    % disp(bad_vector_area);
    
    param_dbv = struct();
    param_dbv.factor = 2; % by what factor must a vector diverge to be considered bad?
    
    if param.do_debug
        param_dbv.do_show = true; % shows result
        [d, ~] = delete_bad_vectors_3d(d, bad_vector_area, param_dbv); % (d, win, method, factor, interp_mode, graphics)
        % fprintf('in black, vectors that were deleted in red, and vectors that didnt have enough neighbros (3) in green.\n');
    else
        param_dbv.do_show = false;
        [d, ~] = delete_bad_vectors_3d(d, bad_vector_area,  param_dbv); % does not open a figure which shows the result
    end
end



%% Interpolate displacement field onto regular grid
[X,Y,dx] = create_grid(d,param.max_nr_pts);
grid_size = length(X(:,1));
param.dx = dx;
i=2;


for j=1:3
    F=scatteredInterpolant(d(1).r(:,1), d(1).r(:,2), d(2).dr(:,j)); % From tp = 1, with displacement of tp = 2
    d(2).u(j).u = F(X,Y);
    
end



% zero displacement interpolated
ind2=find(isnan(d(i).u(1).u));
d(i).fov=ones(size(d(i).u(1).u));
d(i).fov(ind2)=0;

%make field of view cumulative
d(i).fov=d(i).fov.*d(i).fov;

if param.do_debug
    figure('Name','X Y and Z magnitude of displacement');
    imagesc([d(i).u(1).u,d(i).u(2).u,d(i).u(3).u]);axis image;
    
    figure
    for j=1:3
        subplot(2,2,j)
        imagesc('XData',X(:),'YData',Y(:),'CData',[d(i).u(j).u]);axis image; hold on
        scatter(d(i).r(:,1), d(i).r(:,2),[], d(i).dr(:,j),'x');
    end
end



d(i).X = X;
d(i).Y = Y;
d(i).dx = dx;

% show interpolated data
if param.do_debug
    %3d quiver
    %%%%% this is not to be trusted, quiver3 does something weird! check
    %%%%% with histogram or scatter3 instead!
    figure
    scatter3(d(2).r(:,1), d(2).r(:,2), d(2).dr(:,3),'r.');
    hold on
    scatter3(X(:)+d(2).u(1).u(:),Y(:)+d(2).u(2).u(:),d(2).u(3).u(:),'bo');
    subtitle('red raw, blue interpolated');
    
        figure; hold on
    scatter(d(2).r(:,1), d(2).r(:,2), [],d(2).dr(:,3),'.');
    figure; hold on
    scatter(X(:)+d(2).u(1).u(:),Y(:)+d(2).u(2).u(:),[],d(2).u(3).u(:),'.');
end



%% calculate M matrix
if false % old way of doing it, might run out of RAM
    warning('Using old way of TFM analysis.');
    if ~exist('MM','var')
        %%MM = tfm_ln2_regularization(grid_size,param);
        MM = tfm_ln2_regularization(grid_size,param);
    end
    
    % Replace Nan's with zeros
    for i=1:3
        d(2).u(i).u(isnan(d(2).u(i).u)) = 0;
    end
    
    Ux=fftshift(fft2(d(2).u(1).u)); % FT of ux
    Uy=fftshift(fft2(d(2).u(2).u)); % FT of uy
    Uz=fftshift(fft2(d(2).u(3).u)); % FT of uz
    
    
    %- Apply Window
    if true
        % Prepare window
        [szr,szc]=size(d(2).u(1).u); % IS THIS RIGHT?
        w_c = tukeywin(szc,0.25)';
        w_r = tukeywin(szr,0.25)';
        wn=w_c'*w_r;
        
        % apply
        Ux2=fft2(ifft2(Ux).*wn);
        Uy2=fft2(ifft2(Uy).*wn);
        Uz2=fft2(ifft2(Uz).*wn);
    end
    
    %- Calculate stresses
    if true % old, MM as one huge (3*grid_size^2)^2 matrix, which can easily run out of RAM
        S_diag = MM\[Ux2(:);Uy2(:);Uz2(:)];
        [sx,sy,sz] =local_reshape_S(S_diag,grid_size);
    else % try to omit diagonal matrix with 3 grid_size^2 matrices
        warning('Experimental!');
        sx = MM(1:grid_size^2,1:grid_size^2)\Ux2(:);
        sy = MM(1*grid_size^2+1:2*grid_size^2,1*grid_size^2+1:2*grid_size^2)\Uy2(:);
        sz = MM(2*grid_size^2+1:3*grid_size^2,2*grid_size^2+1:3*grid_size^2)\Uz2(:);
        
        sx = reshape(sx,grid_size,grid_size);
        sy = reshape(sy,grid_size,grid_size);
        sz = reshape(sz,grid_size,grid_size);
    end
    
    
else % do it after disp2stress.m (From Eric 2009)
    
    % - calculate Q matrix
    thick=param.substrate_thickness;
    EM=param.Ey ;
    nu=param.nu;
    dx=param.dx;
    nx = grid_size;
    
    Q = calcQ(thick,thick,EM,nu,nx,dx,3);
    
    
    u = d(2).u;
    
    % Apply turkey window on displacement data edge
    if true
        % Prepare window
        [szr,szc]=size(u(1).u); % IS THIS RIGHT?
        w_c = tukeywin(szc,0.25)';
        w_r = tukeywin(szr,0.25)';
        wn=w_c'*w_r;
        
        % apply
        for i=1:3
            u(i).u = u(i).u.*wn;
        end
    end
    
    % Wiener filter on displacement data?
    if false
        for i=1:3
            u(i).u =wiener2(u(i).u);
        end
    end
    
    
    
    % Take the FFT of each component of the displacement field
    for i=1:3
        u(i).U =fftshift(fft2(u(i).u));
    end
    
    
    % Gaussian low-pass filter
    if param.fourier_filter ~= 0
        if param.do_debug; fprintf('Fourier filtering with min featuer size %d \n', param.fourier_filter); end
        min_feature_size =  param.fourier_filter;
        
        nr2 = length(u(1).u(:,1));
        qmax=nr2/(pi*min_feature_size);
        
        % Get distance from of a grid point from the centre of the array
        y=repmat((1:nr2)'-nr2/2,1,nr2);
        x=y';
        q=sqrt(x.^2+y.^2);
        
        % Make the filter
        qmsk=exp(-(q./qmax).^2);
        qmsk=(qmsk); % or  ifftshift ?
        
        % apply
        for i=1:3
            u(i).U = u(i).U .* qmsk;
        end
    end
    
    
    % Calculate the FFT of each component of the stress field
    for i=1:3
        s(i).S = 0;
        for j=1:3
            s(i).S = s(i).S + (Q(i,j).Q).*(u(j).U);
        end
    end
    
    % IFFT the stresses
    for i=1:3
        s(i).s =ifft2(ifftshift(s(i).S),'symmetric');
    end
    
    
    sx = s(1).s;
    sy = s(2).s;
    sz = s(3).s;
    
end


    f_noLN2 = figure('Name','Without LN2','Visible',figure_visible);
    subplot(2,2,1)
    imagesc(real(sx));
    title('sx');
    colorbar
    
    subplot(2,2,2)
    imagesc(real(sy));
    title('sy');
    colorbar
    
    subplot(2,2,3)
    imagesc(real(sz));
    title('sz');
    colorbar
    drawnow





%% subtract mean value in water area

if ~exist('roi_ref','var')
    fprintf('Select where stresses are 0 as reference.\n');
    [sz_water,roi_ref] = cut_selection([X(:),Y(:),real(sz(:))]);
else
    sz_water = cut_selection([X(:),Y(:),real(sz(:))],roi_ref);
end

mn_sz_water = mean(sz_water(:,3));

sz = sz-mn_sz_water;

f_sz = figure('Visible',figure_visible);
imagesc('XData',X(:),'YData',Y(:),'CData',real(sz)); colorbar
title('sz after water subtracted');

if param.do_debug
    figure;
    surf(sz);
    
    figure;
    plot(sz(10,:),'-');
end




%% export

sigma = struct();
sigma.d = d; % copy d() over

% results in matrix form
sigma.Sx = sx;
sigma.Sy = sy;
sigma.Sz = sz;
sigma.X = X;
sigma.Y = Y;

% results in vector form
sigma.x = X(:);
sigma.y = Y(:);
sigma.sx = real(sx(:));
sigma.sz = real(sy(:));
sigma.sz = real(sz(:));

% note down infos of tfm
sigma.tfm.param = param;
sigma.tfm.date_done = datetime('now');
if ~isstruct(PATH_DISPLACEMENTS)
    sigma.tfm.PATH_DISPLACEMENTS = PATH_DISPLACEMENTS;
else
    sigma.tfm.PATH_DISPLACEMENTS = d(2).path_stressed;
end





if param.do_save
    if ~exist('ffname','var') % timeseries mode has no file input
        ffname = d(2).path_stressed;
        ffolder = pwd; % current folder
    end
    
    core_name = extractAfter(ffname,'rel_');
    export_name = strcat('sigma_', core_name,'.mat');
    export_full_path = fullfile(ffolder,export_name);
    save(export_full_path,'sigma');
    
    if false
    savefig(f_noLN2,fullfile(ffolder,...
        strcat('raw_sigma_overview',core_name,'.fig')));
    
    savefig(f_sz,fullfile(ffolder,...
        strcat('sz_', core_name,'.fig')));
    end
    
    export_fig(f_sz,...
        fullfile(ffolder,...
        strcat('sz_', core_name,'.png')),...
        '-r300','-dpng');
    
    if param.do_talk; fprintf('Succesfully saved to \n \t %s.\n',export_full_path); end
end

if param.do_debug; toc; end

end % function




%% ------- Sub-Functions -------


function [sx,sy,sz] =local_reshape_S(S,grid_size)
i=1;
sx= reshape(S(1:1*grid_size^2,i),grid_size,grid_size);
sy= reshape(S(1*grid_size^2+1:2*grid_size^2,i),grid_size,grid_size);
sz=reshape(S(2*grid_size^2+1:3*grid_size^2,i),grid_size,grid_size);
end


