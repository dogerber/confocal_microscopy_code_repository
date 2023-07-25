function [good_disp, varargout] = delete_bad_vectors_3d(d, win, param)
%[good_disp, bad_disp] = delete_bad_vectors(d, win, param)
%DELETES ABERRANT VECTORS FROM A SAMPLED VECTOR FIELD, BY COMPARING EACH
%VECTOR TO A VECTOR INTERPOLATED FROM THE SURROUNDING VECTORS.
%
%
% see also Rob's find_outliers.m
%
%REQUIRED INPUTS
%d: is structure of displacements that is output from track2disp.m.
% This is a Nx1 structure, where N is the number of time points. Each time
% point has two fields, the first is a two column array where the first
% colunm is x positions and the second column is y positions, the third z. The second
% field is a two column array where the first column is displacements in x
% and the second column is displacements in y and the third z. Each row corresponds to a
% different particle
%
%win: is a number specifying the size of the region of interest used to
% interpolate the vector field about each sampled location. the region of
% interest is square and this number is the width of the square. this
% number should be large enough
% to include several other beads but not so large that it covers regions
% that vary greatly in displacement magnitude and direction. If it's not big
% enough there will not be enough points to interpolate.  use the optional
% graphical out put, described below, to help you pick a good value of win
%
%OPIONAL INPUTS in structure param
%param.comparing_method is a string that specifies the criterion to distiguish a bad
% vector.  The choices are:
%       'dif' or 'difference':  the magnitude of the vector difference
%       'mag' or 'magnitude':  difference of the magnitude
%       'angle'  the difference between the angles
%       'manual' calls the the manual delete function and bypasses the
%       automatic deletion code
% If you don't specify a method it will do both 'magnitude' and 'difference'.
%
%param.factor is a number specifying how much the vector is allowed to deviate
% from the interpolated vector (default is 2). If you select 'angle' as the
% method, then the factor should be expressed as a multiples of 45 degrees.
% (doesn't have to be an integer!)
%
%param.do_show if true, this code will generate
% a handy figure (default is no figure). The figure shows vectors that were kept
% in black, vectors that were deleted in red, and vectors that didn't have
% enough neighbros (3) in green.  This is a great option
% to turn on to help you determine which method, factor and win to use.
%
%OUTPUTS
% good_disp: all that vectors that were kept.  the structure has the same
%     format as the output of track2disp
% bad_disp: all that vectors that were deleted.  the structure has the same
%     format as the output of track2disp
%
%MODIFICATION HISTORY
% Created by Callen Hyland, Yale University, August 2010
% Incorporated manual delete code written by Guy German and modified to
% allow multiple time points and multiple fields
%
% Dominic Gerber, September 2020: Changed it for 3D capability
%
%KNOWN ISSUES
% 1. Does not yet know how to handle vectors at the edge of the field of view
% 2. small good vectors near large bad vectors sometimes get deleted
% because the large bad vector dominates the interpolation.
%
% DETAILS OF THE MANUAL DELETE CODE:
% 1. quiver plot of input is displayed
% 2. zoom and pan to find desired field
% 3. when you press enter it will give you crosshairs to select vectors-
% click near the base of vectors
% 4. when you press enter again it will display red circles on selected vectors
% 5. Now you can pan and zoom again. press enter to get cross hairs tool
% and enter to select.
% 6. After dispaying red circles, press ESC to end or go to the next frame.


%% Options and Parameters



% Parameters default values
if ~exist('param','var'); param = {}; end
if ~isfield(param,'do_talk'); param.do_talk = false; end % progress to command window
if ~isfield(param,'do_show'); param.do_show = false; end % figure with results displayed in the end


if ~isfield(param,'min_nr_neighbours'); param.min_nr_neighbours = 3; end %% Minimum number of neighbours a particle needs to be compared to. Particles with fewer neighbour are ignored
if ~isfield(param,'keep_not_enough_neighbours'); param.keep_not_enough_neighbours = false; end % What should happen with particles that do not have enough neighbours?
if ~isfield(param,'comparing_method'); param.comparing_method = ''; end % metric to judge the vectors 'angle','dif,'mag' or '' (for mag and angle mixture)
if ~isfield(param,'factor'); param.factor = 2; end %




%% Input handling

if size(d(1).r,2) ~=3
    error('This code is for 3D data ONLY. Please use delete_bad_vectors.m instead.\n');
end

% convert to work with this code
method = param.comparing_method;
factor = param.factor;


%% Main part

time = length(d); % how many timepoints are there
num_beads=length(d(1).r);

to_delete = [];     %holds indices of vectors to delete
ignored = [];       %holds indices of vectors that do not have enough neighbors
win2 =(win/2);

%ignore the first time point because it is always zeros.
for i = 2:time %loop through times
    
    switch lower(method)
        
        case{'manual'}
            %This is the manual delete code modified from G.G.
            cee=[];ree=[]; scl=10;
            quiver(d(i).r(:,1),d(i).r(:,2),d(i).dr(:,1),d(i).dr(:,2),1,'k');
            axis equal; hold on; h1=gcf;
            
            while ishandle(h1)
                zoom on; % use mouse button to zoom in or out
                waitfor(h1,'CurrentCharacter',13)
                [x,y] = ginput; zoom off; plot(x,y,'ro');
                cee = vertcat(cee,x); ree = vertcat(ree,y);
                title('Selected points: ENTER to select more, ESC to accept')
                k = waitforbuttonpress;
                if get(h1,'CurrentCharacter')==27
                    break
                else
                    pause;
                end
            end
            
            ind = [];
            for ii=1:length(cee); %find vectors closest to selection
                distances = (d(i).r(:,1)-cee(ii,1)).^2 + (d(i).r(:,2)-ree(ii,1)).^2;
                [a,ind_current] = min(distances);
                ind = [ind;ind_current];
            end
            
            ind = sort(ind);
            close all
            quiver(d(i).r(:,1),d(i).r(:,2),d(i).dr(:,1),d(i).dr(:,2),1,'k');
            axis equal; hold on
            plot(d(i).r(ind,1),d(i).r(ind,2),'bx')
            title('Closest vectors to selections: Space to continue')
            pause; close all
            to_delete = [to_delete;ind];
            
        otherwise % not manual method
            
            for j = 1:num_beads %loop through all particles
                
                %construct window of specified size around the particle
                pos = d(i).r(j,:);
                x_win = [pos(1)-win2, pos(1)+win2];
                y_win = [pos(2)-win2, pos(2)+win2];
                z_win = [pos(3)-win2, pos(3)+win2];
                
                idwin = find(abs(d(i).r(:,1)-pos(1))<win2 & ... % indices of particles in window
                    abs(d(i).r(:,2) - pos(2))<win2 & ... % COULD BE MADE FASTER WITH RADIUSSEARCH
                    abs(d(i).r(:,3) - pos(3))<win2 & ...
                    (d(i).r(:,1)-pos(1)).^2 + (d(i).r(:,2) - pos(2)).^2 + (d(i).r(:,3) - pos(3)).^2 ~=0); % ignore the point i itself
                
                field = [d(i).r(idwin,1:3),d(i).dr(idwin,1:3)]; % coordinates (and directions) of particles in window
                % field is [x,y,z,dx,dy,dz; *next particle*]
                
                if size(field,1) > param.min_nr_neighbours % check if there are enough points in this window
                    
                    
                    % new for 3d
                    interp3d_dx = scatteredInterpolant(field(:,1:3),field(:,4));
                    interp3d_dy = scatteredInterpolant(field(:,1:3),field(:,5));
                    interp3d_dz = scatteredInterpolant(field(:,1:3),field(:,6));
                    
                    interp_x = interp3d_dx(pos); % evaluate interpolant of dx at pos
                    interp_y = interp3d_dy(pos);
                    interp_z = interp3d_dz(pos);
                    
                    interp_disp = [interp_x, interp_y, interp_z];
                    
                    
                    %apply a criterion to determine if there is a discrepency
                    real_disp = d(i).dr(j,:);
                    angle_between = calc_angle(real_disp,interp_disp); % angle between actual vector and interpolated vector
                    
                    
                    %compare vectors and delete some
                    switch lower(method)
                        case {'mag','magnitude'}
                            if norm(real_disp)> factor*norm(interp_disp)
                                to_delete = [to_delete;j];
                            end
                            
                        case {'dif', 'difference'}
                            if norm(real_disp-interp_disp)> factor*norm(interp_disp)
                                to_delete = [to_delete;j];
                            end
                            
                        case {'angle'}
                            if angle_between>factor*45
                                to_delete = [to_delete;j];
                            end
                        otherwise %if no method is specified do both mag and diff and angle
                            if norm(real_disp)> factor*norm(interp_disp)||...
                                    norm(real_disp-interp_disp)> factor*norm(interp_disp) ||...
                                    angle_between > factor*20
                                to_delete = [to_delete;j];
                            end
                    end
                else % particle did not have enough neighbours to judge
                    ignored = [ignored;j];
                end
                
            end % for j-loop, over all particles
    end % method manual or not
    
end % for i, timepoint loop

if param.do_talk
    fprintf('Out of a total of %.0f vectors, %.0f were bad and %.0f did not have enough neighbours.\n',...
        num_beads, length(to_delete), length(ignored) );
end

%% Apply changes



%make a structure of the vectors that did not have enough neighbors
ignored = unique(ignored);
ign = struct([]);

to_delete = unique(to_delete);
deleted = struct([]);

for j = 1:time
    % note down which were ignored and which are to delete
    ign(j).r = d(j).r(ignored,:);
    ign(j).dr = d(j).dr(ignored,:);
    
    deleted(j).r = d(j).r(to_delete,:);
    deleted(j).dr = d(j).dr(to_delete,:);
    
    % deleting
    if ~param.keep_not_enough_neighbours
        %d(j).r(ignored,:) = nan(length(ignored),3);
        %d(j).dr(ignored,:) = nan(length(ignored),3);
        to_delete = [to_delete; ignored];
        to_delete = unique(to_delete); % remove dublicates
    end
    
    %d(j).r(to_delete,:) = nan(length(to_delete),3);
    %d(j).dr(to_delete,:) = nan(length(to_delete),3);
    
    d(j).r = local_remove_idx(d(j).r,to_delete);
    d(j).dr = local_remove_idx(d(j).dr,to_delete);
end

%Now display the figure if the user has specified 'graphics'
if param.do_show
    figure
    
   warning('This display is bugged!');
    quiver3((d(1).r(:,1)),(d(1).r(:,2)),d(1).r(:,3),(d(2).dr(:,1)),(d(2).dr(:,2)),d(2).dr(:,3),0,...
        'black','DisplayName','passed');
    hold on
    if ~isempty(deleted)
        quiver3(deleted(1).r(:,1),deleted(1).r(:,2),deleted(1).r(:,3),...
            deleted(2).dr(:,1),deleted(2).dr(:,2),deleted(2).dr(:,3),0,...
            'red','DisplayName','deleted');
    end
    if ~isempty(ign)
        quiver3((ign(1).r(:,1)),(ign(1).r(:,2)),ign(1).r(:,3),...
            (ign(2).dr(:,1)),(ign(2).dr(:,2)),ign(2).dr(:,3),0,...
            'green','DisplayName','too few neighbours');
    end
    legend;
    hold off
end





%OUTPUTS


good_disp = d;
if nargout>1
    varargout{1,1} = deleted;
end

end % function


%% ----- Sub-Functions -----
function angle_out = calc_angle(a,b)
% angle between two vectors (2D and 3D)
c = dot(a,b)/(norm(a)*norm(b));
angle_out = real(acos(c))*180/pi;

end

function mat_out = local_remove_idx(mat,idx_nr)

idx = ones(length(mat(:,1)),1);
idx(idx_nr) = 0; % remove all with these numbers

mat_out = mat ( logical ( idx ),: );


end
