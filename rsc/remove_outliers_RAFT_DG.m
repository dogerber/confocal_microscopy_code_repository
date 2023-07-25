function tks2=remove_outliers_RAFT_DG(tks,varargin)
%
% e.g. tks2=remove_outliers_RAFT(tks,'n',20,'abs_diff',5);
% tks is [x,y,*z,time,particleID]

% code to identify outliers by comparing the movement of points to that of
% nearby points. If these significantly differ - either by 3 standard
% deviations, or by a fixed amount - then the track is split at the point
% where the deviation happens.
% If you find this useful please consider citing DOI: 10.5281/zenodo.4884065
% INPUT
% tks is the output of track_RAFT
% OPTIONS:
% 'min_track_length': the minimum number of timepoints that individual
% particles should contain. Split tracks that are shorter than this will be
% thrown away. Must be 2 or greater. Default is 2
%
% 'dim' is the number of dimensions the data has. Default is the number of
% columns in the tracking structure - 2
%
% 'radius' is the radius around a particle to compare it to neighbours 2 *
% average spacing
%
% 'abs_diff' is the maximum difference that a particles displacement can
% have from the average displacement of its nearest neighbours in any
% cartesian direction. Set to nan if you don't wish to use this option.
% Default is nan
%
% 'angle_diff' maximum allowed angle difference between a displacement
% vector and its neighbours
%
% 'min_number_neighbours' minimum number of neighbours in distance 'radius'
% to not throw a particle away
%
% 'show_result' displays all passed and removed vectors for each timepoint
%
%
% Warning! If using this on data where you are measuring something like
% diffusivity, make sure it is working like it is meant to (i.e. it is not
% breaking up real tracks), otherwise it will bias your diffusivity data,
% by throwing away bigger jumps.
%
% Outputs
% trks2: the trks data with the particle tracks split at points where bad
% tracking is detected
%
% Created by Rob Style in lockdown, completed 31/05/2021
% Updated by Dominic Gerber 02/08/2022
%     added: angle and min number of neighbours criterion;
%     changed: neighbours by distance 'radius' and not by n_number nearest



%% Input handling

inputExist = find(cellfun(@(x) strcmpi(x, 'min_track_length') , varargin));
if inputExist
    min_track_length = varargin{inputExist+1};
else
    min_track_length=2;
end

inputExist = find(cellfun(@(x) strcmpi(x, 'dim') , varargin));
if inputExist
    dim = varargin{inputExist+1};
else
    dim = size(tks,2)-2; % tks is [x,y,*z,time,particlenumber]
end

inputExist = find(cellfun(@(x) strcmpi(x, 'abs_diff') , varargin));
if inputExist
    abs_diff = varargin{inputExist+1};
else
    abs_diff=0; % if zero, it is not applied
end

inputExist = find(cellfun(@(x) strcmpi(x, 'angle_diff') , varargin));
if inputExist
    angle_diff = varargin{inputExist+1};
else
    angle_diff=nan;
end

inputExist = find(cellfun(@(x) strcmpi(x, 'radius') , varargin));
if inputExist
    radius = varargin{inputExist+1};
else
    radius=2* calculate_average_particle_spacing(tks(tks(:,end-1)==min(tks(:,end-1)),1:3)); % 3* average spacing of first timepoint
end

inputExist = find(cellfun(@(x) strcmpi(x, 'min_number_neighbours') , varargin));
if inputExist
    min_number_neighbours = varargin{inputExist+1};
else
    min_number_neighbours= 5;
end

inputExist = find(cellfun(@(x) strcmpi(x, 'show_result') , varargin));
if inputExist
    show_result = varargin{inputExist+1};
else
    show_result=false;
end

ts=unique(tks(:,end-1));


%% Main loop over timepoints


for i=1:length(ts)
    if ~isempty(tks(tks(:,end-1)==ts(i)+1,:))
        if false % old
        pno_diff=tks(:,end) - circshift(tks(:,end),-1); % same particle number
        t_diff=tks(:,end-1) - circshift(tks(:,end-1),-1); % same timepoint
        inds_pair_ts1= pno_diff==0 & t_diff==-1 & tks(:,end-1)==ts(i);
        
        pno_diff=tks(:,end) - circshift(tks(:,end),1);
        t_diff=tks(:,end-1) - circshift(tks(:,end-1),1);
        inds_pair_ts2=pno_diff==0 & t_diff==1 & tks(:,end-1)==ts(i)+1;
        else
           idx_t1 = find(tks(:,end-1)== ts(i));
           idx_t2= find(tks(:,end-1)== ts(i)+1);
            [~,inds_pair_ts1,inds_pair_ts2] = intersect(tks(idx_t1, end),tks( idx_t2, end));
            inds_pair_ts1=idx_t1(inds_pair_ts1);
            inds_pair_ts2=idx_t2(inds_pair_ts2);
            
        end
        
        inds_bad = local_find_outliers_RAFT(tks(inds_pair_ts1,1:dim),tks(inds_pair_ts2,1:dim),'radius',radius,'abs_diff',abs_diff,...
            'min_number_neighbours',min_number_neighbours,...
            'angle_diff',angle_diff,...
            'show_result',show_result);
        
      %old  inds=inds_pair_ts1;
      %old  inds(inds_pair_ts1)=inds_bad; % This has a list of rows which correspond to bad pairs
        
        % Now take bad tracks and split the tracks there and allocate the new split track a new number
       %old  bad_pnos=tks(inds,end); % read out track-id
        bad_pnos=tks(inds_pair_ts1(inds_bad),end);
        
        if length(ts) > min_track_length  % do try and split [slow for large data set]
            
            for j=1:length(bad_pnos)
                if (sum(tks(:,end)==bad_pnos(j))-1<min_track_length) % not enough tp of this track left, no need to split it
                    tks(tks(:,end)==bad_pnos(j),:)=[];
                else % spilit the track
                    if sum(tks(:,end)==bad_pnos(j) & tks(:,end-1)>ts(i),'all')>=min_track_length % if there are less than min_track_length particles left in the track after the split, delete. Otherwise assign the rest of the track a new number
                        tks(tks(:,end)==bad_pnos(j) & tks(:,end-1)>ts(i),end)=max(tks(:,end))+1;
                    else % to short, completly delete this track
                        tks(tks(:,end)==bad_pnos(j) & tks(:,end-1)>ts(i),:)=[];
                    end
                    
                    % if there are less than min_track_length particles left in the track before the split, delete
                    if sum(tks(:,end)==bad_pnos(j) & tks(:,end-1)<=ts(i),'all')<min_track_length
                        tks(tks(:,end)==bad_pnos(j) & tks(:,end-1)<=ts(i),:)=[];
                    end
                end
            end
            %tks(ismember(tks(:,end),tks(inds,end)),:)=[];
        else % do not split
            warning('splitting of tracks disabled');
            tks(tks(bad_pnos,end),:)=[];
        end
        
    end
end

tks2=tks;
end

%% --- Sub-Functions ---
function inds_bad=local_find_outliers_RAFT(pks1,pks2,varargin)

% code to identify outliers. Idea is that for every point, find the
% n nearest points, and work out what the average displacement vector is for
% them. Then compare and see if they're within 3 standard deviations of
% each other
% INPUTS:
% pks1,pks2 are two arrays that are the same size, and have matched data
% sets. I.e. pks1 contains the positions of points at one time, and pks2
% contains points at the second time, where points in the same row
% correspond to each other
% OPTIONS:
% 'abs_diff' is the maximum difference that a particles displacement can
% have from the average displacement of its nearest neighbours in any
% cartesian direction
inputExist = find(cellfun(@(x) strcmpi(x, 'radius') , varargin));
if inputExist
    radius = varargin{inputExist+1};
else
    error('Need radius for this.');
end

inputExist = find(cellfun(@(x) strcmpi(x, 'min_number_neighbours') , varargin));
if inputExist
    min_number_neighbours = varargin{inputExist+1};
else
    min_number_neighbours= 5;
end

inputExist = find(cellfun(@(x) strcmpi(x, 'abs_diff') , varargin));
if inputExist
    abs_diff = varargin{inputExist+1};
else
    abs_diff=nan;
end

inputExist = find(cellfun(@(x) strcmpi(x, 'angle_diff') , varargin));
if inputExist
    angle_diff = varargin{inputExist+1};
else
    angle_diff=nan;
end

inputExist = find(cellfun(@(x) strcmpi(x, 'show_result') , varargin));
if inputExist
    show_result = varargin{inputExist+1};
else
    show_result=false;
end

% Start by working out who to compare to
%near_neighb_inds=zeros(length(pks1),n_consider+1);

near_neighb_inds = rangesearch(pks1,pks1,radius);


% Get rid of first column, as this is just the particle itself.
%[~,ids]=mink(D11,n+1,2);
%ids=near_neighb_inds;
%ids(:,1)=[];

% Measure metrics
disp=pks2-pks1; % displacements

disp_mean=zeros(size(pks1,1),size(pks1,2));
disp_std=zeros(size(pks1,1),size(pks1,2));

for i=1:size(pks1,1)
    disp_temp=pks2(near_neighb_inds{i}(2:end),:)-pks1(near_neighb_inds{i}(2:end),:);
    
    disp_mean(i,:)=mean(disp_temp,1);
    disp_std(i,:)=std(disp_temp,1);
    
    avg_disp_vector_i(i,:) = mean(disp_temp,1);
    try
        angle_difference(i,:) = calc_angle(disp(i,:),avg_disp_vector_i(i,:));
    catch
        beep;
    end
end

[~,nr_neighbours] = cellfun(@size,near_neighb_inds);

% apply cirterion what is considered bad
bads_abs_diff = zeros(length(pks1(:,1)),1);
bads_nr_neighbours = zeros(length(pks1(:,1)),1);
bads_angle_diff = zeros(length(pks1(:,1)),1);

if ~isnan(abs_diff)
    bads_abs_diff=[abs(disp_mean-disp)>3*disp_std] | [abs(disp_mean-disp)>abs_diff];
end
bads_abs_diff=any(bads_abs_diff,2); %bads(:,1) | bads(:,2) | bads(:,3);

if ~isnan(min_number_neighbours)
    bads_nr_neighbours = nr_neighbours < min_number_neighbours;
end

if ~isnan(angle_diff)
    bads_angle_diff= angle_difference > angle_diff;
end

inds_bad = any([bads_abs_diff,bads_nr_neighbours, bads_angle_diff],2);

if show_result % for debugging
    figure('Units','normalized','Position',[0.6,0.1,0.3,0.8]); hold on
    subplot(3,1,1); hold on
    local_d_plotter(pks1,pks2,logical(bads_abs_diff),'r','bads_abs_diff',1)
    local_d_plotter(pks1,pks2,logical(bads_nr_neighbours),'g','bads_nr_neighbours',1)
    local_d_plotter(pks1,pks2,logical(bads_angle_diff),'b','bads_angle_diff',1)
    
    local_d_plotter(pks1,pks2,logical(~inds_bad),'black','passed',1)
    
    
    
    fprintf('removing... \n \t %d for abs diff, \n \t %d for nr_neighbours, \n \t %d for angle diff, \n \t %d passed\n',...
        sum(bads_abs_diff),sum(bads_nr_neighbours),sum(bads_angle_diff), sum(~inds_bad));
     view(20,10);
    title('Result of remove_outliers_RAFT_DG','Interpreter','none');
    legend('Interpreter','none','location','best');
    subtitle('autoscaled');
    
    subplot(3,1,2); hold on
    local_d_plotter(pks1,pks2,logical(bads_abs_diff),'r','bads_abs_diff',0)
    local_d_plotter(pks1,pks2,logical(bads_nr_neighbours),'g','bads_nr_neighbours',0)
    local_d_plotter(pks1,pks2,logical(bads_angle_diff),'b','bads_angle_diff',0)
    
    local_d_plotter(pks1,pks2,logical(~inds_bad),'black','passed',0)
    
        % show one example of a particle and it's neighbours
        if true
    part_id_for_neighbours = 10;
    inds_near_id_for_neighbours = near_neighb_inds{part_id_for_neighbours};

    plot3(pks1(inds_near_id_for_neighbours,1),pks1(inds_near_id_for_neighbours,2),pks1(inds_near_id_for_neighbours,3),...
        'bo','DisplayName','example neighbours');
        plot3(pks1(part_id_for_neighbours,1),pks1(part_id_for_neighbours,2),pks1(part_id_for_neighbours,3),...
            'ro','DisplayName','example center');
        end
 
    view(20,10); shg;
    xlabel('x'); ylabel('y'); zlabel('z');
    subtitle('not autoscaled','Interpreter','none');
    legend('Interpreter','none','location','best');
    axis equal
    
        subplot(3,1,3); hold on
    local_d_plotter(pks1,pks2,logical(~inds_bad),'black','bads_angle_diff',0)
    xlabel('x'); ylabel('y'); zlabel('z');
     view(20,10); shg;
    title('Passed arrows')
    subtitle('not autoscaled','Interpreter','none');
    legend('Interpreter','none','location','best');
    axis equal
   
end



end


function local_d_plotter(pks1,pks2,idx,colorspec,displayName,scale_factor)
if sum(idx) >0
    %scale_factor = 1; % make 0 for actaul length of quiver
    diff = pks2(idx,:)-pks1(idx,:);
    
    quiver3(pks1(idx,1),pks1(idx,2),pks1(idx,3),...
        diff(:,1),diff(:,2),diff(:,3), scale_factor,...
        colorspec,'DisplayName',displayName);
end

end



function angle_out = calc_angle(a,b)
% angle between two vectors (2D and 3D)
c = dot(a,b)/(norm(a)*norm(b));
angle_out = real(acos(c))*180/pi;

end

function average_spacing = calculate_average_particle_spacing(pks)

if length(pks(1,:))==3
    vol = (max(pks(:,1))-min(pks(:,1)))*((max(pks(:,2))-min(pks(:,2))))*((max(pks(:,1))-min(pks(:,2))));
    
    average_spacing = (vol/length(pks))^(1/3);
else
    vol = (max(pks(:,1))-min(pks(:,1)))*((max(pks(:,2))-min(pks(:,2))));
    
    average_spacing = (vol/length(pks))^(1/2);
    
end
end
