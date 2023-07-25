function viz_tracks(trks,varargin)
% to explore the trks format produced e.g. in track_RAFT
% track is [x,y, t, id] or [x,y,z,t,id]
%

[unique_ids,idx_ids,idx_inverse] = unique(trks(:,end));
unique_t = unique(trks(:,end-1));
 
fprintf('trks has %d particles in %d unique tracks and %d timepoints. \n',length(trks(:,1)),length(unique_ids),length(unique_t));

if unique_t>2
counts_ids_occurence = accumarray(idx_inverse,1);
figure;
histogram(counts_ids_occurence);
xlabel('length of track'); ylabel('# tracks with this length');
end

if false
colors_tp = colors_fader([1,0,0],[0,0,1],length(unique_t));
else
    colors_tp = parula(length(unique_t));
end



%% quiver3 with color for timepoint
%if nargin<2; figure;end

tp_idx = 1;
for tp=unique_t(2:end)'
    
    fprintf('timepoint %d \n', tp);
    
    trks0 = trks( trks(:,end-1)==tp-1, :);
    trks1= trks( trks(:,end-1)==tp, :);
    
    common_idx = intersect(trks0(:,end),trks1(:,end));
    
    % only keep common idx
    trks0 = trks0( any(trks0(:,end) == common_idx',2),:);
    trks1 = trks1( any(trks1(:,end) == common_idx',2),:);
    
    %figure; plot(trks0(:,end),'.'); hold on; plot(trks1(:,end),'o'); % must be on top of each other
    
    if length(trks0)~=length(trks1); error('They must have the same length now, unless something is wrong.'); end
    
    pos = trks0(:,[1:3]);
    dif = trks1(:,[1:3])-trks0(:,[1:3]);
    
    hold on;
    if nargin<2
    quiver3(pos(:,1),pos(:,2),pos(:,3),...
        dif(:,1),dif(:,2),dif(:,3),0,...
        'Color',colors_tp(tp_idx,:));
    else
            quiver3(pos(:,1),pos(:,2),pos(:,3),...
        dif(:,1),dif(:,2),dif(:,3),0,...
       varargin{:});
    end
    
    % plot all trks xyz (to make very small or 0 displacements visible)
    if true
       hold on;
       mrker_size = 2;
       plot3(trks0(:,1),trks0(:,2),trks0(:,3),'black.','MarkerSize',mrker_size);
       plot3(trks1(:,1),trks1(:,2),trks1(:,3),'blackx','MarkerSize',mrker_size);
    end
    
    tp_idx = tp_idx+1;
    
end




%% ----- old -----
if false
%% one track at a time
 

figure; hold on
for i = 1:1:length(unique_ids)
    
    
    sub_trks = trks(trks(:,end)==unique_ids(i),:);
    
    if length(trks(1,:))== 4 % 2d
        if length(sub_trks)>1
            for j=1:length(sub_trks)
                plot(sub_trks(:,1),sub_trks(:,2),'.-','Color',plot_random_color);
            end
        end
        
    else % 3D
         plot3(sub_trks(:,1),sub_trks(:,2),sub_trks(:,3),...
             '.-','Color',plot_random_color);
    end    
end


%% one timepoint at a time

figure
if length(trks(1,:))== 4 % 2d
scatter(trks(:,1),trks(:,2),[],trks(:,3),'.');
else %3D
    scatter3(trks(:,1),trks(:,2),trks(:,3),[],trks(:,end-1),'.');
    title('color for timepoint');
end
colorbar

%% two timepoints at a time
for t=1:length(unique_t)
    
    
end
end