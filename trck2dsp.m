function disp=trck2dsp(tracks)

% PURPOSE
% Creates a useful structure for particle positions and displacements from 
% the output of track.m (by Blair and Dufresne) or Tracker.m 
% (by Boltyanskiy and Dufresne)

% INPUT

% tracks        4xN matrix where N is the number of particles found at both
                % time points such that tks(i,1) is the x-coordinate of the 
                % ith particle; tks(i,2) is the y-coordinate of the ith 
                % particle; tks(i,3) is the image number, and tks(i,4) is v
                % the particle identification number (partcle ID) 

% OUTPUT  

% disp:         % a structure array of N elements, where disp(t).r(i,:) are 
                % the coordinates of the ith particle at time t, 
                % disp(t).dr(i,:) is the displacment of the ith particle at 
                % time t, relative to the first time point.


% NOTE:
% this code rejects all particles that are not present at all time points.


[r, c] = size(tracks);

if r == 0
    error('No particles were tracked')
end

if c==4
    d=2;     % number of dimensions   
else
    d=3;     % number of dimensions
end

ti=min(tracks(:,d+1));  % first frame number

% make first frame number = 1
if ti~=1
    tracks(:,d+1)=tracks(:,d+1)-ti+1;
end
tf=max(tracks(:,d+1));  % final frame number

np = max(tracks(:,d+2));    % number of particles   


% reject all particles not present in all time steps (new DG 220203)
idx_all_tp = zeros(np,1);
for n_part = 1:np
    if sum(tracks(:,d+2)==n_part)==tf
      %  tracks(tracks(:,d+2)==n_part,:)
%         idx_all_tp(1:20)
%         n_part
         idx_all_tp(n_part) = 1;
%         beep;
    end
end
idx_keep= ismember(tracks(:,d+2),find(idx_all_tp));
tracks = tracks(logical(idx_keep),:);


for t=1:tf
    tks = tracks(tracks(:,d+1)==t,:);
    if t==1
        r0=tks(:,1:d);              % particle positions in the first frame
        tks0 = tks;
    end
   disp(t).r = tks(:,1:d);          % particle positions in frame t
   disp(t).dr = tks(:,1:d)-r0;      % particle displacements in frame t
end
