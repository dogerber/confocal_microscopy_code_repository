function [pts_out,roi_out] = cut_selection(pts_in, roi)
%  [pts_out,roi_out] = cut_selection(pts_in, roi)
% Function that lets the user select an xy area who's points shall be kept
% pts_in is [x,y,z]
% or applies a roi (if you want to do the same to two sets of pts)
%



if nargin <2 % let user select roi
f1 = figure;
scatter(pts_in(:,1),pts_in(:,2),[],pts_in(:,3)); % we only cut in 2d

roi = drawpolygon();
%roi = drawrectangle();
end

% apply Roi
idx = roi.inROI(pts_in(:,1),pts_in(:,2));
roi_out = copy(roi);
pts_out = pts_in(idx,:);



if false % debugging
hold on
plot(pts_out(:,1),pts_out(:,2),'bo');
pause;
end

if exist('f1','var'); close(f1);end


end