function scatter3m(mat,varargin)
% like scatter3 but works with matrixinput [x,y,z]
%figure;

if nargin ==0
    load(uigetfile_to_fullpath());
    mat=pks;
end

if nargin <2
scatter3(mat(:,1),mat(:,2),mat(:,3),[],mat(:,3),...
    '.','MarkerFaceAlpha',0.6);
else
  scatter3(mat(:,1),mat(:,2),mat(:,3),[],mat(:,3),...
    varargin{:});
end

xlabel('x'); ylabel('y');
axis equal
shg



end