function slice_midplane(V,varargin)
V = single(V);
sz = size(V);
slice(V,...
    round(sz(1)/2),...
    round(sz(2)/2),...
    round(sz(3)/2));
xlabel('x'); ylabel('y'); zlabel('z');

end