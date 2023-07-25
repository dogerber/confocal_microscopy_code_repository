function mat_out = yflip(mat_in)
% takes a [x,y,z] matrix and flips the y-axis

y = mat_in(:,2);
maxy = max(y); miny = min(y);

yt = y - maxy;
yt = abs(yt);
yt = yt + miny; 


mat_out = [mat_in(:,1),yt, mat_in(:,3:end)];

end