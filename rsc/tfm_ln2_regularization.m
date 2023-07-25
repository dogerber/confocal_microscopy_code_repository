function MM = tfm_ln2_regularization(nx,param)
%%tfm_ln2_regularization
%
% Returns the M matrix, needs interpolation grid size nx
%
%
% original file: regularising_stresse_in_regions_image_mk2data.m
%% Dominic Gerber, 2020, original by Dr. Robert Style

% Code to calculate the matrix MM such that if the stresses in real space are
% written in one vector:
%
% S=[sx_11; sx_12;...;sx_1N;...;sx_N1;...;sx_NN;sy_11;
% sy_12;...;sy_1N;...;sy_N1;...;sy_NN;sz_11;
% sz_12;...;sz_1N;...;sz_N1;...;sz_NN];
%
% and the fourier transformed displacements are written in a single vector:
%
% Uf=[uf_11; uf_12;...;uf_1N;...;uf_N1;...;uf_NN;vf_11;
% vf_12;...;vf_1N;...;vf_N1;...;vf_NN;wf_11;
% wf_12;...;wf_1N;...;wf_N1;...;wf_NN];
%
% Then Uf=M*S. By inverting this, we can get the stresses from the fourier
% transformed displacements. We can also easily apply L2 regularisation to
% the solution, to set parts of the stress field to zero.
%
% NB in fourier space, we start from the fftshift of fft2(u), so Uf is
% [fftshift(fft2(u)(:);fftshift(fft2(v)(:),fftshift(fft2(w)(:)]

param.do_debug = true;
param.do_show = true;

% First calculate the matrix that gives the fourier transform
%nx=45; % Number of points in x. Must be odd
ny=nx; % Number of points in y. Only checked when nx=ny; Must be square
if ~mod(nx,2); error('nx must be odd!'); end
x=[1:1:nx]-1;
y=[1:1:ny]-1;
[X,Y]=meshgrid(x,y);
X=X(:)';
Y=Y(:)';

kx=[-(nx-1)/2:1:(nx-1)/2]/nx;
ky=[-(nx-1)/2:1:(nx-1)/2]/nx;
[Kx,Ky]=meshgrid(kx,ky);
Kx=Kx(:);
Ky=Ky(:);

F=exp(-2*pi*1i*Kx*X).*exp(-2*pi*1i*Ky*Y); % This does a straight fourier transform (have double checked)
% To do a 2d FFT on a matrix a,
% A=F*a(:);
% A=reshape(A,nx,ny)

% Now we need to make a block version of this, because we need the fourier
% transform of S
F3=blkdiag(F,F,F);
% And rotate the rows, so that the output is sf=(sx_11,sy_11,sz_11,...), as
% this is easier to apply the Qinv matrix to
for i=1:3*nx*ny
    if mod(i,3)==1
        inds(i)=(i-1)/3+1;
    else
        inds(i)=inds(i-1)+nx*ny;
    end
end
f=F3(inds,:);

%%
% Now make the Q inv matrix to go from sf to uf=(ux_11,uy_11,uz_11,...)

thick=param.substrate_thickness;
EM=param.substrate_youngs_modulus ;
nu=param.nu;
dx=param.dx;


Q = calcQ(thick,thick,EM,nu,nx,dx,3);

qhold=zeros(3,3,nx^2);
for i=1:3
    for j=1:3
        qhold(i,j,:)=[Q(i,j).Qinv(:)];
    end
end
qinv=num2cell(qhold,[1,2]);
qinv=blkdiag(qinv{:});
%switch the rows so that it corresponds to Uf instead of uf
Qinv=qinv([1:3:3*nx*ny-2,2:3:3*nx*ny-1,3:3:3*nx*ny],:);

%% Now make an inverse fourier transform matrix in case we need it
if false
%nx=45; % Number of points in x. Must be odd
%ny=45; % Number of points in y. Only checked when nx=ny;
x=[1:1:nx]-1;
y=[1:1:ny]-1;
[X,Y]=meshgrid(x,y);
X=X(:);
Y=Y(:);

kx=[-(nx-1)/2:1:(nx-1)/2]/nx;
ky=[-(nx-1)/2:1:(nx-1)/2]/nx;
[Kx,Ky]=meshgrid(kx,ky);
Kx=Kx(:)';
Ky=Ky(:)';

Finv=1/nx^2*exp(2*pi*1i*X*Kx).*exp(2*pi*1i*Y*Ky); % Checked, and looks good
end


%% Finally make M matrix
MM=Qinv*f;
%MMinv = inv(MM);



end


