function[E,psi]=Schroed3D_PWE_f(x,y,z,V0,Mass,n,Nx,Ny,Nz,NGx,NGy,NGz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant [J.s]
hbar=h/(2*pi);
e=1.602176487E-19;              %% electron charge [C]
m0=9.10938188E-31;              %% electron mass [kg]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Interpolation on a grid that have 2^N points %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NGx = 2*floor(NGx/2);           %% round to lower even number
NGy = 2*floor(NGy/2);           %% round to lower even number
NGz = 2*floor(NGz/2);           %% round to lower even number

[X,Y,Z] = meshgrid(x,y,z);
xx=linspace(x(1),x(end),Nx);
yy=linspace(y(1),y(end),Ny);
zz=linspace(z(1),z(end),Nz);
[XX,YY,ZZ] = meshgrid(xx,yy,zz);

V=interp3(X,Y,Z,V0,XX,YY,ZZ);

dx=x(2)-x(1);
dxx=xx(2)-xx(1);
dy=y(2)-y(1);
dyy=yy(2)-yy(1);
dz=z(2)-z(1);
dzz=zz(2)-zz(1);
Ltotx=xx(end)-xx(1);
Ltoty=yy(end)-yy(1);
Ltotz=zz(end)-zz(1);

[XX,YY,ZZ] = meshgrid(xx,yy,zz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Building of the potential in Fourier space %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vk = fftshift(fftn(V))*dxx*dyy*dzz/Ltotx/Ltoty/Ltotz;
Vk =Vk(Ny/2-NGy+1:Ny/2+NGy+1 , Nx/2-NGx+1:Nx/2+NGx+1 , Nz/2-NGz+1:Nz/2+NGz+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Reciprocal lattice vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gx = (-NGx/2:NGx/2)'*2*pi/Ltotx;
Gy = (-NGy/2:NGy/2)'*2*pi/Ltoty;
Gz = (-NGz/2:NGz/2)'*2*pi/Ltotz;

NGx=length(Gx);
NGy=length(Gy);
NGz=length(Gz);
NG=NGx*NGy*NGz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Building Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx_x = reshape(1:NGx, [1 NGx 1]);
idx_x = repmat(idx_x, [NGy 1 NGz]);
idx_x = idx_x(:);

idx_y = reshape(1:NGy, [NGy 1 1]);
idx_y = repmat(idx_y, [1 NGx NGz]);
idx_y = idx_y(:);

idx_z = reshape(1:NGz, [1 1 NGz]);
idx_z = repmat(idx_z, [NGy NGx 1]);
idx_z = idx_z(:);

%idx_X = (idx_x-idx_x') + NGx;      %% work only in Octave
%idx_Y = (idx_y-idx_y') + NGy;      %% work only in Octave
%idx_Z = (idx_z-idx_z') + NGz;      %% work only in Octave
idx_X = (repmat(idx_x,[1 NG])-repmat(idx_x',[NG 1])) + NGx;     %% work in Octave and Matlab
idx_Y = (repmat(idx_y,[1 NG])-repmat(idx_y',[NG 1])) + NGy;     %% work in Octave and Matlab
idx_Z = (repmat(idx_z,[1 NG])-repmat(idx_z',[NG 1])) + NGz;     %% work in Octave and Matlab

idx = sub2ind(size(Vk), idx_Y(:), idx_X(:), idx_Z(:));
idx = reshape(idx, [NG NG]);

GX = diag(Gx(idx_x));
GY = diag(Gy(idx_y));
GZ = diag(Gz(idx_z));

D2 = GX.^2 + GY.^2 + GZ.^2 ;

H =  hbar^2/(2*m0*Mass)*D2  +  Vk(idx)*e ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Solving Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = sparse(H);
[psik, Ek] = eigs(H,n,'SM');
E = diag(Ek)  / e;
%E=abs(E);
E=real(E);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Transforming & Scaling the waves functions %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:n
    PSI = reshape(psik(:,j),[NGy,NGx,NGz]);
    PSI = invFFT3D(PSI,Ny,Nx,Nz)/(dxx*dyy*dzz) ;
    psi_temp = interp3(XX,YY,ZZ,PSI,X,Y,Z);
    psi(:,:,:,j) = psi_temp / max(psi_temp(:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here is a small patch due to differences between Octave and Matlab
% Matlab order the eigen values while Octave reverse it

if E(1)>E(2)
  psi=psi(:,:,:,end:-1:1);
  E=E(end:-1:1);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Vxyz] = invFFT3D(Vk3D,Ny,Nx,Nz)

Nkx=length(Vk3D(1,:,1));
Nky=length(Vk3D(:,1,1));
Nkz=length(Vk3D(1,1,:));

Nx1=Nx/2-floor(Nkx/2);
Nx2=Nx/2+ceil(Nkx/2);
Ny1=Ny/2-floor(Nky/2);
Ny2=Ny/2+ceil(Nky/2);
Nz1=Nz/2-floor(Nkz/2);
Nz2=Nz/2+ceil(Nkz/2);

Vk3D00=zeros(Ny,Nx,Nz);
Vk3D00( Ny1+1:Ny2 , Nx1+1:Nx2 , Nz1+1:Nz2)=Vk3D;
Vxyz=ifftn(ifftshift(Vk3D00));

end
