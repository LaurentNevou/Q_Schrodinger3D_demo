%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% last update 9Jan2018, lne %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% model activation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 for turn off
% 1 for turn on

FE_Method=0;            % Diagonalization of the Hamiltonian (FEM)
PWE_Method=1;           % Plane Wave Expansion (PWE)

saveXY=0;
saveV=0;
savePSI=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=6;                   %% number of solution asked 
Mass = 0.043;           %% effective mass, constant over all the structure...
Fx=0;%5e7;              %% Electric field [V/m] in the x-direction
Fy=0;%5e6;              %% Electric field [V/m] in the y-direction
Fz=0;%5e7;              %% Electric field [V/m] in the z-direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Potential definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Three vectors (x, y and z) and one 3D-matrix V0 must be defined with homogeneous grid
% x, y and z [meter]
% V0 [eV]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx=64;                  %% Meshing point in x-direction
Ny=64;                  %% Meshing point in y-direction
Nz=32;                  %% Meshing point in z-direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose your between the next 3 potentials or build your own!

Pot_InAs_GaAs
%Pot_InAs_AlInAs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vb=0.55;                 %% Potential barrier height[eV]
V0=(idx)*0 + (1-idx)*Vb ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% NOTHING TO CHANGE ANYMORE!!! %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V0=(Fx*X)+V0;        % adding the electric field Fx to the potential in the x-direction
V0=(Fy*Y)+V0;        % adding the electric field Fy to the potential in the y-direction
V0=(Fz*Z)+V0;        % adding the electric field Fz to the potential in the y-direction
V0=V0-min(min(min(V0)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Selection of the model %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E1=[];E2=[];
display('=======================================')

if FE_Method==1
    tic
    if length(x)*length(y)*length(z)>1e4
      N=length(x)*length(y)*length(z);
      display(strcat('Warning: Take care, H=',num2str(N),'x',num2str(N),'elements'))
    end
    [E1,psi1] = Schroed3D_FEM_f(x,y,z,V0,Mass,n);  % m=cste
    display(strcat('-> Finite Elements method =',num2str(toc),'sec'))
end

if PWE_Method==1
    Nx = 64 ;        % number of points on the x grid % has to be a power of 2 (32,64,128,256,512,...) (smaller => faster)
    Ny = 64 ;        % number of points on the y grid % has to be a power of 2 (32,64,128,256,512,...) (smaller => faster)
    Nz = 64 ;        % number of points on the y grid % has to be a power of 2 (32,64,128,256,512,...) (smaller => faster)
    NGx = 9;%Nx/2-1  ;    % number of harmonics % must be at least 2 times -1 smaller than Nz (smaller => faster)
    NGy = 7;%Ny/2-1  ;    % number of harmonics % must be at least 2 times -1 smaller than Nz (smaller => faster)
    NGz = 11;%Ny/2-1  ;    % number of harmonics % must be at least 2 times -1 smaller than Nz (smaller => faster)
    
    tic
    [E2,psi2] = Schroed3D_PWE_f(x,y,z,V0,Mass,n,Nx,Ny,Nz,NGx,NGy,NGz);
    display(strcat('-> PWE method =',num2str(toc),'sec'))
    if Fx~=0 || Fy~=0 || Fz~=0
      display('Warning: The PWE method is not the best for non-periodic potential')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E=nan(n,2);
E(1:length(E1),1)=E1;
E(1:length(E2),2)=E2;

display('=======================================')
display('Results:')
display('=======================================')
display(strcat('E(eV)='))
display(strcat(num2str(E)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure('Name','Potential','position',[100 100 1200 400])
figure('Name','Potential','position',[-1900 50 1600 500])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,1,'fontsize',15)
hold on;grid on;

pcolor(y*1e9,z*1e9,squeeze(V0(:,round(end/2),:))')
colormap(cool)
colorbar

xlim([-1 1]*My/2*1e9)
ylim([-1 1]*My/2*1e9)

xlabel('x (nm)')
ylabel('z (nm)')
zlabel('Energy (eV)')
%title(strcat('Potential, Fx=',num2str(Fx,'%.1e'),'[V/m], Fy=',num2str(Fy,'%.1e'),'[V/m]'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,2,'fontsize',15)
hold on;grid on;

pcolor(x*1e9,y*1e9,squeeze(V0(:,:,round(end/2))))
colormap(cool)
colorbar

xlim([-1 1]*My/2*1e9)
ylim([-1 1]*My/2*1e9)

xlabel('x (nm)')
ylabel('y (nm)')
%title(strcat('Potential, Fx=',num2str(Fx,'%.1e'),'[V/m], Fy=',num2str(Fy,'%.1e'),'[V/m]'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FE_Method==1
c=0;
ii=0;
for i=1:n
    if i>45
      break
    end
    if i==1 || i==16 || i==31
      figure('Name','FEM method','position',[100 100 1600 900])
      c=c+1;
      ii=0;
    end
    ii=ii+1;
    
    subplot(3,5,ii,'fontsize',10)
    hold on
    
    pcolor(x*1e9,y*1e9,(psi1(:,:,i)) )  
    contour(x*1e9,y*1e9,V0,1,'linewidth',3,'linecolor','w')
    
    xlabel('x (nm)')
    ylabel('y (nm)')
    title(strcat('E',num2str(i),'=',num2str(E(i,1)*1000,'%.1f'),'meV'))
    %axis equal
    shading flat
    colormap(jet)
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PWE_Method==1
c=0;
ii=0;
for i=1:n
    if i>18
      break
    end
    if i==1 || i==7 || i==13
      figure('Name','PWE method','position',[-1900 50 1850 1000])
      c=c+1;
      ii=0;
    end
    ii=ii+1;
    
    subplot(2,3,ii,'fontsize',10)
    hold on;grid on;view (-38, 20);
    
    xslice = 0; yslice =0; zslice = 0.4;
    slice(X*1e9,Y*1e9,Z*1e9,V0,[],[],zslice)
    colormap(cool)
    %shading flat
    
    
    p = patch(isosurface(x*1e9,y*1e9,z*1e9,psi2(:,:,:,i),max(psi2(:))/6));
    isonormals(x*1e9,y*1e9,z*1e9,psi2(:,:,:,i), p)
    set(p, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceLighting', 'gouraud')
    daspect([1,1,1])
    light ("Position", [1 1 5]);
    M=max([Mx My]);
    xlim([-1 1]*M/3*1e9)
    ylim([-1 1]*M/3*1e9)
    zlim([-1 1]*M/3*1e9)
    
    xlabel('x (nm)')
    ylabel('y (nm)')
    title(strcat('E',num2str(i),'-E1=',num2str((E(i,2)-E(1,2))*1000,'%.1f'),'meV'))

    end
end

break
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Data save %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveXYZ==1;
    x=x';y=y';z=z';
    save('data_x.txt','x','-ascii')
    save('data_y.txt','y','-ascii')
    save('data_z.txt','z','-ascii')
end

if saveV==1;
    save('data_V.txt','V0','-ascii')
end

if savePSI==1;
    
  if FE_Method==1
    for i=1:n
      M1 = psi1(:,:,:,i);
      save(strcat('data_psi',num2str(i),'_FEM.txt'),'M1','-ascii')
    end
  end
  if PWE_Method==1
    for i=1:n
      M2 = real(psi2(:,:,:,i));
      save(strcat('data_psi',num2str(i),'_PWE.txt'),'M2','-ascii')
    end
  end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%