%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% last update 11Jan2018, lne %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% model activation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 for turn off
% 1 for turn on

FE_Method=0;            % Diagonalization of the Hamiltonian (FEM)
PWE_Method=1;           % Plane Wave Expansion (PWE)

saveXYZ=0;
saveV=0;
savePSI=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=12;                   %% number of solution asked 

Fx=0;%5e7;              %% Electric field [V/m] in the x-direction
Fy=0;%5e6;              %% Electric field [V/m] in the y-direction
Fz=0;%5e7;              %% Electric field [V/m] in the z-direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Potential definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Three vectors (x, y and z) and one 3D-matrix V0 must be defined with homogeneous grid
% x, y and z [meter]
% V0 [eV]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx=64;                  %% Meshing point in x-direction
Ny=64;                  %% Meshing point in y-direction
Nz=64;                  %% Meshing point in z-direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose your between the next 3 potentials or build your own!

%Pot_InAs_GaAs
%Pot_InAs_AlInAs
Pot_GaN_AlN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% NOTHING TO CHANGE ANYMORE!!! %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Selection of the model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    NGx = 13;%Nx/2-1  ;    % number of harmonics % must be at least 2 times -1 smaller than Nz (smaller => faster)
    NGy = 15;%Ny/2-1  ;    % number of harmonics % must be at least 2 times -1 smaller than Nz (smaller => faster)
    NGz = 17;%Ny/2-1  ;    % number of harmonics % must be at least 2 times -1 smaller than Nz (smaller => faster)
    
    tic
    [E2,psi2] = Schroed3D_PWE_f(x,y,z,V0,Mass,n,Nx,Ny,Nz,NGx,NGy,NGz);
    display(strcat('-> PWE method =',num2str(toc),'sec'))
    if Fx~=0 || Fy~=0 || Fz~=0
      display('Warning: The PWE method is not the best for non-periodic potential')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E=nan(n,2);
E(1:length(E1),1)=E1;
E(1:length(E2),2)=E2;

display('=======================================')
display('Results:')
display('=======================================')
display(strcat('E(eV)='))
display(strcat(num2str(E)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FE_Method==1
c=0;
ii=0;
for i=1:n
    if i>18
      break
    end
    if i==1 || i==7 || i==13
      figure('Name','FEM method','position',[100 100 1600 800])
      c=c+1;
      ii=0;
    end
    ii=ii+1;
    
    subplot(2,3,ii,'fontsize',10)
    hold on;grid on;view (-38, 20);
    
    idz=find(z>0);idz=idz(1);
    pcolor(x*1e9,y*1e9,squeeze(V0(:,:,idz)))
    colormap(cool)
    %shading flat
    
    PSI=abs(psi1(:,:,:,i)).^2;
    
    p = patch(isosurface(x*1e9,y*1e9,z*1e9,PSI,max(PSI(:))/6));
    isonormals(x*1e9,y*1e9,z*1e9,PSI, p)
    set(p, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceLighting', 'gouraud')
    daspect([1,1,1])
    light ('Position', [1 1 5]);
    M=max([Mx My]);
    xlim([-1 1]*M/3*1e9)
    ylim([-1 1]*M/3*1e9)
    zlim([-1 1]*M/3*1e9)
    
    xlabel('x (nm)')
    ylabel('y (nm)')
    zlabel('z (nm)')
    title(strcat('E',num2str(i),'-E1=',num2str((E(i,1)-E(1,1))*1000,'%.1f'),'meV'))

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PWE_Method==1
c=0;
ii=0;
for i=1:n
    if i>18
      break
    end
    if i==1 || i==7 || i==13
      figure('Name','PWE method','position',[-1900 50 1850 1050])
      c=c+1;
      ii=0;
    end
    ii=ii+1;
    
%    subplot(2,3,ii,'fontsize',10)
    x0=0;y0=0;
    ww=0.35;
    hh=0.35;
    
    if i<4
      subplot('position',[x0+(i-1)*0.35 y0+0.6 ww hh])
    elseif i<7
      subplot('position',[x0+(i-4)*0.35 y0+0.05 ww hh])
    elseif i<10
      subplot('position',[x0+(i-7)*0.35 y0+0.6 ww hh])
    elseif i<13
      subplot('position',[x0+(i-10)*0.35 y0+0.05 ww hh])
    end
    hold on;%grid on;
    %axis off
    view (-38, 30);
    
    
    aa=[
    -sqrt(3)/2*a  -a/2   0
    -sqrt(3)/2*a  +a/2   0
    0             +a     0
    +sqrt(3)/2*a  +a/2   0
    +sqrt(3)/2*a  -a/2   0
    0             -a     0
    -sqrt(3)/2*a  -a/2   0
    ]*1e9;
    
    plot3(aa(:,1),aa(:,2),aa(:,3),'k','linewidth',2)
    
    bb=[
    +(QDh-h)/h*sqrt(3)/2*a  -(h-QDh)*a/h/2   QDh
    +(QDh-h)/h*sqrt(3)/2*a  +(h-QDh)*a/h/2   QDh
    0                       +(h-QDh)*a/h     QDh
    
    -(QDh-h)/h*sqrt(3)/2*a  +(h-QDh)*a/h/2   QDh
    -(QDh-h)/h*sqrt(3)/2*a  -(h-QDh)*a/h/2   QDh
    0                       -(h-QDh)*a/h     QDh
    +(QDh-h)/h*sqrt(3)/2*a  -(h-QDh)*a/h/2   QDh
    ]*1e9;
    
    plot3(bb(:,1),bb(:,2),bb(:,3),'k','linewidth',2)
    
    
    for j=1:length(aa(:,1))
      
      plot3([aa(j,1) bb(j,1)],[aa(j,2) bb(j,2)],[aa(j,3) bb(j,3)],'k','linewidth',2)
    
    end
    
    
    %contour(x*1e9,y*1e9,squeeze(V0(:,:,idz)),1,'linewidth',3,'linecolor','k')
    %idz=find(z>QDh);idz=idz(1)-1;
    %contour(x*1e9,y*1e9,squeeze(V0(:,:,idz)),1,'linewidth',3,'linecolor','k')
    
    %shading flat
    
    PSI=abs(psi2(:,:,:,i)).^2;
    
    p = patch(isosurface(x*1e9,y*1e9,z*1e9,PSI,max(PSI(:))/6));
    isonormals(x*1e9,y*1e9,z*1e9,PSI, p)
    set(p, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceLighting', 'gouraud')
    daspect([1,1,1])
    light ('Position', [1 1 5]);
    M=max([Mx My]);
    xlim([-1 1]*M/3*1e9)
    ylim([-1 1]*M/3*1e9)
    zlim([-1 1]*M/3*1e9)
    
    xlabel('x (nm)')
    ylabel('y (nm)')
    zlabel('z (nm)')
    title(strcat('E',num2str(i),'-E1=',num2str((E(i,2)-E(1,2))*1000,'%.1f'),'meV'))

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data save %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%