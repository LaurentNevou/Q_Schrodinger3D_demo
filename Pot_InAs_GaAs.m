%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vb=0.55;                %% Potential barrier height[eV]
Mass = 0.043;           %% effective mass, constant over all the structure...
Mx=50e-9;               %% map X [m]
My=50e-9;               %% map Y [m]
Mz=10e-9;               %% map Z [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WLt = 0.5e-9;           %% Wetting Layer thickness [m]
Rx = 12E-9;             %% radius in the the x-direction of the ellipse [m]
Ry = 14E-9;             %% radius in the the y-direction of the ellipse [m]
Rz = 2.25E-9;           %% radius in the the z-direction of the ellipse [m]
x0=0;y0=0;z0=-WLt;      %% center position of the ellipse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=linspace(-Mx/2,Mx/2,Nx);
y=linspace(-My/2,My/2,Ny);
z=linspace(-Mz/2,Mz/2,Nz)+1e-9;

[X,Y,Z]=meshgrid(x,y,z);

idx_QD = ((X-x0)/Rx).^2 + ((Y-y0)/Ry).^2 + ((Z-z0)/Rz).^2 < 1 ;   %% elipse index
idx_QD(Z < z0) = 0;                                               %% cut the ellipse in half

idx_WL = (Z < 0) & (Z > -WLt);      %% Wetting Layer index

idx = idx_QD | idx_WL;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V0 = (idx)*0 + (1-idx)*Vb ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V0 = Fx*X + V0;        % adding the electric field Fx to the potential in the x-direction
V0 = Fy*Y + V0;        % adding the electric field Fy to the potential in the y-direction
V0 = Fz*Z + V0;        % adding the electric field Fz to the potential in the y-direction
V0 = V0-min(min(min(V0)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','Potential','position',[100 100 1600 800])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,1,'fontsize',15)
hold on;grid on;

pcolor(x*1e9,z*1e9,squeeze(V0(round(end/2),:,:))')
plot([-1 1]*Mx/2*1e9,[0 0],'w','linewidth',2)
plot(+[1 1]*Rx*1e9,[-1 1]*Mx*1e9,'w','linewidth',2)
plot(-[1 1]*Rx*1e9,[-1 1]*Mx*1e9,'w','linewidth',2)

colormap(cool)
colorbar

xlim([-1 1]*Mx/2*1e9)
ylim([-1 1]*My/2*1e9)

xlabel('x (nm)')
ylabel('z (nm)')
title(strcat('Potential (eV): Vxz @y=0nm'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,4,'fontsize',15)
hold on;grid on;

pcolor(y*1e9,z*1e9,squeeze(V0(:,round(end/2),:))')
plot([-1 1]*Mx/2*1e9,[0 0],'w','linewidth',2)
plot(+[1 1]*Ry*1e9,[-1 1]*My*1e9,'w','linewidth',2)
plot(-[1 1]*Ry*1e9,[-1 1]*My*1e9,'w','linewidth',2)

colormap(cool)
colorbar

xlim([-1 1]*Mx/2*1e9)
ylim([-1 1]*My/2*1e9)

xlabel('y (nm)')
ylabel('z (nm)')
title(strcat('Potential (eV): Vyz @x=0nm'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,2,'fontsize',15)
hold on;grid on;

idz=find(z>0);
idz=idz(1);

pcolor(x*1e9,y*1e9,squeeze(V0(:,:,idz)))

plot([-1 1]*My*1e9,[1 1]* Ry*1e9,'w','linewidth',2)
plot([-1 1]*Mx*1e9,[1 1]*-Ry*1e9,'w','linewidth',2)
plot([1 1]* Rx*1e9,[-1 1]*My*1e9,'w','linewidth',2)
plot([1 1]*-Rx*1e9,[-1 1]*My*1e9,'w','linewidth',2)

colormap(cool)
colorbar

xlim([-1 1]*My/2*1e9)
ylim([-1 1]*My/2*1e9)

xlabel('x (nm)')
ylabel('y (nm)')
title(strcat('Potential (eV): Vxy @z=',num2str(z(idz)*1e9,'%.2f'),'nm'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,5,'fontsize',15)
hold on;grid on;

idz=find(z>Rz-WLt);
idz=idz(1)-1;
pcolor(x*1e9,y*1e9,squeeze(V0(:,:,idz)))

plot([-1 1]*Mx*1e9,[1 1]* Ry*1e9,'w','linewidth',2)
plot([-1 1]*Mx*1e9,[1 1]*-Ry*1e9,'w','linewidth',2)
plot([1 1]* Rx*1e9,[-1 1]*My*1e9,'w','linewidth',2)
plot([1 1]*-Rx*1e9,[-1 1]*My*1e9,'w','linewidth',2)

colormap(cool)
colorbar

xlim([-1 1]*My/2*1e9)
ylim([-1 1]*My/2*1e9)

xlabel('x (nm)')
ylabel('y (nm)')
title(strcat('Potential (eV): Vxy @z=',num2str(z(idz)*1e9,'%.2f'),'nm'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,3,'fontsize',15)
hold on;grid on;
plot(z*1e9,squeeze(V0(round(end/2),round(end/2),:)) ,'b.-')
plot(z*1e9,squeeze(V0(round(end/2),end,:)) ,'r.-')

xlabel('z (nm)')
ylabel('Energy (eV)')
title(strcat('Potential (eV): Vz @y=0nm and \color{blue}x=0nm ; \color{red}x=',num2str(x(end)*1e9),'nm'))

