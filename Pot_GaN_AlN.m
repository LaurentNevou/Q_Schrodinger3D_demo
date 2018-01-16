%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vb=2;                   %% Potential barrier height[eV]
Mass = 0.23;            %% effective mass, constant over all the structure...
Mx=15e-9;               %% map X [m]
My=15e-9;               %% map Y [m]
Mz=4e-9;                %% map Z [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WLt = 0.5e-9;           %% Wetting Layer thickness [m]
QDh = 2e-9;           %% Quantum dot height [m]
QDd = 10e-9;            %% Quantum dot height [m]
E   = 10e8;             %% Internal Electrical Field (10MV/cm in GaN/AlN interface)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=linspace(-Mx/2,Mx/2,Nx);
y=linspace(-My/2,My/2,Ny);
z=linspace(-Mz,Mz,Nz)+0.5e-9;

[X,Y,Z]=meshgrid(x,y,z);


a=QDd/2;                      %% side length of the hexagon [m]
h=tan(pi/6) * a*sqrt(3)/2;    %% height of the pyramid with 30deg angle bw side face and base


% Hexagonal base
idx1=  abs(X)<a*sqrt(3)/2;
idx2=(tan(pi/6)*X+a>Y) .* (tan(pi/6)*X-a<Y) .* (-tan(pi/6)*X-a<Y) .* (-tan(pi/6)*X+a>Y);
idx_QD=idx1.*idx2;
idx_QD(Z < 0) = 0;
idx_QD(Z > QDh) = 0;

% Pyramid facet
idx_QD(Z > -2*h/a/sqrt(3)*X + h) = 0;
idx_QD(Z > +2*h/a/sqrt(3)*X + h) = 0;

idx_QD(Z > +h/a/sqrt(3)*X -h/a*Y + h) = 0;
idx_QD(Z > -h/a/sqrt(3)*X -h/a*Y + h) = 0;

idx_QD(Z > +h/a/sqrt(3)*X +h/a*Y + h) = 0;
idx_QD(Z > -h/a/sqrt(3)*X +h/a*Y + h) = 0;

% Wetting Layer index
idx_WL = (Z < 0) & (Z > -WLt);

% Potential index
idx = idx_QD | idx_WL;
V0=(idx)*0 + (1-idx)*Vb ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(x)
  for j=1:length(y)
    
    idz=find(idx(j,i,:)==1);
    z1=z(idz(1));
    z2=z(idz(end));
    
    Fbz= +E*         (z2-z1)/(z2-z(1));       %% Field in the barrier
    Fwz= -E*(z2-z(1)-(z2-z1))/(z2-z(1));     %% Field in the well
  
    V0(j,i,:) = (Fwz*Z(j,i,:)+Fwz*WLt).* idx(j,i,:) + V0(j,i,:) ;        % adding the electric field Fwell-z to the potential in the z-direction
  
    V0(j,i,:) = (Fbz*Z(j,i,:)+Fbz*WLt)   .* (Z(j,i,:)<-WLt)              +V0(j,i,:); % adding the electric field Fbarrier-z to the potential in the z-direction
    V0(j,i,:) = (Fbz*Z(j,i,:)-Fbz*z2+V0(j,i,1)-V0(j,i,:) ).* (Z(j,i,:)>0).*(1-idx(j,i,:)) +V0(j,i,:); % adding the electric field Fbarrier-z to the potential in the z-direction
  
  end
end


V0 = Fx*X+V0;        % adding the electric field Fx to the potential in the x-direction
V0 = Fy*Y+V0;        % adding the electric field Fy to the potential in the y-direction
V0 = Fz*Z+V0;        % adding the electric field Fz to the potential in the z-direction

V0=V0-min(min(min(V0)));

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
plot([-1 1]*10,[0 0],'w','linewidth',2)
plot([1 1]*5,[-1 1]*10,'w','linewidth',2)
plot([1 1]*-5,[-1 1]*10,'w','linewidth',2)

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
plot([-1 1]*10,[0 0],'w','linewidth',2)
plot([1 1]*5,[-1 1]*10,'w','linewidth',2)
plot([1 1]*-5,[-1 1]*10,'w','linewidth',2)

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

plot([-1 1]*10,[1 1]*5,'w','linewidth',2)
plot([-1 1]*10,[1 1]*-5,'w','linewidth',2)

plot([1 1]*5,[-1 1]*10,'w','linewidth',2)
plot([1 1]*-5,[-1 1]*10,'w','linewidth',2)

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

idz=find(z>QDh);
idz=idz(1)-1;
pcolor(x*1e9,y*1e9,squeeze(V0(:,:,idz)))

plot([-1 1]*10,[1 1]*5,'w','linewidth',2)
plot([-1 1]*10,[1 1]*-5,'w','linewidth',2)

plot([1 1]*5,[-1 1]*10,'w','linewidth',2)
plot([1 1]*-5,[-1 1]*10,'w','linewidth',2)

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

