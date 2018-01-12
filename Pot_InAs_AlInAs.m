Mass = 0.043;           %% effective mass, constant over all the structure...
Mx=50e-9;               %% map X [m]
My=100e-9;              %% map Y [m]
Mz=10e-9;               %% map Z [m]

x=linspace(-Mx/2,Mx/2,Nx);
y=linspace(-My/2,My/2,Ny);
z=linspace(-Mz/2.5,Mz/1.5,Nz);

[X,Y,Z]=meshgrid(x,y,z);

Rx = 8e-9;            %% radius in the the x-direction of the ellipse [m]
Rz = 2.5e-9;          %% radius in the the z-direction of the ellipse [m]
Ly = 32e-9;           %% length of the cylinder in the the y-direction
x0=0;y0=0;z0=0e-9;    %% center position of the ellipse & cylinder

WLt = 0.5e-9;         %% Wetting Layer thickness [m]


idx_QD = ( ((X-x0)/Rx).^2 + ((Z-z0)/Rz).^2 < 1 ) & ( abs(Y-y0)<Ly );  %% elipse & cylinder index
idx_QD(Z < z0) = 0;                                                   %% cut the ellipse in half
idx_WL = (Z < WLt/2) & (Z > -WLt/2) ;                                 %% Wetting Layer index

idx = idx_QD | idx_WL;




