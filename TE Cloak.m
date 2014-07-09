%***********************************************************************
% 2D TE Code with PML
%***********************************************************************
clear all; close all;

%***********************************************************************
% Fundamental constants
%***********************************************************************
c_0   = 3.0e8;                      % Speed of light in free space
mu_0  = 4.0*pi*1.0e-7;              % Permeability of free space
eps_0 = 8.8542e-12;                 % Permittivity of free space

losstangent = 0.0;                        %losses
%***********************************************************************
% Grid parameters
%***********************************************************************
ie = 401;                           % # of grid cells in x-direction
je = 400;                           % # of grid cells in y-direction

ib = ie + 1;
jb = je + 1;

npmls = 20;                         % Depth of PML region in # of cells
ip = ie - npmls;
jp = je - npmls;

xc = round(ie/2);
yc = round(je/2);

is = 25;                            % Location of plane-wave source
js = 25;

it = 30;                            % Top/bottom total field/scattered field region
jt = 30;                            % Left/right total/scattered field region

freq = 6.0e14;
omega = 2*pi*freq;
k_0 = omega/c_0;

dx = c_0/freq/400;                         % (lamda/30)Space increment of square lattice
dt = dx/c_0/sqrt(2);              % Time step

R1 = c_0/freq/8/dx;                          %cloaking dimensions
R2 = 1.43*R1;                          

nmax = 400000;                       % Total number of time steps

aimp = sqrt(mu_0/eps_0);            % Wave impedance in free space

threshold = 0.1;                  % FDTD error in percentage

% threshold = -40;                    % FDTD error in dB

step = ceil(nmax/10);

%***********************************************************************
% Wave excitation
%***********************************************************************
% freq = 2.0e9;
% omega = 2*pi*freq;
tau  = 1/freq;
delay = 3*tau/dt;
N = round(tau/dt); M = round(c_0/freq/dx);
source = zeros(1,nmax);
clear j;
ST = 300;
st = 20e3;

for n=1:nmax
    if n < ST*N
        x = 1.0 - (ST*N-n)/(ST*N);
        g = 10.0*x^3 - 15.0*x^4 + 6.0*x^5;
%         source(n) = g * sin(2*pi*freq*n*dt);
        source(n) = g * exp(j*2*pi*freq*n*dt);
    else
%         source(n) = sin(2*pi*freq*n*dt);
       source(n) = exp(j*2*pi*freq*n*dt);
    end
%    source(n) = exp(-(((n-delay)*dt)/tau)^2) .* exp(j*2*pi*freq*n*dt);
end

%***********************************************************************
% Initialise matrices for field components
%***********************************************************************
Dx = zeros(ie,jb); Dx_h1 = zeros(ie,jb); Dx_h2 = zeros(ie,jb); 
Dxsynch = zeros(ie,jb); Dx_h1synch = zeros(ie,jb); Dx_h2synch = zeros(ie,jb); 
caDx = ones(ie,jb); cbDx = ones(ie,jb).*(dt/eps_0/dx);
Ex = zeros(ie,jb); Ex_h1 = zeros(ie,jb); Ex_h2 = zeros(ie,jb); 
Ex_h1synch = zeros(ie,jb); Ex_h2synch = zeros(ie,jb); 
Dy = zeros(ib,je); Dy_h1 = zeros(ib,je); Dy_h2 = zeros(ib,je); 
Dysynch = zeros(ib,je); Dy_h1synch = zeros(ib,je); Dy_h2synch = zeros(ib,je);
caDy = ones(ib,je); cbDy = ones(ib,je).*(dt/eps_0/dx);
Ey = zeros(ib,je); Ey_h1 = zeros(ib,je); Ey_h2 = zeros(ib,je); 
Ey_h1synch = zeros(ib,je); Ey_h2synch = zeros(ib,je); 
Bz = zeros(ie,je); Bz_h1 = zeros(ie,je); Bz_h2 = zeros(ie,je);
Bzx = zeros(ie,je); Bzy = zeros(ie,je); Hz_p = zeros(ie,je);
daBzx = ones(ie,je); dbBzx = ones(ie,je).*(dt/mu_0/dx);
daBzy = ones(ie,je); dbBzy = ones(ie,je).*(dt/mu_0/dx);
daBz = ones(ie,je); dbBz = ones(ie,je).*(dt/mu_0/dx);
Hz = zeros(ie,je); Hz_h1 = zeros(ie,je); Hz_h2 = zeros(ie,je); 

% HzsourceTF = zeros(length(it:ie-it),nmax);
% HztransTF = zeros(length(it:ie-it),nmax);
% HzsourceSF = zeros(length(it:ie-it),nmax);
% HztransSF = zeros(length(it:ie-it),nmax);

% Coefficients for cloak
a0x = ones(ie,jb); b0x = ones(ie,jb); b0xy = zeros(ie,jb); 
a1x = zeros(ie,jb); b1x = zeros(ie,jb); a1xy = zeros(ie,jb); b1xy = zeros(ie,jb);
a2x = zeros(ie,jb); b2x = zeros(ie,jb); a2xy = zeros(ie,jb); b2xy = zeros(ie,jb);

a0y = ones(ib,je); b0y = ones(ib,je); b0yx = zeros(ib,je);
a1y = zeros(ib,je); b1y = zeros(ib,je); a1yx = zeros(ib,je); b1yx = zeros(ib,je);
a2y = zeros(ib,je); b2y = zeros(ib,je); a2yx = zeros(ib,je); b2yx = zeros(ib,je);

a0z = ones(ie,je); b0z = ones(ie,je);
a1z = zeros(ie,je); b1z = zeros(ie,je);
a2z = zeros(ie,je); b2z = zeros(ie,je);

% Incident field initialisation
Ex_inc    = zeros(1,jb);   Hzx_inc = zeros(1,jb);
Ey_inc    = zeros(1,ib);   Hzy_inc = zeros(1,ib);
caEx_inc  =  ones(1,jb);  cbEx_inc =  ones(1,jb).*(dt/eps_0/dx);
caEy_inc  =  ones(1,ib);  cbEy_inc =  ones(1,ib).*(dt/eps_0/dx);
daHzx_inc =  ones(1,je); dbHzx_inc =  ones(1,je).*(dt/mu_0/dx);
daHzy_inc =  ones(1,ie); dbHzy_inc =  ones(1,ie).*(dt/mu_0/dx);

%***********************************************************************
% cloaking device
%***********************************************************************
  
for i = 1 : ie
    for j = 1 : je

       r = sqrt((i - xc + 0.5).^2 + (j - yc).^2);
        if (r<=R2) && (r>=R1);  
            
              erx = ((r-R1) / (r));
            ethetax =  ((r) / (r-R1)); 

            sigma_px = ethetax*losstangent*((tan((omega*dt)/2))/(dt/2));

            gamma_px = (2*losstangent*erx*sin((omega*dt)/2))/((1-erx)*dt*cos((omega*dt)/2));

          omega_px = sqrt((2*sin((omega*dt)/2)*(-2*(erx-1)*sin((omega*dt)/2)+...
              losstangent*erx*gamma_px*dt*cos((omega*dt)/2)))/((dt^2)*((cos((omega*dt)/2))^2)));


            sinx = (j - yc)/r;
            cosx = (i - xc + 0.5)/r;
            
            ax = ((cosx.^2) + (ethetax.*(sinx.^2))) ./ (dt^2); 
            bx = ((gamma_px.*(cosx.^2)) + ((sigma_px + (ethetax*gamma_px)).*(sinx.^2))) ./ (2*dt);
            cx = (((cosx.^2).*((omega_px)^2)) + (sigma_px*gamma_px.*(sinx.^2))) ./ (4);
            wx = ((1-ethetax).*sinx.*cosx) ./ (dt^2); 
            fx = ((gamma_px-sigma_px-(ethetax*gamma_px)).*sinx.*cosx) ./ (2*dt);
            vx = ((((omega_px)^2)-(sigma_px*gamma_px)).*sinx.*cosx) ./ (4);
            kx = 1 / (dt^2);
            lx = gamma_px / (2*dt);
            tx = ((sinx.^2) + (ethetax.*(cosx.^2))) ./ (dt^2); 
            qx = ((gamma_px.*(sinx.^2)) + ((sigma_px + (ethetax*gamma_px)).*(cosx.^2))) ./ (2*dt);
            px = (((sinx.^2).*((omega_px)^2)) + (sigma_px*gamma_px.*(cosx.^2))) ./ (4);
            
            Ax1 = ax+bx+cx; Ax2 = (2.*ax)-(2.*cx); Ax3 = ax-bx+cx; 
            Bx1 = wx+fx+vx; Bx2 = (2.*wx)-(2.*vx); Bx3 = wx-fx+vx;
            Cx1 = kx+lx; Cx2 = (2.*kx); Cx3 = kx-lx;
            Dx1 = tx+qx+px; Dx2 = (2.*tx)-(2.*px); Dx3 = tx-qx+px;
            
          
            a0x(i, j) = Ax1 - (((Bx1).^2) ./ Dx1);
            a1x(i, j) = Ax2 - ((Bx1.*Bx2) ./ Dx1);
            a2x(i, j) = ((Bx1.*Bx3) ./ Dx1) - (Ax3);
          
            
            a1xy(i, j) = Bx2 - ((Bx1.*Dx2) ./ Dx1);
            a2xy(i, j) = ((Bx1.*Dx3) ./ Dx1) - (Bx3);
            
            
            b0x(i, j) = Cx1;
            b1x(i, j) = -Cx2;
            b2x(i, j) = Cx3;
            
            
            b0xy(i, j) = -(Bx1.*Cx1) ./ Dx1;
            b1xy(i, j) = (Bx1.*Cx2) ./ Dx1;
            b2xy(i, j) = -(Bx1.*Cx3) ./ Dx1;
           
            
            
        elseif (r<R1) 

            cbDx(i, j) = (dt/(3*eps_0)/dx);
                     
        end

        r = sqrt((i - xc).^2 + (j - yc + 0.5).^2);
        
        if (r<=R2) && (r>=R1)
            
            ery =  ((r-R1) / (r));
            ethetay =  ((r) / (r-R1)); 

            sigma_py = ethetay*losstangent*((tan((omega*dt)/2))/(dt/2));

            gamma_py = (2*losstangent*ery*sin((omega*dt)/2))/((1-ery)*dt*cos((omega*dt)/2));

            omega_py = sqrt((2*sin((omega*dt)/2)*(-2*(ery-1)*sin((omega*dt)/2)+...
              losstangent*ery*gamma_py*dt*cos((omega*dt)/2)))/((dt^2)*((cos((omega*dt)/2))^2)));

            siny = (j - yc + 0.5)/r;
            cosy = (i - xc)/r;
            
            ay = ((cosy.^2) + (ethetay*(siny.^2))) ./ (dt^2); 
            by = ((gamma_py*(cosy.^2)) + ((sigma_py + (ethetay*gamma_py)).*(siny.^2))) ./ (2*dt);
            cy = (((cosy.^2).*((omega_py)^2)) + (sigma_py*gamma_py.*(siny.^2))) ./ (4);
            wy = ((1-ethetay).*siny.*cosy) ./ (dt^2); 
            fy = ((gamma_py-sigma_py-(ethetay*gamma_py)).*siny.*cosy) ./ (2*dt);
            vy = ((((omega_py)^2)-(sigma_py*gamma_py)).*siny.*cosy) ./ (4);
            ky = 1 / (dt^2);
            ly = gamma_py / (2*dt);
            ty = ((siny.^2) + (ethetay*(cosy.^2))) ./ (dt^2); 
            qy = ((gamma_py*(siny.^2)) + ((sigma_py + (ethetay*gamma_py)).*(cosy.^2))) ./ (2*dt);
            py = (((siny.^2).*((omega_py)^2)) + (sigma_py*gamma_py.*(cosy.^2))) ./ (4);
                      
            
            Ay1 = ay+by+cy; Ay2 = (2.*ay)-(2.*cy); Ay3 = ay-by+cy; 
            By1 = wy+fy+vy; By2 = (2.*wy)-(2.*vy); By3 = wy-fy+vy;
            Cy1 = ky+ly; Cy2 = (2*ky); Cy3 = ky-ly;
            Dy1 = ty+qy+py; Dy2 = (2.*ty)-(2.*py); Dy3 = ty-qy+py;
            
          
            a0y(i, j) = Dy1 - (((By1).^2) ./ Ay1);
            a1y(i, j) = Dy2 - ((By1.*By2) ./ Ay1);
            a2y(i, j) = ((By1.*By3) ./ Ay1) - (Dy3);
             
            
            a1yx(i, j) = By2 - ((By1.*Ay2) ./ Ay1);
            a2yx(i, j) = ((By1.*Ay3) ./ Ay1) - (By3);
            
            
            b0y(i, j) = Cy1;
            b1y(i, j) = -Cy2;
            b2y(i, j) = Cy3;
           
            b0yx(i, j) = -(By1.*Cy1) ./ Ay1;
            b1yx(i, j) = (By1.*Cy2) ./ Ay1;
            b2yx(i, j) = -(By1.*Cy3) ./ Ay1;
            
                 
            
        elseif (r<R1) 

            cbDy(i, j) = (dt/(3*eps_0)/dx);
            
        end
        
        r = sqrt((i - xc).^2 + (j - yc).^2);
        
        if (r<=R2) && (r>=R1) 
            
         muz =  (((r-R1) / (r)) * (((R2) / (R2-R1))^2));
         
         if (muz < 1)
            
            gamma_pz = (2*losstangent*muz*sin((omega*dt)/2))/((1-muz)*dt*cos((omega*dt)/2));
            omega_pz = sqrt((2*sin((omega*dt)/2)*(-2*(muz-1)*sin((omega*dt)/2)+...
              losstangent*muz*gamma_pz*dt*cos((omega*dt)/2)))/((dt^2)*((cos((omega*dt)/2))^2)));         
         
            a0z(i, j) = 1/dt^2 + gamma_pz/(2*dt) + ((omega_pz^2))/4;
            a1z(i, j) = -2/dt^2 + ((omega_pz^2))/2;
            a2z(i, j) = 1/dt^2 - gamma_pz/(2*dt) + ((omega_pz^2))/4;
            b0z(i, j) = 1/dt^2 + gamma_pz/(2*dt);
            b1z(i, j) = -2/dt^2;
            b2z(i, j) = 1/dt^2 - gamma_pz/(2*dt);
            
         elseif (muz >= 1)
             
            sigma_pz = mu_0*muz*losstangent*((tan((omega*dt)/2))/(dt/2));
           
            daBzx(i, j) = (1-((sigma_pz*dt)/(2*mu_0*muz))) / (1+((sigma_pz*dt)/(2*mu_0*muz))); 
            dbBzx(i, j) = (dt/(mu_0*muz)/dx) / (1+((sigma_pz*dt)/(2*mu_0*muz)));
            daBzy(i, j) = (1-((sigma_pz*dt)/(2*mu_0*muz))) / (1+((sigma_pz*dt)/(2*mu_0*muz))); 
            dbBzy(i, j) = (dt/(mu_0*muz)/dx) / (1+((sigma_pz*dt)/(2*mu_0*muz)));
            daBz(i, j) = (1-((sigma_pz*dt)/(2*mu_0*muz))) / (1+((sigma_pz*dt)/(2*mu_0*muz))); 
            dbBz(i, j) = (dt/(mu_0*muz)/dx) / (1+((sigma_pz*dt)/(2*mu_0*muz)));
    
             
         end
        end
    end
end


%***********************************************************************
% Setup Berenger's PML material constants
%***********************************************************************
sigmax = -3.0*eps_0*c_0*log(1.0e-5)/(2.0*dx*npmls);
rhomax = sigmax*(aimp^2);
for m=1:npmls
    sig(m) = sigmax*((m-0.5)/(npmls+0.5))^2;
    rho(m) = rhomax*(m/(npmls+0.5))^2;
end

%***********************************************************************
% Setup Berenger's PML coefficients
%***********************************************************************
for m=1:npmls
    re = sig(m)*dt/eps_0;
    rm = rho(m)*dt/mu_0;
    ca(m) = exp(-re);
    cb(m) = -(exp(-re)-1.0)/sig(m)/dx;
    da(m) = exp(-rm);
    db(m) = -(exp(-rm)-1.0)/rho(m)/dx;
end

%***********************************************************************
% Initialize all of the matrices for the Berenger's PML
%***********************************************************************
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Hz Fields >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for j=1:je                  % Front & back
    for i=1:npmls
        m = npmls+1-i;
        daBzx(i,j) = da(m);
        dbBzx(i,j) = db(m);
    end
    for i=ip+1:ie
        m = i-ip;
        daBzx(i,j) = da(m);
        dbBzx(i,j) = db(m);
    end
end
for i=1:ie                  % Left & right
    for j=1:npmls
        m = npmls+1-j;
        daBzy(i,j) = da(m);
        dbBzy(i,j) = db(m);
    end
    for j=jp+1:je
        m = j-jp;
        daBzy(i,j) = da(m);
        dbBzy(i,j) = db(m);
    end
end
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Ex Fields >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for i=1:ie                  % Left & right
    for j=2:npmls+1
        m = npmls+2-j;
        caDx(i,j) = ca(m);
        cbDx(i,j) = cb(m);
    end
    for j=jp+1:je
        m = j-jp;
        caDx(i,j) = ca(m);
        cbDx(i,j) = cb(m);
    end
end
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Ey Fields >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for j=1:je                  % Front & back
	for i=2:npmls+1
		m = npmls+2-i;
        caDy(i,j) = ca(m);
        cbDy(i,j) = cb(m);
    end
    for i=ip+1:ie
        m = i-ip;
        caDy(i,j) = ca(m);
        cbDy(i,j) = cb(m);
    end
end


% PML regions for the 1D FDTD algorithm
%***********************************************************************

for j=1:npmls
  m = npmls+1-j;
        daHzx_inc(j) = da(m);
        dbHzx_inc(j) = db(m);
end
for j=jp+1:je
        m = j-jp;
        daHzx_inc(j) = da(m);
        dbHzx_inc(j) = db(m);
end

for j=2:npmls+1
  m = npmls+2-j;
        caEx_inc(j) = ca(m);
        cbEx_inc(j) = cb(m);
end
for j=jp+1:je
        m = j-jp;
        caEx_inc(j) = ca(m);
        cbEx_inc(j) = cb(m);
end
%-----------------------------------------------
for i=1:npmls
  m = npmls+1-i;
        daHzy_inc(i) = da(m);
        dbHzy_inc(i) = db(m);
end
for i=ip+1:ie
        m = i-ip;
        daHzy_inc(i) = da(m);
        dbHzy_inc(i) = db(m);
end

for i=2:npmls+1
  m = npmls+2-i;
        caEy_inc(i) = ca(m);
        cbEy_inc(i) = cb(m);
end
for i=ip+1:ie
        m = i-ip;
        caEy_inc(i) = ca(m);
        cbEy_inc(i) = cb(m);
end



%***********************************************************************
% Movie initialisation
%***********************************************************************
figure('position',[200 200 700 520]); set(gcf,'color','white');
rect = get(gcf,'position'); rect(1:2) = [0 0];

%***********************************************************************
% Store last field distribution
%***********************************************************************
Dx_h1 = Dx; Dy_h1 = Dy; Bz_h1 = Bz; Ex_h1 = Ex; Ey_h1 = Ey; Hz_h1 = Hz;

clear j;
%***********************************************************************
% Main FDTD time-stepping loop
%***********************************************************************
% See "A radially-dependant dispersive FDTD Method" for explaination of 
%time stepping formulae. 

% sy = 500;       % Sigma of gaussian propagating in y direction
% sx = 500;       % Sigma of gaussian propagating in x direction



% n = 1; ne = 1; err = 0;
n = 1; ne = 1; err = 100;

% while (n<nmax) && (err(ne)>threshold)
while (n<ST*N) || (n<nmax) && (err(ne)>threshold)


    %*******************************************************************
    % Store last field distribution
    %*******************************************************************
    Dx_h2 = Dx_h1; Dx_h1 = Dx; Dy_h2 = Dy_h1; Dy_h1 = Dy; Bz_h2 = Bz_h1; Bz_h1 = Bz;
    Ex_h2 = Ex_h1; Ex_h1 = Ex; Ey_h2 = Ey_h1; Ey_h1 = Ey; Hz_h2 = Hz_h1; Hz_h1 = Hz;
    
    
%     source(n) = exp(-(((n-delay)*dt)/(0.8*tau))^2) .*sin(2*pi*freq*n*dt);

    %*******************************************************************
    % Update incident Ex_inc (1-D FDTD)
    %*******************************************************************
    Ex_inc(2:jb) = caEx_inc(2:jb).*Ex_inc(2:jb) + cbEx_inc(2:jb).*(Hzx_inc(1:je)-Hzx_inc(2:jb));

    %*******************************************************************
    % Dx and Dy field update
    %*******************************************************************
    Dx( :  ,2:je) = caDx( :  ,2:je).*Dx( :  ,2:je) + cbDx( :  ,2:je).*(Hz( :   ,2:je) - Hz( :  ,1:je-1));
    Dy(2:ie, :  ) = caDy(2:ie, :  ).*Dy(2:ie, :  ) + cbDy(2:ie, :  ).*(Hz(1:ie-1,:  ) - Hz(2:ie, :    ));

    %*******************************************************************
    % space synchronized D and E
    %******************************************************************* 
    
    Dysynch(1:ie, 2:je) = 0.25*(Dy(1:ie, 2:je)+Dy(2:ib, 2:je)+Dy(1:ie, 1:je-1)+Dy(2:ib, 1:je-1));
    Dy_h1synch(1:ie, 2:je) = 0.25*(Dy_h1(1:ie, 2:je)+Dy_h1(2:ib, 2:je)+Dy_h1(1:ie, 1:je-1)+Dy_h1(2:ib, 1:je-1));
    Dy_h2synch(1:ie, 2:je) = 0.25*(Dy_h2(1:ie, 2:je)+Dy_h2(2:ib, 2:je)+Dy_h2(1:ie, 1:je-1)+Dy_h2(2:ib, 1:je-1));
    Ey_h1synch(1:ie, 2:je) = 0.25*(Ey_h1(1:ie, 2:je)+Ey_h1(2:ib, 2:je)+Ey_h1(1:ie, 1:je-1)+Ey_h1(2:ib, 1:je-1));
    Ey_h2synch(1:ie, 2:je) = 0.25*(Ey_h2(1:ie, 2:je)+Ey_h2(2:ib, 2:je)+Ey_h2(1:ie, 1:je-1)+Ey_h2(2:ib, 1:je-1));
    
    Dxsynch(2:ie, 2:je) = 0.25*(Dx(2:ie, 2:je)+Dx(2:ie, 3:jb)+Dx(1:ie-1, 2:je)+Dx(1:ie-1, 3:jb));
    Dx_h1synch(2:ie, 2:je) = 0.25*(Dx_h1(2:ie, 2:je)+Dx_h1(2:ie, 3:jb)+Dx_h1(1:ie-1, 2:je)+Dx_h1(1:ie-1, 3:jb));
    Dx_h2synch(2:ie, 2:je) = 0.25*(Dx_h2(2:ie, 2:je)+Dx_h2(2:ie, 3:jb)+Dx_h2(1:ie-1, 2:je)+Dx_h2(1:ie-1, 3:jb));
    Ex_h1synch(2:ie, 2:je) = 0.25*(Ex_h1(2:ie, 2:je)+Ex_h1(2:ie, 3:jb)+Ex_h1(1:ie-1, 2:je)+Ex_h1(1:ie-1, 3:jb));
    Ex_h2synch(2:ie, 2:je) = 0.25*(Ex_h2(2:ie, 2:je)+Ex_h2(2:ie, 3:jb)+Ex_h2(1:ie-1, 2:je)+Ex_h2(1:ie-1, 3:jb));

    %*******************************************************************
    % Ex and Ey field update
    %*******************************************************************
    
    Ex(1:ie, 2:je) = ((b0x(1:ie,2:je).*Dx(1:ie,2:je) + b0xy(1:ie,2:je).*Dysynch(1:ie,2:je) + ...
                      b1x(1:ie,2:je).*Dx_h1(1:ie,2:je) + b1xy(1:ie,2:je).*Dy_h1synch(1:ie,2:je) + ...
                      a1x(1:ie,2:je).*Ex_h1(1:ie,2:je) + a1xy(1:ie,2:je).*Ey_h1synch(1:ie,2:je) + ...
                      b2x(1:ie,2:je).*Dx_h2(1:ie,2:je) + b2xy(1:ie,2:je).*Dy_h2synch(1:ie,2:je) + ...
                      a2x(1:ie,2:je).*Ex_h2(1:ie,2:je) + a2xy(1:ie,2:je).*Ey_h2synch(1:ie,2:je)) ...
                      ./ (a0x(1:ie,2:je)));
    Ey(2:ie, 1:je) = ((b0y(2:ie, 1:je).*Dy(2:ie, 1:je) + b0yx(2:ie, 1:je).*Dxsynch(2:ie, 1:je) + ...
                      b1y(2:ie, 1:je).*Dy_h1(2:ie, 1:je) + b1yx(2:ie, 1:je).*Dx_h1synch(2:ie, 1:je) + ...
                      a1y(2:ie, 1:je).*Ey_h1(2:ie, 1:je) + a1yx(2:ie, 1:je).*Ex_h1synch(2:ie, 1:je) + ...
                      b2y(2:ie, 1:je).*Dy_h2(2:ie, 1:je) + b2yx(2:ie, 1:je).*Dx_h2synch(2:ie, 1:je) + ...
                      a2y(2:ie, 1:je).*Ey_h2(2:ie, 1:je) + a2yx(2:ie, 1:je).*Ex_h2synch(2:ie, 1:je)) ...
                      ./ (a0y(2:ie, 1:je)));
                
    %*******************************************************************
    % Incident Ex values (1-D & 2-D interface)
    %*******************************************************************
    Ex(it+1:ie-it,   jt+1) = Ex(it+1:ie-it,   jt+1) + (dt/eps_0/dx).*Hzx_inc(   jt+1);%.*exp(-((it+1:ie-it)-ic).^2./(2*sy^2))';     % Left
    Ex(it+1:ie-it,je-jt+1) = Ex(it+1:ie-it,je-jt+1) - (dt/eps_0/dx).*Hzx_inc(je-jt+1);%.*exp(-((it+1:ie-it)-ic).^2./(2*sy^2))';     % Right

    %*******************************************************************
    % Incident Ey values (1-D & 2-D interface)
    %*******************************************************************
    Ey(   it+1,jt+1:je-jt) = Ey(   it+1,jt+1:je-jt) - (dt/eps_0/dx).*Hzx_inc(jt+1:je-jt);%.*exp(-((   it+1)-ic).^2./(2*sy^2))';  % Bottom
    Ey(ie-it+1,jt+1:je-jt) = Ey(ie-it+1,jt+1:je-jt) + (dt/eps_0/dx).*Hzx_inc(jt+1:je-jt);%.*exp(-((ie-it+1)-ic).^2./(2*sy^2))';  % Top
    
    %*******************************************************************
    % Calculate incident Hzx_inc (1-D FDTD)
    %*******************************************************************
    Hzx_inc(1:je) = daHzx_inc(1:je).*Hzx_inc(1:je) + dbHzx_inc(1:je).*(Ex_inc(1:je)-Ex_inc(2:jb));

    %*******************************************************************
    %          ..... Incident Source Excitation (1-D FDTD).....
    %*******************************************************************
    Hzx_inc(js) = Hzx_inc(js) + 1*source(n);  

    %*******************************************************************
    % Bzx and Bzy (Bz) field update
    %*******************************************************************
    Bzx(1:ie,:  ) = daBzx(1:ie,:  ).*Bzx(1:ie,:  ) + dbBzx(1:ie,:  ).*(Ey(1:ie,:  ) - Ey(2:ib,  :));
    Bzy( :, 1:je) = daBzy( :, 1:je).*Bzy( :, 1:je) + dbBzy( :, 1:je).*(Ex( :, 2:jb) - Ex( :, 1:je));
    Bz = Bzx + Bzy;

    %*******************************************************************
    % Hz field update
    %*******************************************************************
    Hz = (b0z.*Bz + b1z.*Bz_h1 + b2z.*Bz_h2 - (a1z.*Hz_h1 + a2z.*Hz_h2)) ./ a0z;

    %*******************************************************************
    % Incident Hz values (1-D & 2-D interface)
    %*******************************************************************
    Hz(it+1:ie-it,   jt+1) = Hz(it+1:ie-it,   jt+1) - (dt/mu_0/dx).*Ex_inc(   jt+1);    % Left
    Hz(it+1:ie-it,je-jt+1) = Hz(it+1:ie-it,je-jt+1) + (dt/mu_0/dx).*Ex_inc(je-jt+1);    % Right  
    
    
%     HzsourceTF(1:length(it:ie-it),n)=Hz(it:ie-it,jt+1);
%     HztransTF(1:length(it:ie-it),n)=Hz(it:ie-it,je-jt-1);
%     HzsourceSF(1:length(it:ie-it),n)=Hz(it:ie-it,js);
%     HztransSF(1:length(it:ie-it),n)=Hz(it:ie-it,je-js);
    
    
    %*******************************************************************
    % Visualise fields
    %*******************************************************************

    if mod(n,10)==0
        timestep = num2str(n);
        pcolor((1:je)./M,(1:ie)./M,real(Hz)); shading interp; axis image; colorbar; axis([1/M je/M 1/M ie/M]);
%         pcolor((1:je)./M,(1:ie)./M,real(Hz)/max(max(real(Hz)))); shading interp; axis image; colorbar; axis([1/M je/M 1/M ie/M]);
        line((R1*sin((0:360)*pi/180)+yc)./M,(R1*cos((0:360)*pi/180)+xc)./M,'color','k');
        line((R2*sin((0:360)*pi/180)+yc)./M,(R2*cos((0:360)*pi/180)+xc)./M,'color','k');
%         xlabel('y/\lambda'); ylabel('x/\lambda'); title(['amplitude of H_z']);
        xlabel('y/\lambda'); ylabel('x/\lambda'); title(['H_z  at timestep  ', timestep]);
        getframe(gcf,rect);
    end

%     
       
    %*******************************************************************
    % Calculate errors for convergence pulse
    %*******************************************************************
    if mod(n,N)==0
        ne = n/N;
        if (20*log10(sum(sum(abs(Hz(it:ie-it,jt:je-jt)))))) > (20*log10(sum(sum(abs(Hz_p(it:ie-it,jt:je-jt))))))
         err(ne) = 0;
        else
         err(ne) = (20*log10(sum(sum(abs(Hz(it:ie-it,jt:je-jt)))))) - (20*log10(sum(sum(abs(Hz_p(it:ie-it,jt:je-jt))))));
        end
        Hz_p = Hz;
    end

      
    %*******************************************************************
    % Calculate errors for convergence
    %*******************************************************************
    if mod(n,N)==0
        ne = n/N;
        if max(max(abs(Hz_p(npmls+10:ip-10,npmls+10:jp-10))))==0
            err(ne) = 100;
        else
            err(ne) = max(max(abs(abs(Hz(npmls+10:ip-10,npmls+10:jp-10))-abs(Hz_p(npmls+10:ip-10,npmls+10:jp-10))))) / max(max(abs(Hz_p(npmls+10:ip-10,npmls+10:jp-10)))) * 100;
        end
        Hz_p = Hz;
    end

    %*******************************************************************
    % Display progress
    %*******************************************************************
    if mod(n,step)==0
        %save results_2D.mat dx dt freq nmax ie je is js err N M Ex Ey Hz;
        c = round(clock);
        disp(['Current time step: ',int2str(n-mod(n,step)),', ',int2str(n/nmax*100),'% completed at ',num2str(c(4)),':',num2str(c(5)),':',num2str(c(6))]);
    end
    
    
n = n + 1;
end

%***********************************************************************
% Save final results
%***********************************************************************
save results_2DTEidealcloakdielcore3.mat dx dt freq n nmax ie je xc yc R1 R2 is js err N M Ex Ey Hz;



