%***********************************************************************
% 2D TE Code with PML
%***********************************************************************
clear all; close all; clc;

% eps_temp=eps_map;
% close;

% Physical Constants
c0   = 3e8;                                 % Light speed  in free space
eps0 = 8.8546e-12;                          % Permittivity of free space
mu_0 = pi*4e-7;                             % Permeability of free space

% Material Parameters
% cloak    = 0;                               % Set to 1 to activate the cloak
epsr     = 1.00;                            % Relative permittivity
mu_r     = 1.00;                            % Relative permeability
epsc     = eps0*2.6257;
sigma0   = 5.8e7;                           % Conductivity of the reflecting surface (copper=5.8e7)
eps_back = eps0*epsr;                       % Permittivity of background material
mu_back  = mu_0*mu_r;                       % Permeability of background material
aimp     = sqrt(mu_back/eps_back);          % Wave impedance of the background material
losstangent = 0;

% The location of the file that contains the skewed mesh info
% s0='\\staff-homes\homes$\themosk\MATLAB\Carpet\grids\themos1.GRD';

%***********************************************************************
% FDTD Grid parameters
%***********************************************************************
freq  = 5.0e8;                                % Operating frequency
omega = 2*pi*freq;
k0    = 2*pi*freq/c0*sqrt(epsr*mu_r);       % Wavevector inside the background material
lamda = 2*pi/k0;                            % Wavelength inside the background material
dx = (1/15)*2*pi/(2*pi*8e8/c0*sqrt(1.00*mu_r))  % Space increment of square lattice
dy = dx;

% Geometry of the simulation region
Ly1 = -20.0;                                 % Start length of simulation box in x [um]
Ly2 = +60.0;                                 % End   length of simulation box in x [um]
Lx1 = -15.0;                                 % Start length of simulation box in y [um]
Lx2 = +15.0;                                 % End   length of simulation box in y [um]
ie = ceil((Lx2-Lx1)/dx);                    % Number of cells in the x-dimension for the simulation box
je = ceil((Ly2-Ly1)/dy);                    % Number of cells in the y-dimension for the simulation box

Xa = linspace(Lx1,Lx2,ie);                  % Define x-coordinate
Ya = linspace(Ly1,Ly2,je);                  % Define y-coordinate
[xa,ya] = meshgrid(Xa',Ya');                % Define the 2D simulation grid
xa=xa'; ya=ya';
[Nxa,Nya] = size(xa);  
% Calculate the permitivity and conductivity tables for every point in the
% FDTD domain
eps   = ones(Nxa,Nya)*eps0;
mu    = ones(Nxa,Nya)*mu_0;
sigma = zeros(Nxa,Nya);
% sigma = rot90(sigma); 
 load 'eps_fullmap_dis_100by60_corrected.mat';
% if cloak==0, eps = ones(ie,je)*eps_back; end
% mu=eps/eps0*mu_0;

%if eps<eps0, disp('There are \epsilon<1 values'); end
%if  mu<mu_0, disp('There are \mu<1 values'); end
% eps(eps<eps0) = eps0;
%  mu(mu<mu_0) = mu_0;
%eps = eps/eps0; eps=round(eps*1e3)/1e3; eps=eps*eps0;
% % % %*********Load the original lens************************%
% fid=fopen('SingleLens.txt','r');
% data = textscan(fid, '%f %f  ');
% xl=data{1};
% yl=data{2};
% close;
% 
% Nxl=length(xl);
% flag=0;
% for i=1:Nxa
%     for j=1:Nya         
%         if xa(i,j)>-9.51 && xa(i,j)<9.51
% %             if xa(i,j)>-1.148 && xa(i,j)<1.148
%             flag=1;
%             delta=10;
%             for k=1:Nxl
%             temp=abs(xl(k)-xa(i,j));
%             if temp<delta
%                 delta=temp;
%                  ymax=yl(k);
%                  ymin=-yl(k)-3;
% %                ymin=yl(k);
%             end
%             end
%             if ya(i,j)<ymax && ya(i,j)>ymin
%                 eps(i,j)=epsc;
%             end
%         end
%     end
 % end
% 
% save '/homes/sjw30/MATLAB/cooketriplet/full_map/test';


ib = ie + 1;                                % end-index for magnetic fields
jb = je + 1;                                % end-index for magnetic fields

npmls = 20;                                 % Depth of PML region in # of cells
ip = ie - npmls;
jp = je - npmls;

it = 30;
jt = 30;
%***********************************************************************
% Wave excitation
%***********************************************************************

nmax = 100000;                               % Total number of time steps
n    = 1:nmax;
dt   = dx/c0/sqrt(2);                      % Time step and currant condition
source = zeros(1,nmax);


pulse = 0;                                  % Setup a temporal pulse? (0 = no pulse)

if pulse == 1                               % If there is a pulse
    st = 5e-9/dt;                           % Pulse time sigma in units of dt
    T0 = 4e-8/dt;                           % Pulse center in units of dt
    w0 = 2*pi*freq;                         % Cyclic frequency
%     %source = sin(2*pi*freq*n*dt).*exp(-(n-T0).^2./(2*st^2));
    A1 = 0.5*exp(-w0^2*(st*dt)^2/2)*sqrt(pi/2)*(st*dt);
    A2 = sqrt(-1)*exp(-sqrt(-1)*w0*T0*dt)*erfz((sqrt(-1)*(st*dt)^2*w0 + (n*dt-T0*dt))/(sqrt(2)*(st*dt)));
    A3 = sqrt(-1)*exp(+sqrt(-1)*w0*T0*dt)*erfz((sqrt(-1)*(st*dt)^2*w0 - (n*dt-T0*dt))/(sqrt(2)*(st*dt)));
    source = real(A1.*(A2+A3)); source=source/max(abs(source)); % Amplitude = 1

else                                        % Setup a plane wave with a small switch-on time
    N = round(1/dt/freq);                   % M = round(cc/freq/dx);
    ST = 10;                                % Wave switch-on time
    x = 1.0 - (ST*N-n)/(ST*N);
    A = 10.0*x.^3 - 15.0*x.^4 + 6.0*x.^5;
    
%     source(n< ST*N) = A(n<ST*N).*sin(2*pi*freq*n(n< ST*N)*dt);   % In early times switch-on
%     source(n>=ST*N) =            sin(2*pi*freq*n(n>=ST*N)*dt);   % Later keep constant

    source(n< ST*N) = A(n<ST*N).*exp(sqrt(-1)*2*pi*freq*n(n< ST*N)*dt);   % In early times switch-on
    source(n>=ST*N) =            exp(sqrt(-1)*2*pi*freq*n(n>=ST*N)*dt);   % Later keep constant

end
clear x;

% Setup the space-dependence of the sources

% phi = 2*pi-90*(pi/180);                     % The angle of incidence
% kx = k0*cos(phi);                           % Define the k-vector for the incident plane wave
% ky = k0*sin(phi);                           % Define the k-vector for the incident plane wave
% s  = 1/dx;
% sx = s/sin(phi);                            % The sigma_x of the Gaussian modulation in x (units of dx)
% sy = s/cos(phi);                            % The sigma_y of the Gaussian modulation in y (units of dy)
% 
% delta = 5;                                  % Number of cells to leave gap between the PML and the sources
% xv_x = (npmls+delta):(ip-delta);            % The locations of the x-positioned sources in x
% yv_x = (jp -delta+1);                       % The location  of the x-positioned sources in y
% xv_y = (npmls+delta-1);                     % The location  of the y-positioned sources in x
% yv_y = (npmls+delta):(jp-delta+1);          % The locations of the y-positioned sources in y


disp('Initializing FDTD matrices...');
%***********************************************************************
% Initialise matrices for field components
%***********************************************************************
caHx  =  ones(ie,jb); cbHx  =  ones(ie,jb).*(dt/mu_back/dx);
Hx    = zeros(ie,jb);
caHy  =  ones(ib,je); cbHy  =  ones(ib,je).*(dt/mu_back/dx);
Hy    = zeros(ib,je);
Ezx   = zeros(ie,je); Ezy   = zeros(ie,je);
Ez    = zeros(ie,je);
daEzx =  ones(ie,je); dbEzx =  ones(ie,je).*(dt/eps_back/dx);
daEzy =  ones(ie,je); dbEzy =  ones(ie,je).*(dt/eps_back/dx);
daEz  =  ones(ie,je); dbEz  =  ones(ie,je).*(dt/eps_back/dx);


caBx  =  ones(ie,jb); cbBx  =  ones(ie,jb).*(dt/mu_back/dx);
Bx    = zeros(ie,jb);
caBy  =  ones(ib,je); cbBy  =  ones(ib,je).*(dt/mu_back/dx);
By    = zeros(ib,je);
Dzx   = zeros(ie,je); Dzy   = zeros(ie,je);
Dz    = zeros(ie,je);
daDzx =  ones(ie,je); dbDzx =  ones(ie,je).*(dt/eps_back/dx);
daDzy =  ones(ie,je); dbDzy =  ones(ie,je).*(dt/eps_back/dx);
daDz  =  ones(ie,je); dbDz  =  ones(ie,je).*(dt/eps_back/dx);


a0x = ones(ie,jb); b0x = ones(ie,jb);
a1x = zeros(ie,jb); b1x = zeros(ie,jb);
a2x = zeros(ie,jb); b2x = zeros(ie,jb);

a0y = ones(ib,je); b0y = ones(ib,je); 
a1y = zeros(ib,je); b1y = zeros(ib,je);
a2y = zeros(ib,je); b2y = zeros(ib,je);

a0z = ones(ie,je); b0z = ones(ie,je);
a1z = zeros(ie,je); b1z = zeros(ie,je);
a2z = zeros(ie,je); b2z = zeros(ie,je);

omega_pz = zeros(ie,je); gamma_pz = zeros(ie,je); 

% Incident field initialisation
Ex_inc    = zeros(1,jb);   Hzx_inc = zeros(1,jb);
Ey_inc    = zeros(1,ib);   Hzy_inc = zeros(1,ib);
caEx_inc  =  ones(1,jb);  cbEx_inc =  ones(1,jb).*(dt/eps0/dx);
caEy_inc  =  ones(1,ib);  cbEy_inc =  ones(1,ib).*(dt/eps0/dx);
daHzx_inc =  ones(1,je); dbHzx_inc =  ones(1,je).*(dt/mu_0/dx);
daHzy_inc =  ones(1,ie); dbHzy_inc =  ones(1,ie).*(dt/mu_0/dx);



%***********************************************************************
% Setup space parameters
%***********************************************************************
for i = 1:ie
for j = 1:je
    
   if eps(i,j)>=eps0

    
            %%%% CONVENTIONAL
            daDzx(i,j) = (1-((sigma(i,j).*dt)./(2*eps(i,j)))) ./ (1+((sigma(i,j).*dt)./(2*eps(i,j))));
            dbDzx(i,j) =                 (dt ./(eps(i,j))/dx) ./ (1+((sigma(i,j).*dt)./(2*eps(i,j))));

            daDzy(i,j) = (1-((sigma(i,j).*dt)./(2*eps(i,j)))) ./ (1+((sigma(i,j).*dt)./(2*eps(i,j))));
            dbDzy(i,j) =                  (dt./(eps(i,j))/dx) ./ (1+((sigma(i,j).*dt)./(2*eps(i,j))));  
            
            
            %%%% DISPERSIVE
            
            
        
%             a0x(i, j) = 1/dt^2 + gamma_px(i,j)/(2*dt) + (omega_px(i,j).^2)/4;
%             a1x(i, j) = -2/dt^2 + (omega_px(i,j).^2)/2;
%             a2x(i, j) = 1/dt^2 - gamma_px(i,j)/(2*dt) + (omega_px(i,j).^2)/4;
%             b0x(i, j) = 1/dt^2 + gamma_px(i,j)/(2*dt);
%             b1x(i, j) = -2/dt^2;
%             b2x(i, j) = 1/dt^2 - gamma_px(i,j)/(2*dt);
%             
%             
%             
%             a0y(i, j) = 1/dt^2 + gamma_py(i,j)/(2*dt) + (omega_py(i,j).^2)/4;
%             a1y(i, j) = -2/dt^2 + (omega_py(i,j).^2)/2;
%             a2y(i, j) = 1/dt^2 - gamma_py(i,j)/(2*dt) + (omega_py(i,j).^2)/4;
%             b0y(i, j) = 1/dt^2 + gamma_py(i,j)/(2*dt);
%             b1y(i, j) = -2/dt^2;
%             b2y(i, j) = 1/dt^2 - gamma_py(i,j)/(2*dt); 

    else
        
%         gamma_pz(i,j) = 0;
            gamma_pz(i,j) = (2*losstangent*(eps(i,j)/eps0).*sin((omega*dt)/2))./((1-(eps(i,j)/eps0)).*dt*cos((omega*dt)/2));
            
            
%             omega_pz(i,j) = omega.*sqrt(1-(eps(i,j)/eps0));
            omega_pz(i,j) = sqrt((2*sin((omega*dt)/2).*(-2*((eps(i,j)/eps0)-1).*sin((omega*dt)/2)+...
              losstangent.*(eps(i,j)/eps0).*gamma_pz(i,j).*dt*cos((omega*dt)/2)))./((dt^2)*((cos((omega*dt)/2))^2)));
               
            a0z(i, j) = 1/dt^2 + gamma_pz(i,j)/(2*dt) + ((omega_pz(i,j).^2))/4;
            a1z(i, j) = -2/dt^2 + ((omega_pz(i,j).^2))/2;
            a2z(i, j) = 1/dt^2 - gamma_pz(i,j)/(2*dt) + ((omega_pz(i,j).^2))/4;
            b0z(i, j) = 1/dt^2 + gamma_pz(i,j)/(2*dt);
            b1z(i, j) = -2/dt^2;
            b2z(i, j) = 1/dt^2 - gamma_pz(i,j)/(2*dt);
   
            
    end   
end
end
 

% ***********************************************************************
% Setup Berenger's PML material constants
%***********************************************************************
sigmax = -3.0*eps_back*c0*log(1.0e-5)/(2.0*dx*npmls);
rhomax = sigmax*(aimp^2);
for m=1:npmls
    sig(m) = sigmax*((m-0.5)/(npmls+0.5))^2;
    rho(m) = rhomax*(m/(npmls+0.5))^2;
end

%***********************************************************************
% Setup Berenger's PML coefficients
%***********************************************************************
for m=1:npmls
      re  =   sig(m)*dt/eps_back;
      rm  =   rho(m)*dt/mu_back;
    ca(m) =   exp(-re);
    cb(m) = -(exp(-re)-1.0)/sig(m)/dx;
    da(m) =   exp(-rm);
    db(m) = -(exp(-rm)-1.0)/rho(m)/dx;
end

%***********************************************************************
% Initialize all of the matrices for the Berenger's PML
%***********************************************************************
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Dz Fields >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for j=1:je                  % Front & back
    for i=1:npmls
        m = npmls+1-i;
        daDzx(i,j) = ca(m);
        dbDzx(i,j) = cb(m);
    end
    for i=ip+1:ie
        m = i-ip;
        daDzx(i,j) = ca(m);
        dbDzx(i,j) = cb(m);
    end
end
for i=1:ie                  % Left & right
    for j=1:npmls
        m = npmls+1-j;
        daDzy(i,j) = ca(m);
        dbDzy(i,j) = cb(m);
    end
    for j=jp+1:je
        m = j-jp;
        daDzy(i,j) = ca(m);
        dbDzy(i,j) = cb(m);
    end
end
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Bx Fields >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for i=1:ie                  % Left & right
    for j=2:npmls+1
        m = npmls+2-j;
        caBx(i,j) = da(m);
        cbBx(i,j) = db(m);
    end
    for j=jp+1:je
        m = j-jp;
        caBx(i,j) = da(m);
        cbBx(i,j) = db(m);
    end
end
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<< By Fields >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for j=1:je                  % Front & back
     for i=2:npmls+1
		m = npmls+2-i;
        caBy(i,j) = da(m);
        cbBy(i,j) = db(m);
     end
    for i=ip+1:ie
        m = i-ip;
        caBy(i,j) = da(m);
        cbBy(i,j) = db(m);
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

%*********Test the permittivity tables*********%
my_eps=eps/eps0;
% save '/homes/sjw30/MATLAB/cooketriplet/full_map/test_eps';   
figure();
pcolor(xa*100,ya*100,my_eps); 
shading flat;
colorbar; 
% caxis([0.5 1.5]);
axis image;
% axis([Xr(1) Xr(end) Yr(1) Yr(end)]);
xlabel('(mm)');
ylabel('(mm)');
title('Relative Permittivity Map in FDTD grid');
% %--------------------------------------------------------------------------
% % *********Test the permeability tables*********%
% my_mu=mu/mu_0;
%  
% figure();
% pcolor(xa*100,ya*100,my_mu); 
% shading flat;
% colorbar; 
% % caxis([0.5 1.5]);
% axis image;
% % axis([Xr(1) Xr(end) Yr(1) Yr(end)]);
% xlabel('(mm)');
% ylabel('(mm)');
% title('Relative Permeability Map in FDTD grid');
%--------------------------------------------------------------------------
% figure();
% pcolor(xa,ya,sigma); 
% shading flat;
% colorbar; 
% axis image;
% xlabel('x-distance');
% ylabel('y-distance');
% title('Sigma Map in FDTD grid');
% save '/homes/sjw30/MATLAB/cooketriplet/full_map/sigma';

%***********************************************************************
% Movie initialisation
%***********************************************************************
%figure('position',[200 100 600 520]); set(gcf,'color','white');
%rect = get(gcf,'position'); rect(1:2) = [0 0];
figure(gcf)

%***********************************************************************
% Store last field distribution
%***********************************************************************
Hx_h1 = Hx;
Hy_h1 = Hy;
Ez_h1 = Ez;

Bx_h1 = Bx;
By_h1 = By;
Dz_h1 = Dz;

clear j;
%Efield = zeros(ie,je,nmax);
% Ez_angle = zeros(360,nmax); xx=zeros(1,180); yy=xx;
% Ez_line = zeros(ie,nmax);
% tr = sigma(ie/4,:); y_cent = length(tr(tr>0));

%***********************************************************************
% Main FDTD time-stepping loop
%***********************************************************************
disp('Starting Simulation');
threshold = 0.1;                  % FDTD error in percentage
n = 1; ne = 1; err = 100;
Ez_p  = zeros(ie,je);

while (n<ST*N) || (n<nmax) && (err(ne)>threshold)
    %*******************************************************************
    % Store last field distribution
    %*******************************************************************
    Hx_h2 = Hx_h1; Hx_h1 = Hx;
    Hy_h2 = Hy_h1; Hy_h1 = Hy;
    Ez_h2 = Ez_h1; Ez_h1 = Ez;
    
    Bx_h2 = Bx_h1; Bx_h1 = Bx;
    By_h2 = By_h1; By_h1 = By;
    Dz_h2 = Dz_h1; Dz_h1 = Dz;
  
    %*******************************************************************
    % Update incident Ex_inc (1-D FDTD)
    %*******************************************************************
    Ex_inc(2:jb) = caEx_inc(2:jb).*Ex_inc(2:jb) - cbEx_inc(2:jb).*(Hzx_inc(2:jb) - Hzx_inc(1:je));
    Ey_inc(2:ib) = caEy_inc(2:ib).*Ey_inc(2:ib) - cbEy_inc(2:ib).*(Hzy_inc(2:ib) - Hzy_inc(1:ie));
    
    %*******************************************************************
    % Dzx and Dzy (Dz) field update
    %*******************************************************************
    Dzx(2:ie,:  ) = daDzx(2:ie,:  ).*Dzx(2:ie,:  ) + dbDzx(2:ie,:  ).*(Hy(2:ie,     :) - Hy(1:ie-1,    :));
    Dzy( :, 2:je) = daDzy( :, 2:je).*Dzy( :, 2:je) + dbDzy( :, 2:je).*(Hx( :  ,1:je-1) - Hx( :    , 2:je));
    Dz = Dzx + Dzy;

    %*******************************************************************
    % Hz field update
    %*******************************************************************
    Ez = (b0z.*Dz + b1z.*Dz_h1 + b2z.*Dz_h2 - (a1z.*Ez_h1 + a2z.*Ez_h2)) ./ a0z;
    
    %*******************************************************************
    %          ..... Incident Source Excitation (1-D FDTD).....
%     %*******************************************************************
     Hzx_inc(25) = Hzx_inc(25) + 0*source(n);
     Hzy_inc(25) = Hzy_inc(25) + 1*source(n);
    
    % off axis point source
    % Ez(round(npmls+it),round(je/4)) = Ez(round(nmpls+it),round(je/4)) + source(n);
     % off axis point source at j = 200
    % % % % % %***********************source for convex lenses********************************%
    
   % 0 incidence
%  for j=npmls+jt:je-npmls-jt
%      Ez(npmls+jt,j) = Ez(npmls+jt,j) + source(n);
%  end  % 0 incidence
    
%     
% for j=1:650
%    Ez(298-round(j/5),75+j) = Ez(298-round(j/5),75+j) + source(n);
% end          % 23 degree incidence

%for j=1:650
%    Ez(298-round(j/4),75+j) = Ez(298-round(j/4),75+j) + source(n);
%end          % 14 degree incidence

% for j=1:650
%     Ez(428-round(j/3),75+j) = Ez(428-round(j/3),75+j) + source(n);
% end          % 18 degree incidence

% for j=1:500
%     Ez(481-round(j/2),75+j) = Ez(481-round(j/2),75+j) + source(n);
% end          % 27 degree incidence
    %*******************************************************************
    % Incident Hz values (1-D & 2-D interface)
    %*******************************************************************
    Ez(it+1:ie-it,   jt+1) = Ez(it+1:ie-it,   jt+1) - (dt/eps0/dx).*Hzx_inc(   jt+1);    % Left
    Ez(it+1:ie-it,je-jt+1) = Ez(it+1:ie-it,je-jt+1) + (dt/eps0/dx).*Hzx_inc(je-jt+1);    % Right 
    
    Ez(   it,jt+2:je-jt+1) = Ez(   it,jt+2:je-jt+1) - (dt/eps0/dx).*Hzy_inc(   it);%.*exp(-((jt+1:je-jt)-jc).^2./(2*sx^2));    % Bottom/Front
    Ez(ie-it,jt+2:je-jt+1) = Ez(ie-it,jt+2:je-jt+1) + (dt/eps0/dx).*Hzy_inc(ie-it);%.*exp(-((jt+1:je-jt)-jc).^2./(2*sx^2));    % Top/Back
    
    %*******************************************************************
    % Bx and By field update
    %*******************************************************************    
    Bx(:,1:je-1) = caBx(:,1:je-1).*Bx(:,1:je-1) + cbBx(:,1:je-1).*(Ez(:,1:je-1) - Ez(:,2:je));
    By(1:ie-1,:) = caBy(1:ie-1,:).*By(1:ie-1,:) + cbBy(1:ie-1,:).*(Ez(2:ie,:) - Ez(1:ie-1,:));

    %*******************************************************************
    % Hx and Hy field update
    %*******************************************************************
  
    Hx = (b0x.*Bx + b1x.*Bx_h1 + b2x.*Bx_h2 - (a1x.*Hx_h1 + a2x.*Hx_h2)) ./ a0x;
    Hy = (b0y.*By+ b1y.*By_h1  + b2y.*By_h2 - (a1y.*Hy_h1  + a2y.*Hy_h2))./ a0y;
    
    %*******************************************************************
    % Incident Hx values (1-D & 2-D interface)
    %*******************************************************************
    Hx(it+1:ie-it,   jt+1) = Hx(it+1:ie-it,   jt+1) - (dt/mu_0/dx).*Ex_inc(   jt+1);%.*exp(-((it+1:ie-it)-ic).^2./(2*sy^2))';     % Left
    Hx(it+1:ie-it,je-jt+1) = Hx(it+1:ie-it,je-jt+1) + (dt/mu_0/dx).*Ex_inc(je-jt+1);%.*exp(-((it+1:ie-it)-ic).^2./(2*sy^2))';     % Right

    Hx(it+1:ie-it,   jt+1) = Hx(it+1:ie-it,   jt+1) - (dt/mu_0/dx).*Ey_inc(it+1:ie-it)';%.*exp(-((   jt+1)-jc).^2./(2*sx^2))';     % Left
    Hx(it+1:ie-it,je-jt+1) = Hx(it+1:ie-it,je-jt+1) + (dt/mu_0/dx).*Ey_inc(it+1:ie-it)';%.*exp(-((je-jt+1)-jc).^2./(2*sx^2))';     % Right

    %*******************************************************************
    % Incident Hy values (1-D & 2-D interface)
    %*******************************************************************
    Hy(   it,jt+2:je-jt+1) = Hy(   it,jt+2:je-jt+1) + (dt/mu_0/dx).*Ey_inc(   it);%.*exp(-((jt+1:je-jt)-jc).^2./(2*sx^2));  % Bottom
    Hy(ie-it,jt+2:je-jt+1) = Hy(ie-it,jt+2:je-jt+1) - (dt/mu_0/dx).*Ey_inc(ie-it);%.*exp(-((jt+1:je-jt)-jc).^2./(2*sx^2));  % Top
        
    Hy(   it,jt+2:je-jt+1) = Hy(   it,jt+2:je-jt+1) + (dt/mu_0/dx).*Ex_inc(jt+2:je-jt+1);%.*exp(-((   it+1)-ic).^2./(2*sy^2))';  % Bottom
    Hy(ie-it,jt+2:je-jt+1) = Hy(ie-it,jt+2:je-jt+1) - (dt/mu_0/dx).*Ex_inc(jt+2:je-jt+1);%.*exp(-((ie-it+1)-ic).^2./(2*sy^2))';  % Top
    
    %*******************************************************************
    % Calculate incident Hzx_inc (1-D FDTD)
    %*******************************************************************
    Hzx_inc(1:je) = daHzx_inc(1:je).*Hzx_inc(1:je) + dbHzx_inc(1:je).*(Ex_inc(1:je)-Ex_inc(2:jb));
    Hzy_inc(1:ie) = daHzy_inc(1:ie).*Hzy_inc(1:ie) + dbHzy_inc(1:ie).*(Ey_inc(1:ie)-Ey_inc(2:ib));
   
    
    
    % Display progress bar
%     waitbar(n/nmax);
    
    
    %*******************************************************************
    % Visualise fields
    %*******************************************************************
%     % Evaluate the Field at the Circle
%     for agl=1:360
%         rad = round(ie/3);                                  % The radius (in # of cells) where the field should be calculated
%         Ix  = round(ie/2) + round(rad*cos(agl/180*pi));     % Locations of a circle
%         Iy  = round(je/2) + round(rad*sin(agl/180*pi));
% 
%         Ez_angle(agl,n) = (Ez(Ix,Iy));                   % Find the E-field at the circle points   
%         xx(agl) = xa(Ix,Iy);                                % Find the actual x-values at the circle points
%         yy(agl) = ya(Ix,Iy);                                % Find the actual y-values at the circle points
%     end
%     
%     % Evaluate the Field at the Line
% 
%     Ix  = 290;     % Locations of a line
%     Ez_line(:,n) = Ez(Ix,:);                   % Find the E-field at the circle points   
     
    
    
    if mod(n,50)==0
        pcolor(xa,ya,real(Ez)); shading interp; axis image; colorbar; 
%         caxis([-0.15 0.15]);
       % caxis([-max(max(real(Ez)))/4 max(max(real(Ez)))/4]);
%        xlabel('x [um]'); ylabel('y [um]'); 
%         title(['E_z @ n = ', num2str(n)]);
        xlabel('x [um]'); ylabel('y [um]'); title(['E_z  at timestep  ', num2str(n), ' err=',num2str(err(end))]);
        % Draw the PML boundaries
%         line(xa(npmls:ip,npmls   ), ya(npmls:ip,npmls   ));            % Bottom line
%         line(xa(npmls:ip,jp      ), ya(npmls:ip,jp      ));            % Top line
%         line(xa(npmls   ,npmls:jp), ya(npmls   ,npmls:jp));            % Left line
%         line(xa(ip      ,npmls:jp), ya(ip      ,npmls:jp));            % Right line
        
        % Draw field probe circle
%         line(xx,yy);
%         Efield(:,:,n/20) = Ez;
        getframe(gcf);
    end
    %*******************************************************************
    % Calculate errors for convergence
    %*******************************************************************
    if mod(n,N)==0
        ne = n/N; 
        if max(max(abs(Ez_p)))==0, err(ne) = 100;
        else
            err(ne) = max(max(abs(abs(Ez)-abs(Ez_p))))/max(max(abs(Ez_p))) * 100;
        end
        Ez_p = Ez;
 %        save '/homes/sjw30/My Documents/Programming/MATLAB/CookeTriplet/dispersiveFDTD/field_singlelens_0inci';
    end  
    n=n+1;    
end

disp('simulation complete');
 save 'Efield_planewave_0inci_lambdaby15';
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
