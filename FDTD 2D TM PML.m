%***********************************************************************
% FDTD 2D TM Code with PML
%***********************************************************************
clear all; close all;

%***********************************************************************
% Fundamental constants
%***********************************************************************
c_0 = 3.0e8;                        % Speed of light in free space
mu_0 = 4.0*pi*1.0e-7;               % Permeability of free space
eps_0 = 8.8542e-12;                 % Permittivity of free space

%***********************************************************************
% Grid parameters
%***********************************************************************
ie = 201;                           % # of grid cells in x-direction
je = 201;                           % # of grid cells in y-direction
ib = ie + 1;
jb = je + 1;
npmls = 10;                         % Depth of PML region in # of cells
ip = ie - npmls;
jp = je - npmls;

dx = 2.0e-3;                        % Space increment of square lattice
dt = dx/c_0/sqrt(2);                % Time step

is = (ie+1)/2;                      % Location of source
js = (je+1)/2;                      % Location of source

nmax = 2000;                        % Total number of time steps

aimp = sqrt(mu_0/eps_0);            % Wave impedance in free space

step = ceil(nmax/10);

%***********************************************************************
% Wave excitation
%***********************************************************************
freq = 3.0e9;
tau  = 1/freq;
delay = 3*tau/dt;                   % Delay of Gaussian source
N = round(tau/dt);
source = zeros(nmax,1);
ST = 5;                             % Switching time for sinusoidal source
for n=1:nmax
    if n<ST*N
        x = 1.0 - (ST*N-n)/(ST*N);
        g = 10.0*x^3 - 15.0*x^4 + 6.0*x^5;
        source(n) = g * sin(2*pi*freq*n*dt);
    else
        source(n) = sin(2*pi*freq*n*dt);
    end
%     source(n) = exp(-((n-delay)*dt/tau)^2) * sin(2*pi*freq*n*dt);
end

%***********************************************************************
% Initialise all field components
%***********************************************************************
Ezx = zeros(ib,jb); Ezy = zeros(ib,jb);
Ez = zeros(ib,jb); Ez_r = zeros(ib,jb); Ez_i = zeros(ib,jb); Ez_c = zeros(ib,jb);
Hx = zeros(ib,je); Hx_r = zeros(ib,je); Hx_i = zeros(ib,je); Hx_c = zeros(ib,je);
Hy = zeros(ie,jb); Hy_r = zeros(ie,jb); Hy_i = zeros(ie,jb); Hy_c = zeros(ie,jb);
caEzx = ones(ib,jb); cbEzx = ones(ib,jb).*(dt/eps_0/dx);
caEzy = ones(ib,jb); cbEzy = ones(ib,jb).*(dt/eps_0/dx);
daHx = ones(ib,je); dbHx = ones(ib,je).*(dt/mu_0/dx);
daHy = ones(ie,jb); dbHy = ones(ie,jb).*(dt/mu_0/dx);

%***********************************************************************
% Set up the Berenger's PML material constants
%***********************************************************************
sigmax = -3.0*eps_0*c_0*log(1.0e-5)/(2.0*dx*npmls);
rhomax = sigmax*(aimp^2);
for m=1:npmls
    sig(m) = sigmax*((m-0.5)/(npmls+0.5))^2;
    rho(m) = rhomax*(m/(npmls+0.5))^2;
end

%***********************************************************************
% Set up constants for Berenger's PML
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
% Initialise all matrices for the Berenger's PML
%***********************************************************************
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Ez Fields >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                ..... Left and Right PML Regions .....
for i=2:ie
    for j=2:npmls+1
        m = npmls+2-j;
        caEzy(i,j) = ca(m);
        cbEzy(i,j) = cb(m);
    end
    for j=jp+1:je
        m = j-jp;
        caEzy(i,j) = ca(m);
        cbEzy(i,j) = cb(m);
    end
end
%                ..... Front and Back PML Regions .....
for j=2:je
    for i=2:npmls+1
        m = npmls+2-i;
        caEzx(i,j) = ca(m);
        cbEzx(i,j) = cb(m);
    end
    for i=ip+1:ie
        m = i-ip;
        caEzx(i,j) = ca(m);
        cbEzx(i,j) = cb(m);
    end
end

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Hx Fields >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                ..... Left and Right PML Regions .....
for i=2:ie
    for j=1:npmls
        m = npmls+1-j;
        daHx(i,j) = da(m);
        dbHx(i,j) = db(m);
    end
    for j=jp+1:je
        m = j-jp;
        daHx(i,j) = da(m);
        dbHx(i,j) = db(m);
    end
end

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Hy Fields >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                ..... Front and Back PML Regions .....
for j=2:je
	for i=1:npmls
		m = npmls+1-i;
        daHy(i,j) = da(m);
        dbHy(i,j) = db(m);
    end
    for i=ip+1:ie
        m = i-ip;
        daHy(i,j) = da(m);
        dbHy(i,j) = db(m);
    end
end

%***********************************************************************
% Movie initialisation
%***********************************************************************
figure('position',[10 250 940 360]); set(gcf, 'color', 'white');
rect = get(gcf,'position'); rect(1:2) = [0 0];

%***********************************************************************
% Time-stepping loop
%***********************************************************************
for n=1:nmax

    %*******************************************************************
    % Ez field update
    %*******************************************************************
    Ezx(2:ie,2:je) = caEzx(2:ie,2:je).*Ezx(2:ie,2:je) + cbEzx(2:ie,2:je).*(Hy(2:ie,2:je) - Hy(1:ie-1,2:je));
    Ezy(2:ie,2:je) = caEzy(2:ie,2:je).*Ezy(2:ie,2:je) + cbEzy(2:ie,2:je).*(Hx(2:ie,1:je-1) - Hx(2:ie,2:je));
    Ez(2:ie,2:je) = Ezx(2:ie,2:je) + Ezy(2:ie,2:je);

    %*******************************************************************
    % Source excitation
    %*******************************************************************
    Ez(is,js) = source(n);

    %*******************************************************************
    % H fields update
    %*******************************************************************
    Hx(1:ib,1:je) = daHx(1:ib,1:je).*Hx(1:ib,1:je) + dbHx(1:ib,1:je).*(Ez(1:ib,1:je) - Ez(1:ib,2:jb));
    Hy(1:ie,1:jb) = daHy(1:ie,1:jb).*Hy(1:ie,1:jb) + dbHy(1:ie,1:jb).*(Ez(2:ib,1:jb) - Ez(1:ie,1:jb));

    %*******************************************************************
    % Calculate complex value
	%*******************************************************************
    if (n>nmax-N)
        Hx_r = Hx_r + Hx.*cos(2*pi*freq*dt*n);
        Hx_i = Hx_i - Hx.*sin(2*pi*freq*dt*n);
        Hy_r = Hy_r + Hy.*cos(2*pi*freq*dt*n);
        Hy_i = Hy_i - Hy.*sin(2*pi*freq*dt*n);
        Ez_r = Ez_r + Ez.*cos(2*pi*freq*dt*n);
        Ez_i = Ez_i - Ez.*sin(2*pi*freq*dt*n);
    end

    %*******************************************************************
	% Visualise fields
	%*******************************************************************
    frame = 10;
    if mod(n,frame)==0
        timestep = int2str(n);
        subplot(1,3,1); pcolor(Hx); shading interp; axis image; caxis([-3e-4 3e-4]); axis([0 jb 0 ie]);
        title(['Hx']); xlabel('y'); ylabel('x');
        subplot(1,3,2); pcolor(Hy); shading interp; axis image; caxis([-3e-4 3e-4]); axis([0 je 0 ib]);
        title(['Hy at time step ',timestep]); xlabel('y'); ylabel('x');
        subplot(1,3,3); pcolor(Ez); shading interp; axis image; caxis([-0.15 0.15]); axis([0 jb 0 ib]);
        title(['Ez']); xlabel('y'); ylabel('x');
        nn = n/frame;
        getframe(gcf,rect);
    end

    %*******************************************************************
    % Display progress
    %*******************************************************************
    if (mod(n,step)==0)
        c = round(clock);
        disp(['Current time step: ',num2str(n-mod(n,step)),', ',num2str(n/nmax*100),'% completed at ',num2str(c(4)),':',num2str(c(5)),':',num2str(c(6))]);
    end

    %*******************************************************************
    % End time stepping loop
    %*******************************************************************
end

%***********************************************************************
% Calculate complex value
%***********************************************************************
clear j;
Hx_c = (Hx_r + j*Hx_i).*(2/N);
Hy_c = (Hy_r + j*Hy_i).*(2/N);
Ez_c = (Ez_r + j*Ez_i).*(2/N);

%***********************************************************************
% Save data to file
%***********************************************************************
save results_2D_TM.mat dx dt nmax ie je is js Hx_c Hy_c Ez_c;


