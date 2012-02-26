%------------------------------------------------------
% Diod junction current-voltage curve evaluation script

%------------------------------------------------------
% Initialization block
close all;
clear all;
clc;

T = 300; % temperature, K
kB = 1.38*1e-23; % Bolzman constant, J/K

epp=11.7; % silicon dielectric permittivity (constant), nil
ep0=8.85e-12; % void dielectric permittivity, F/m
el=1.6e-19; % electron charge, C

phit=kB*T/el; % temperature voltage, V
ni=2e16; % intrinsic carrier concentration, m^-3

mup=450e-4; % mobility of holes, m^2/V/s
mun=1400e-4; % mobility of electrons, m^2/V/s

%------------------------------------------------------
% Norming coefficients
L0=(epp*ep0*phit/el/ni)^0.5; % space, m
tau0=1e-6; % time, s
D0=L0^2/tau0; % diffusion coefficient, m^2/s
mu0=D0/phit; % mobility m^2/V/s

%------------------------------------------------------
% Gauge coefficients and assumptions
k=0.007;
A=75;
ERR=1e-6;

%------------------------------------------------------
% Space and voltage domains block
UL = 0.0; % current-voltage curve domain left boundary, V
UR = 1.3; % current-voltage curve domain right boundary, V
UP = 10; % number of current-voltage curve evaluation points, nil

XL = 0; % left space domain boundary, m
XR = 1e-6; % right space domain boundary, m
XP = 100; % number of mesh points, nil

NA = 1e21; % acceptor density, m^-3
ND = 1e23; % donor density, m^-3

%------------------------------------------------------
% Discretization block
X = XL:(XR-XL)/(XP-1):XR; % space domain mesh generation

U = UL:(UR-UL)/(UP-1):UR; % voltage domain mesh generation

%------------------------------------------------------
% Normalization block
mun = mun./mu0;

mup = mup./mu0;

x = X./L0;

U = U./phit;

UR = UR/phit;
UL = UL/phit;

XL = XL/L0;
XR = XR/L0;

NA = NA/ni;
ND = ND/ni;

%------------------------------------------------------
% Mesh functions generation block
Nd = ones(size(x)).*ND;
Na = ones(size(x)).*NA;

mun = ones(size(x)).*mun;
mup = ones(size(x)).*mup;

for i = 1:length(x)
    N(i) = Nd(i)*erfc(x(i)/k) - Na(i); % effective impurity concentration mesh function generation
end

%------------------------------------------------------
% Main loop initialization block

for i=1:1 % voltage domain cycle
    ULC = U(1);    
    URC = U(9);
    
    ULC*phit
    URC*phit
    pause

    [PHI] = EPID(x, N, ULC, URC); % initial electrostatic potential distribution

    [phi] = zeros(XP,1); % current electrostatic potential distribution

    err1 = 1; % initial error value
       
    while err1>ERR
        
        if (ULC == 0 & URC == 0)
            PHIn = ones(XP, 1);
            PHIp = ones(XP, 1);
        else
            phi = PHI; % test 2
            [PHIn, PHIp] = CESF(x, mun, mup, phi, ULC, URC); % concentration
        end
                     
        [phi] = PESF(x, N, PHIn, PHIp, PHI);
        
        err1=max(abs(PHI-phi))/max(abs(PHI))
        
        pause(0.1)
                
        % New electrostatic potential approximation
        PHI=phi;
        
        phin = log(PHIn);
        phip = log(PHIp);
        n = exp(phin+PHI);
        p = exp(phip-PHI);
        
        clc
    end

end



%------------------------------------------------------
% Visualization block

figure
semilogy(x.*L0*1e6, abs(N.*ni*1e-6)); % effective impurity concentration
xlabel('x, mkm');
ylabel('N, cm^{-3}');
grid on;

figure
plot(x.*L0*1e6, phi.*phit); % electrostatic potential
xlabel('x, mkm');
ylabel('\phi, V');
grid on;

figure
plot(x.*L0*1e6, n.*ni*1e-6); % electron density
xlabel('x, mkm');
ylabel('n, cm^{-3}');
grid on;

figure
plot(x.*L0*1e6, p.*ni*1e-6); % hole density
xlabel('x, mkm');
ylabel('p, cm^{-3}');
grid on;

figure
semilogy(x.*L0*1e6, n.*ni*1e-6, x.*L0*1e6, p.*ni*1e-6); % electron-hole density
xlabel('x, mkm');
ylabel('n, cm^{-3}, p, cm^{-3}');
grid on;