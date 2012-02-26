function [phi] = PESF(x, N, PHIn, PHIp, PHI)
% Poisson's equation solving function

XP = length(x);
dx = diff(x);

%-------------------------------------------------------------------------
a=zeros(XP,XP); % potential [fi]
b=zeros(1,XP);

%-------------------------------------------------------------------------
% Boundary conditions block

J(1,1)=1;
J(XP,XP)=1;

b(1)=PHI(1);
b(XP)=PHI(XP);

%-------------------------------------------------------------------------
% Newton-Raphson equation solver
F = zeros(1, XP);
J = zeros(XP);

% First and last grid points F vector and J matrix values
F(1) = PHI(1) - b(1);
F(XP) = PHI(XP) - b(XP);
J(1,1) = 1;
J(XP,XP) = 1;

% Inner grid points F vector and J matrix values
for i=2:XP-1
    D2PHI = 2/(dx(i)+dx(i-1))*((PHI(i+1)-PHI(i))/dx(i) - (PHI(i)-PHI(i-1))/dx(i-1));
    F(i) = D2PHI - PHIn(i)*exp(PHI(i)) + PHIp(i)*exp(-PHI(i)) + N(i);
    J(i,i) = -2/(dx(i)+dx(i-1))*(1/dx(i) + 1/dx(i-1)) - PHIn(i)*exp(PHI(i)) - PHIp(i)*exp(-PHI(i));
    J(i,i-1) = 2/(dx(i)+dx(i-1))*(1/dx(i-1));
    J(i,i+1) = 2/(dx(i)+dx(i-1))*(1/dx(i));
end

phi=PHI-J^(-1)*F';