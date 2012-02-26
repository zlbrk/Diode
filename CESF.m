function [PHIn, PHIp] = CESF(x, mun, mup, phi, ULC, URC)
% Continuity equation solving function
% Transposition is needed

XP = length(x);
dx = diff(x);

%-------------------------------------------------------------------------
a=zeros(XP,XP); % [PHIn]
b=zeros(1,XP);

%-------------------------------------------------------------------------
% Boundary conditions block

a1(1,1) = 1;
a1(XP,XP) = 1;

b1(1) = exp(-ULC);
b1(XP) = exp(-URC);

a2(1,1) = 1;
a2(XP,XP) = 1;

b2(1) = exp(ULC);
b2(XP) = exp(URC);
%-------------------------------------------------------------------------
% Coefficient matrix evaluation block

for i=2:XP-1
    a1(i,i+1) = 2*mun(i)*exp(phi(i))/(dx(i) + dx(i-1))*1/dx(i);
    a1(i,i) = -2/(dx(i) + dx(i-1))*((mun(i)*exp(phi(i)))/dx(i) + (mun(i-1)*exp(phi(i-1)))/dx(i-1));
    a1(i,i-1) = 2*mun(i-1)*exp(phi(i-1))/(dx(i) + dx(i-1))*1/dx(i-1);
    b1(i) = 0;
    
    a2(i,i+1) = 2*mup(i)*exp(-phi(i))/(dx(i) + dx(i-1))*1/dx(i);
    a2(i,i) = -2/(dx(i) + dx(i-1))*((mup(i)*exp(-phi(i)))/dx(i) + (mup(i-1)*exp(-phi(i-1)))/dx(i-1));
    a2(i,i-1) = 2*mup(i-1)*exp(-phi(i-1))/(dx(i) + dx(i-1))*1/dx(i-1);
    b2(i) = 0;
end

%-------------------------------------------------------------------------
% Gauss-Seidel solver call
PHIn0 = b1/a1';
PHIp0 = b2/a2';

ERR = 1e-6;
Imax = 100;

PHIn = seidel(a1, b1', PHIn0', ERR, Imax);
PHIp = seidel(a2, b2', PHIp0', ERR, Imax);

%-------------------------------------------------------------------------
% Gauss-Seidel solver subfunction
function X = seidel(A, B, X0, ERR, Imax)
n = length(B);
err = 1;
X = X0;
ct = 0;

while err>ERR
    Xp = X;
    for i = 1:n
        X(i) = (B(i) - A(i,[1:i-1, i+1:n])*X([1:i-1, i+1:n]))/A(i,i);
    end
    if max(abs(X0))==0
        error ('Initial approximation is zero! You should try different initial approximation.');
        break
    end
    err = max(abs(X-Xp))/max(abs(X0));
    ct = ct+1;

    if ct > Imax
        error ('Too many iterations!');
        break
    end
end