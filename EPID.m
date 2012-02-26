function [PHI] = EPID(x, N, UL, UR)
% Electrostatic potential initial distribution

XP=length(x);

for i=1:XP
    PHIL = log((N(1)/2 + ((N(1)/2)^2 +1)^0.5)/1)+UL;
    PHIR = log((N(XP)/2 + ((N(XP)/2)^2 +1)^0.5)/1)+UR;
    kappa = (PHIR - PHIL)/(x(XP)-x(1));
    betta = PHIL;
    PHI(i)=kappa*(x(i)-x(1)) + betta;
end
PHI = PHI';