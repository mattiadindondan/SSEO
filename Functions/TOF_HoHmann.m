function T = TOF_HoHmann(P1,P2,muS)
% 
% PROTOTYPE: T = TOF_Hohmann(P1,P2,muS)
%
% Assumption: the orbits of the planets are collinear and circular. 
% 
% The function determine the Time Of Flight of the elliptical 
% standard Hohmann transfer 
%
% P1,P2    Integer number identifying the celestial body (< 11)
%                   1:   Mercury
%                   2:   Venus
%                   3:   Earth
%                   4:   Mars
%                   5:   Jupiter
%                   6:   Saturn
%                   7:   Uranus
%                   8:   Neptune
%                   9:   Pluto
%                   10:  Sun
% INPUT
%
% P1                        [1x1]          1        ibody                   [-]
% P2                        [1x1]          2        ibody                   [-]
% muS                       [1x1]          Planet constant                  [-]
%
% OUTPUT 
%
% T                         [1x1]          Hohmann Period                   [-]
% CONTRIBUTORS
%
% Monai Francesco
% Dora Campana
% Arda Varlı
% Marco Barbieri
% Versions: 2023-10-01 First version

[kep,~] = uplanet(0,P1);
a1 = kep(1);

[kep,~] = uplanet(0,P2);
a2 = kep(1);
clear kep
T = pi*sqrt((a1+a2)^3/(8*muS));
end