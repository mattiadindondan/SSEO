function T = Tsyn(P1,P2,muS)
% 
% PROTOTYPE: T = Tsyn(P1,P2,muS)
%
% DESCRIPTION:
%   The function determines the synodic period between two planets
% 
% P1 P2    Integer number identifying the celestial body (< 11)
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
% 
% INPUT
%
% P1                        [1x1]          1        ibody                  [-]
% P2                        [1x1]          2        ibody                  [-]
% muS                       [1x1]          Planet constant                 [-]
%
% OUTPUT 
%
% T                         [1x1]          Synodic Period                  [-]
% CONTRIBUTORS
%
% Monai Francesco
% Dora Campana
% Arda Varlı
% Marco Barbieri
% Versions: 2023-10-01 First version

[Kep1,~]=uplanet(0,P1);
T1 = 2*pi*sqrt(Kep1(1)^3/muS);

[Kep2,~]=uplanet(0,P2);
T2 = 2*pi*sqrt(Kep2(1)^3/muS);

T = T1*T2/abs(T1-T2);
end
