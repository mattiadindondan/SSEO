function dt = Fb_tof(m,p,fb_mjd2000,Pfb)
% 
% Function to compute the duration of the flyby around Earth considering
% a finite sphere of influence
% 
% PROTOTYPE:
% dt = Fb_tof(m,p,fb_mj2000,Pfb)
%
% INPUT:
%
% m                 [2x1]          semi-major axis and eccentricity        [km,-]
%                                  of incoming hyperbola          
% p                 [2x1]          semi-major axis and eccentricity        [km,-]
%                                  of incoming hyperbola          
% fb_mj2000         [1x1]          Time of flyby in mjd2000                [d]
% Pfb               [1x1]          ID flyby planet                         [-]
% OUTPUT:
%
% dt                [1x1]          duration of the flyby                   [s]
% 
% CONTRIBUTORS
%
% Monai Francesco
% Dora Campana
% Arda Varlı
% Marco Barbieri
% Versions: 2023-10-01 First version

% Constants:
muS = astroConstants(4);                   % Sun gravitational constant
mu = astroConstants(Pfb+10);                  % Earth gravitational constant

[kep,~] = uplanet(fb_mjd2000, Pfb); % muS is the gravitational parameter of the Sun

% position and velocity of the planet in the heliocentric frame
[rP,~] = kep2car(kep,muS); % rP: RadiusPlanet

% Parameter definition
am = m(1);
em = m(2);

ap = p(1);
ep = p(2);
% Radius of Earth's SOI
r_SOI = norm(rP)*(mu/muS)^(2/5); 

% Apply hyperbola time law:
% Incoming hyperbola:
th_m = acos((am*(1-em^2)-r_SOI)/em/r_SOI);
F_m = 2*atanh(tan(th_m/2)/sqrt((1+em)/(em-1)));
dt_m = sqrt((abs(am))^3/mu) * (em*sinh(F_m)-F_m);

% Outgoing hyperbola:
th_p = acos((ap*(1-ep^2)-r_SOI)/ep/r_SOI);
F_p = 2*atanh(tan(th_p/2)/sqrt((1+ep)/(ep-1)));
dt_p = sqrt((abs(ap))^3/mu) * (ep*sinh(F_p)-F_p);

dt = dt_m + dt_p;
dt = dt/3600;
end