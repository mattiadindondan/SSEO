function [Dvtot,fb,dv,eps,m,p] = DvTot(DEP,TOF1,TOF2,LArr,P1,P2,P3,muS,orbitType,R)

%  DvTot.m - 
%
% PROTOTYPE:
% [Dv,rp] = DvTot(DEP,TOF1,ARR,P1,P2,P3,muS,orbitType)
%
% DESCRIPTION:
% the function compute the total Delta v of an ınterplanetray transfer wıth
% one fly-by around a generıc planet (n2) using Lambert Transfer method and
% patch conic method.
%
% INPUT:
%
% DEP                     [1x1]           Departure time                      [km]
% TOF1                    [1x1]           Time of flight                      [-]
% TOF2                    [1x1]           Time of flight                      [-]
% ARR                     [1x1]           Latest Arrival time                        [d]
% n1                      [1x1]           deparure ibody                      [rad]
% n2                      [1x1]           flyby ibody                         [rad]
% n3                      [1x1]           arrival ibody                       [rad]
% mu                      [1x1]           Gravitational parameter             [km^3/s^2]
% OUTPUT:
%
% r                       [3x1]           Position vector                      [km]
% v                       [3x1]           Velocity vector                      [km/s]
% m,p                     [2x1]           vectors for flybyTof function
% CONTRIBUTORS
%
% Monai Francesco
% Dora Campana
% Arda Varlı
% Marco Barbieri
% Versions: 2023-10-01 First version


% Define initial value
Dvtot = NaN;
    fb.vp_m = NaN;
    fb.vp_p = NaN;
   dv.lambert1 = NaN;
   dv.lambert2 = NaN;
m=NaN;
p = NaN;
eps = NaN;

if(DEP+TOF1+TOF2 < LArr)


            [kep1,~] = uplanet(DEP, P1);  % Mercury at departure time
            [r1, v1] = kep2car(kep1, muS); % position and velocity
            [kep2,~] = uplanet(DEP+TOF1, P2);  % Earth at departure + 1st transfer time
            [r2, ~] = kep2car(kep2, muS); % position and velocity
            [kep3,~] = uplanet(DEP+TOF1, P3);  % Saturn at arrival time
            [r3, v3] = kep2car(kep3, muS); % position and velocity

            TOF1_sec = TOF1*3600*24;
            [~,~,~,err,VI1,VF1,~,~] = lambertMR(r1, r2, TOF1_sec, muS,orbitType,0,0,0); % Note: time in [s], V is [1x3]
           if err ==0
             
                TOF2_sec = (TOF2)*3600*24;
                [~,~,~,err,VI2, VF2,~,~] = lambertMR(r2, r3, TOF2_sec, muS,orbitType,0,0,0); % Note: time in [s], V is [1x3]
                if err==0
                   [fb.rp, dv.flyby,fb.vp_m,fb.vp_p,eps,m,p] = flyby(VF1', VI2', P2, DEP+TOF1,R); % Note: time in mjd2000 [days]
                
                            
                   if isnan(dv.flyby)
                      Dvtot = NaN;
                      dv.lambert1 = NaN;
                      dv.lambert2 = NaN;
                   else
                      dv.lambert1 = norm(v1 - VI1'); % V at the beginning of the leg and V of the Earth
                      dv.lambert2 = norm(v3 - VF2');  % V at the end of the leg and V of Saturn
        
                      Dvtot=  dv.lambert1 +  dv.lambert2 + dv.flyby;
                   end
                end
           end
end
end