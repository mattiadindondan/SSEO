function [rp,dV,vp_m,vp_p,eps,m,p] = flyby(V_inf_minus,V_inf_plus,ibody,mjd2000,R)

% INPUT: 
% 
% V_inf_minus       [3x1]        incoming heliocentric velocıty             [km/s]
% V_inf_plus        [3x1]        outgoing heliocentric velocıty             [km/s]
% ibody             [1x1]        ID planet
%                              1:   Mercury
%                              2:   Venus
%                              3:   Earth
%                              4:   Mars
%                              5:   Jupiter
%                              6:   Saturn
%                              7:   Uranus
%                              8:   Neptune
% mjd2000           [1x1] mj2000 of the flyby planet [d]
% 
% 
% 
% OUTPUT:

% dV                 [1x1]      Impulsive maneuvre at the pericenter       [km/s]
% rp                 [1x1]      Radius of percentre                        [km]
% eps                [1x1]      ratio between dV and                       [-]
%                               the total DV provided by
%                                fly by
% Rsoi               [1x1]      Radius of the sphere of influence          [Km]
% 
% CONTRIBUTORS
%
% Monai Francesco
% Dora Campana
% Arda Varlı
% Marco Barbieri
% Versions: 2023-10-01 First version
  vp_m = NaN;
  vp_p = NaN;
%% gravitational parameter

mu = astroConstants(ibody+10);
%% planet data

[kep,muS] = uplanet(mjd2000, ibody); % muS is the gravitational parameter of the Sun

% position and velocity of the planet in the heliocentric frame
[~,VP] = kep2car(kep,muS); % RP: RadiusPlanet, VP: VelocityPlanet

%% planet frame

% hyperbola incoming velocity
v_inf_minus = V_inf_minus - VP;

% hyperbola outgoing velocity
v_inf_plus = V_inf_plus - VP;
% DeltaV provided by the Hyperbola
DV = norm(V_inf_minus-V_inf_plus);
%% finding rp

% angle between incoming velocity of hyperbola1 and 
% outgoing velocity of hyperbola2
delta = acos( dot(v_inf_minus',v_inf_plus') / ( norm(v_inf_minus)*norm(v_inf_plus) ) );


% % fun = hyp1.delta/2 + hyp2.delta/2 - delta = 0
% fun = ; % must be equal to zero
% 
% fun = @(rp) atan(1./sqrt(((1 + rp.*norm(hyp1.v_inf_minus)^2/mu))^2-1))...
%           + atan(1./sqrt(((1 + rp.*norm(hyp2.v_inf_plus)^2/mu))^2-1))...
%           - delta; % must be equal to zero
% 
% fun = @(rp) acos(sqrt(((1 - 1/(1+rp.*norm(vi)^2/mu)^2))))...
%           + acos(sqrt(((1 - 1/(1+rp.*norm(vf)^2/mu)^2))))......
%           - delta; % must be equal to zero
% 
% Uncomment if there are issues with r0
% r = [0:100:10000];
% figure(1)
% plot(r,fun(r))


% while not(isnan(rp))
r0 = 6378; % taken looking at the plot of fun
options = optimset('Display','off');
% options = optimset('MaxIter',1e3,'TolFun',1e-5,'Display','off');
rp = fzero(@(rp) asin(1/(1 + rp.*norm(v_inf_minus)^2/mu))+ asin(1./(1 + rp.*norm(v_inf_plus)^2/mu))- delta,r0,options);
if isnan(rp)
    r0 = [0 6000000];
rp = fzero(@(rp) asin(1/(1 + rp.*norm(v_inf_minus)^2/mu))+ asin(1./(1 + rp.*norm(v_inf_plus)^2/mu))- delta,r0,options);
end
% figure(1)
% plot(step,rp,'Marker','*')
% hold on
%% finding deltaV
if(rp > R)
% semi-major axis
a1 = - mu / norm(v_inf_minus)^2; % [km]
a2 = - mu / norm(v_inf_plus)^2; % [km]

% eccentricity
e1 = rp * norm(v_inf_minus)^2 / mu + 1;
e2 = rp * norm(v_inf_plus)^2 / mu + 1;
%%
m = [a1,e1];
p = [a2,e2];  % vectors for the fb_Tof function
%%

% velocity at pericenter
vp_m = sqrt(mu/a1*(1 + e1)/(1 - e1));
vp_p = sqrt(mu/a2*(1 + e2)/(1 - e2));



% deltaV (is at pericenter)
dV = abs(vp_p-vp_m);


else
    dV = NaN;
    m = NaN;
    p = NaN;
end
% Ratio between dV of the impulsive manouver and DV
eps = dV/DV;
end


