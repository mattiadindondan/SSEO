function [r, v] = ode_orbit1(r0, v0, mu, tspan)

% INPUT:  r0 [3x1]    : initial radius [km]
%         v0 [3x1]    : initial velocity [km/s]
%         mu [1]      : planet's gravitational parameter [km^3/s^2]
%         t0 [1]      : initial time [s]
%         tf [1]      : final time [s]
%         ntspan [1]  : number of elements of tspan
% 
% OUTPUT: r [ntspanx3]: position vectors [rx,ry,rz]
%         v [ntspanx3]: velocity vectors [vx,vy,vz]

%  solves the orbit as an ode for the two body problem


y0 = [r0; v0];



% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration
[~, Y] = ode113(@(t,y) ode_2bp(t,y,mu), tspan, y0, options );

r = [Y(:,1),Y(:,2),Y(:,3)];
v = [Y(:,4),Y(:,5),Y(:,6)];


end