function [r1v,v1,r2v,v2] = LambertArcs(DEP,TOF1,TOF2,muS,n1,n2,n3,orbitType,tvect)

[Kep1,~] = uplanet(DEP,n1);
[r1,~] = kep2car(Kep1,muS);

[Kep2,~] = uplanet(DEP+TOF1,n2);
[r2,~] = kep2car(Kep2,muS);

[Kep3,~] = uplanet(DEP+TOF1+TOF2,n3);
[r3,~] = kep2car(Kep3,muS);

[~,~,~,~,v1,~,~,~] = lambertMR(r1,r2,TOF1*24*3600,muS,orbitType,0,0,0);
[~,~,~,~,v2,~,~,~] = lambertMR(r2,r3,TOF2*24*3600,muS,orbitType,0,0,0);


% integrating first lambert transfer Mercury-Earth
t1 = tvect(tvect<=(DEP+TOF1))*24*3600;
[r1v, ~] = ode_orbit1(r1, v1', muS, t1);

% integrating second lambert transfer Earth-Saturn
t2 = tvect(tvect>=(DEP+TOF1))*24*3600;
[r2v, ~] = ode_orbit1(r2, v2', muS,t2);
end