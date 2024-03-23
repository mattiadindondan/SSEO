function [VTOT, vvect] = DDV(r)

% Departure, flyby, arrival identification number
DID =  3;  % Departure
FBID = 3 ; % Flyby
AID = 5 ;  % Arrival
% Planetary constant
muS = astroConstants(4);
muE = astroConstants(13);

% Departure
Ddate = [2011, 08, 05, 15, 0, 0];
Dmjd2000 = date2mjd2000(Ddate);
[Dkep,~] = uplanet(Dmjd2000, DID); 
[Dr, v1] = kep2car(Dkep, muS); % position and velocity
% Fly-By
FBdate = [2013, 10, 09, 19, 20, 0];
FBmjd2000 = date2mjd2000(FBdate);
[FBkep,~] = uplanet(FBmjd2000, FBID); 
[FBr, v2] = kep2car(FBkep, muS); % position and velocity
% Arrival
Adate = [2016, 07, 05, 02, 30, 0];
Amjd2000 = date2mjd2000(Adate);
[Akep,~] = uplanet(Amjd2000, AID); 
[Ar, v3] = kep2car(Akep, muS); % position and velocity


% DDM
DSMdate = [2012, 09, 7, 0, 0, 0];
DSMmjd2000 = date2mjd2000(DSMdate);




R = 5;
TOF1 = DSMmjd2000-Dmjd2000;
[~,~,~,~,vl1,vl2,~,~] = lambertMR(Dr,r',TOF1*3600*24,muS,0,0,0,0);
V1 = norm(vl1'-v1);

TOF2 = FBmjd2000-DSMmjd2000;
[~,~,~,~,vl3,vl4,~,~] = lambertMR(r',FBr,TOF2*3600*24,muS,0,0,0,0);
V2 = norm(vl3-vl2);

TOF3 = Amjd2000-FBmjd2000;
[~,~,~,~,vl5,vl6,~,~] = lambertMR(FBr,Ar,TOF3*3600*24,muS,0,0,0,0);
V3 = norm(vl6'-v3);

[rp,Vfb,~,~,eps,~,~] = flyby(vl4',vl5',FBID,FBmjd2000,R);

% Time = [Dmjd2000,TOF1,TOF2];
%  norm(vl4'-vl5')
vvect = [V1 V2  Vfb V3];
VTOT = sum([ V2 V3 Vfb]);
end
