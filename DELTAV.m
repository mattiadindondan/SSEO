clc
close all
clear all
%%
% Departure, flyby, arrival identification number
D.ID =  3;  % Departure
FB.ID = 3 ; % Flyby
A.ID = 5 ;  % Arrival
% Planetary constant
muS = astroConstants(4);
muE = astroConstants(13);

% Departure
D.date = [2011, 08, 05, 15, 0, 0];
D.mjd2000 = date2mjd2000(D.date);
[D.kep,~] = uplanet(D.mjd2000, D.ID); 
[D.r, v1] = kep2car(D.kep, muS); % position and velocity
% Fly-By
FB.date = [2013, 10, 09, 19, 20, 0];
FB.mjd2000 = date2mjd2000(FB.date);
[FB.kep,~] = uplanet(FB.mjd2000, FB.ID); 
[FB.r, v2] = kep2car(FB.kep, muS); % position and velocity
% Arrival
A.date = [2016, 07, 05, 02, 30, 0];
A.mjd2000 = date2mjd2000(A.date);
[A.kep,~] = uplanet(A.mjd2000, A.ID); 
[A.r, v3] = kep2car(A.kep, muS); % position and velocity


% DDM
DSM.date = [2012, 09, 7, 0, 0, 0];
DSM.mjd2000 = date2mjd2000(DSM.date);
% [DSM.kep,~] = uplanet(DSM.mjd2000, DSM.ID); 
% [DSM.r, v2] = kep2car(DSM.kep, muS); % position and velocity
%%
R = 5;
TOF1 = FB.mjd2000-D.mjd2000;
[~,~,~,~,vl1,vl2,~,~] = lambertMR(D.r,FB.r,TOF1*3600*24,muS,0,1,1,0);
V1 = norm(vl1'-v1);

TOF2 = A.mjd2000-FB.mjd2000;
[~,~,~,~,vl3,vl4,~,~] = lambertMR(FB.r,A.r,TOF2*3600*24,muS,0,0,0,0);
V2 = norm(vl4'-v3);

[rp,Vfb,~,~,eps,~,~] = flyby(vl2',vl3',FB.ID,FB.mjd2000,R);

Time = [D.mjd2000,TOF1,TOF2]
%%
tvect = linspace(Time(1),sum(Time),2000);
[r1,~,r2,~] = LambertArcs(Time(1),Time(2),Time(3),muS,D.ID,FB.ID,A.ID,0,tvect);
clear rM
rM(:,1) = [r1(1:end-1,1);r2(:,1)];
rM(:,2) = [r1(1:end-1,2);r2(:,2)];
rM(:,3) = [r1(1:end-1,3);r2(:,3)];

ManouverPlot(0,D.ID,FB.ID,A.ID,muS,Time(1),Time(2),Time(3),rM);
a(:) = sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2);
[b,c] = max(a);

aa = r1(c,:);




%%
i =0;
for k = 0.8:0.005:1.3
i = i+1;
R = 5;
ap = aa*k; 
TOF1 = DSM.mjd2000-D.mjd2000;
[~,~,~,~,vl1,vl2,~,~] = lambertMR(D.r,ap',TOF1*3600*24,muS,0,0,0,0);
V1 = norm(vl1'-v1);

TOF2 = FB.mjd2000-DSM.mjd2000;
[~,~,~,~,vl3,vl4,~,~] = lambertMR(ap',FB.r,TOF2*3600*24,muS,0,0,0,0);
V2 = norm(vl3-vl2);

TOF3 = A.mjd2000-FB.mjd2000;
[~,~,~,~,vl5,vl6,~,~] = lambertMR(FB.r,A.r,TOF3*3600*24,muS,0,0,0,0);
V3 = norm(vl6'-v3);

[rp,Vfb,~,~,eps,~,~] = flyby(vl4',vl5',FB.ID,FB.mjd2000,R);

Time = [D.mjd2000,TOF1,TOF2];

dVtot(:,i) = [ V2 V3 Vfb];
end
%%

[d,e] = min(sum(dVtot))
dVtot(:,e)
