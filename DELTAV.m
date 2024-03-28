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
muJ = astroConstants(15);
% Departure
D.date = [2011, 08, 05, 15, 0, 0];
D.mjd2000 = date2mjd2000(D.date);
[D.kep,~] = uplanet(D.mjd2000, D.ID); 
[D.r, v1] = kep2car(D.kep, muS); % position and velocity

% Deep SPace Manouver
DSM.date = [2012, 09, 01, 0, 0, 0];
DSM.mjd2000 = date2mjd2000(DSM.date);

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


%% Interplanetary cruise without DSM
R = astroConstants(23)+ 200;

TOF1 = FB.mjd2000-D.mjd2000;
[~,~,~,~,vl1,vl2,~,~] = lambertMR(D.r,FB.r,TOF1*3600*24,muS,0,1,0,0);
V1 = norm(vl1'-v1);

TOF2 = A.mjd2000-FB.mjd2000;
[~,~,~,~,vl3,vl4,~,~] = lambertMR(FB.r,A.r,TOF2*3600*24,muS,0,0,0,0);
V2 = norm(vl4'-v3);

[~,Vfb,~,~,~,~,~] = flyby(vl2',vl3',FB.ID,FB.mjd2000,R);

VV = [V1 Vfb V2];
clear V1 Vfb V2 
%% Plot 
tvect = linspace(D.mjd2000,A.mjd2000,2000);
[r1,~,r2,~] = LambertArcs(D.mjd2000,TOF1,TOF2,muS,D.ID,FB.ID,A.ID,1,tvect);
clear rM
rM(:,1) = [r1(1:end-1,1);r2(:,1)];
rM(:,2) = [r1(1:end-1,2);r2(:,2)];
rM(:,3) = [r1(1:end-1,3);r2(:,3)];
% 
% ManouverPlot(0,D.ID,FB.ID,A.ID,muS,D.mjd2000,TOF1,TOF2,rM);

% searching for the Apoaxis of the first leg
a(:) = sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2);
[b,c] = max(a);
aa = r1(c,:);
clear a b c
%
% %%
% i =0;
% for k = 1
% i = i+1;
% R = 5;
% ap = aa*k; 
% TOF1 = DSM.mjd2000-D.mjd2000;
% [~,~,~,~,vl1,vl2,~,~] = lambertMR(D.r,ap',TOF1*3600*24,muS,0,0,0,0);
% V1 = norm(vl1'-v1);
% 
% TOF2 = FB.mjd2000-DSM.mjd2000;
% [~,~,~,~,vl3,vl4,~,~] = lambertMR(ap',FB.r,TOF2*3600*24,muS,0,0,0,0);
% 
% V2 = norm(vl3-vl2);
% 
% TOF3 = A.mjd2000-FB.mjd2000;
% [~,~,~,~,vl5,vl6,~,~] = lambertMR(FB.r,A.r,TOF3*3600*24,muS,0,0,0,0);
% V3 = norm(vl6'-v3);
% 
% [rp,Vfb,~,~,eps,~,~] = flyby(vl4',vl5',FB.ID,FB.mjd2000,R);
% 
% Time = [D.mjd2000,TOF1,TOF2];
%  norm(vl4'-vl5');
% dVtot(:,i) = [V1 V2 V3 Vfb];
% end
%%  Interplanetary cruise considering one DSM
clear rM
% where perform the DCM (in space) ---> minimize the total DV
AA = @(r) DDV(r);
c = fminsearch(AA,aa);

[a,b] = DDV(c);

TOF1 = DSM.mjd2000-D.mjd2000;
[~,~,~,~,vl1,vl2,~,~] = lambertMR(D.r,c',TOF1*3600*24,muS,0,0,0,0);
V1 = norm(vl1'-v1);
r = 7000;
vp = sqrt(V1^2+2*muE/r);
vvp = sqrt(muE/r*(1+0.2));
DDDDV1 = abs(vp-vvp)

TOF2 = FB.mjd2000-DSM.mjd2000;
[~,~,~,~,vl3,vl4,~,~] = lambertMR(c',FB.r,TOF2*3600*24,muS,0,0,0,0);

V2 = norm(vl3-vl2);

TOF3 = A.mjd2000-FB.mjd2000;
[~,~,~,~,vl5,vl6,~,~] = lambertMR(FB.r,A.r,TOF3*3600*24,muS,0,0,0,0);
V3 = norm(vl6'-v3);
[rp,Vfb,~,~,~,~,~] = flyby(vl4',vl5',FB.ID,FB.mjd2000,R);
%%
clear T rev ev ddv
rv = [6e+4:5e+2:10e+5];
% r = (muJ*((53.5*24*3600)/2/pi)^2)^(1/3);
vp = @(r) sqrt(V3^2+2*muJ./r);
i = 0;
ev = [0.7:0.001:0.98];
% figure()
for e = ev
%     e = 0.7
i = i+1;
vvp = @(r) sqrt(muJ./r*(1+e));
DV2 = @(r) vp(r)-vvp(r);
%

 plot(rv,DV2(rv),Color=[e,e,0])
 [dvjoi,n2] = min(DV2(rv));
 rpp = rv(n2);
hold on
 altitude = rpp-astroConstants(25);

 a = rpp/(1-e);
 T(i) = 2*pi*sqrt(a^3/muJ)/3600/24;
ddv(i) = dvjoi;
end
figure()

% plot(ev,T)
plot(ev,ddv)
hold on
plot(ev,T/53,Color='r')
plot(ev,ones(size(ev)))

[dvjoi,n]= min(ddv);
T = T(n)
a = (muJ*(T*3600*24/(2*pi))^2)^(1/3);
e = ev(n)
rpp = a*(1-e);
%% Plot 
tvect = [D.mjd2000:0.01:A.mjd2000];
%
t1 = [(D.mjd2000):0.01:(DSM.mjd2000) ]*24*3600;
[r1v, ~] = ode_orbit1(D.r, vl1', muS, t1);

t2 = [(DSM.mjd2000):0.01:(FB.mjd2000) ]*24*3600;
[r2v, ~] = ode_orbit1(c', vl3', muS, t2);

t3 = [(FB.mjd2000) :0.01:(A.mjd2000)]*24*3600;
[r3v, ~] = ode_orbit1(FB.r, vl5', muS, t3);


rM(:,1) = [r1v(1:end-1,1);r2v(:,1);r3v(1:end-1,1)];
rM(:,2) = [r1v(1:end-1,2);r2v(:,2);r3v(1:end-1,2)];
rM(:,3) = [r1v(1:end-1,3);r2v(:,3);r3v(1:end-1,3)];
% plot3(rM(:,1), rM(:,2), rM(:,3),'LineWidth',2)

%% PLot
ManouverPlot(0,D.ID,FB.ID,A.ID,muS,D.mjd2000,TOF1+TOF2,TOF3,rM);

%%
b(end) = dvjoi;
% b(1) = dvstart;


% clear rM a b c eps r1 r2 r3 tvect v1 v2 v3 VV
%clear vl1 vl2 vl3 vl4 vl5 vl6 t1 t2 t3
%clear r1v r2v r3v rp TOF1 TOF2 TOF3 R


%%%% Deltav PRM
T2 = 14*3600*24;
a2 = (muJ*(T2/(2*pi))^2)^(1/3);
e2 = 1-rpp/a2;
vPRD = sqrt(muJ*(1+e)/rpp)-sqrt(muJ*(1+e2)/rpp)
b = [b vPRD];
array2table(b,"VariableNames",["1", "dsm", "fb", "JOI","PRM"],"RowNames",["DeltaV"])


%% Science orbit change 
Tsc = 14*24*3600;

a_sc = (muJ*(Tsc/(2*pi))^2)^(1/3);
e_sc = 0.9;
dom = 0.05/180*pi;
dv_sc = 2 * sqrt(muJ / a_sc / (1 - e_sc ^ 2)) * e_sc * sin(dom / 2)
