%% Main Script
clc
close all
clear
%% SET PATH
% global path
filePath = fileparts(mfilename('fullpath'));
currentPath = pwd;
if not(strcmp(filePath, currentPath))
    cd (filePath);
    currentPath = filePath;
end

addpath(genpath(currentPath));

% functions path
addpath(genpath('functions'))
% Constant
muS = astroConstants(4);
muE = astroConstants(13);
R = astroConstants(23) + 100; % Von Karman altitude;
%%

fprintf('----------------------------------------------\n');
fprintf('     INTERPLANETARY MISSION ASSIGNMENT        \n');
fprintf('----------------------------------------------\n');

%%
% Departure, flyby, arrival identification number
dP.ID = 1 ;  % Departure
fbP.ID = 3 ; % Flyby
aP.ID = 6 ;  % Arrival

% Earliest Departure
Date.ED.date = [2028, 01, 01, 0, 0, 0];
Date.ED.mjd2000 = date2mjd2000(Date.ED.date);

% Latest Departure
Date.LD.date = [2058, 01, 01, 0, 0, 0];
Date.LD.mjd2000 = date2mjd2000(Date.LD.date);

% Latest Arrival
Date.LA.date = [2058, 01, 01, 0, 0, 0];
Date.LA.mjd2000 = date2mjd2000(Date.LA.date);

%% Orbital Periods 

[dP.kep,~]=uplanet(0,dP.ID);                % Departure Planet 
dP.T = 2*pi*sqrt(dP.kep(1)^3/muS)/86400;   

[fbP.kep,~]=uplanet(0,fbP.ID);
fbP.T = 2*pi*sqrt(fbP.kep(1)^3/muS)/86400;  % flyby planet

[aP.kep,~]=uplanet(0,aP.ID);
aP.T = 2*pi*sqrt(aP.kep(1)^3/muS)/86400;    % Arrival planet

% Synodic Period

dP.Tsynfb = Tsyn(dP.ID,fbP.ID,muS)/86400;
fbP.Tsynd = dP.Tsynfb;

fbP.Tsyna = Tsyn(fbP.ID,aP.ID,muS)/86400;
aP.Tsynfb = fbP.Tsyna;

aP.Tsynd = Tsyn(aP.ID,dP.ID,muS)/86400;
dP.Tsyna = aP.Tsynd;


% Hohmann Time of flight

ToF1 = TOF_HoHmann(dP.ID,fbP.ID,muS)/86400;
ToF2 = TOF_HoHmann(fbP.ID,aP.ID,muS)/86400;

%%  multi-step approach
% time: 2521.9 s

dep = [Date.ED.mjd2000,Date.LD.mjd2000]; 
Tof1g = [0.3*ToF1 2*ToF1];
Tof2g = [0.3*ToF2 2*ToF2];

[Gs,tGs] = multiStepapproach(dep,Tof1g,Tof2g,Date.LA.mjd2000,dP.ID,fbP.ID,aP.ID,muS,8,R,1e-6,1000,5,5);
% pause(5)
%close(figure(1))

%% Genetic Algotitm Optimisation
% % time: 1249.5 s
% 
prompt= '\nHello, do you have the "Global optimisation Toolbox" ? \n yes = 1      no = 0\n';
flag = input(prompt);


if flag == 1
[Ga,tGa] = GeneticOptimisation(Date.ED.mjd2000,Date.LD.mjd2000,ToF1,ToF2,Date.LA.mjd2000,muS,dP.ID,fbP.ID,aP.ID,R);
else
    fprintf('never mind, this values will be use this values:\n')
tGa =  1e+4*[1.684519222614853   0.004114686288282   0.183986078594262];
[Ga.Dvtot,Ga.fb,Ga.dv,~,~] = DvTot(tGa(1),tGa(2),tGa(3),Date.LA.mjd2000,dP.ID,fbP.ID,aP.ID,muS,0,R);

fprintf('dep   = %.4f   [days]\n',tGa(1))
fprintf('ToF1  = %.4f  \t [days]\n',tGa(2))
fprintf('ToF2  = %.4f  \t [days]\n',tGa(3))
end


%% Best Dv date
if(Gs.Dvtot>Ga.Dvtot)
[Dep] = (tGa(1));
[Flyby] = (tGa(1)+tGa(2));
[Arr] = (tGa(1)+tGa(2)+tGa(3));
else
[Dep] = (tGs(1));
[Flyby] = (tGs(1)+tGs(2));
[Arr] = (tGs(1)+tGs(2)+tGs(3));
end

[Best.Dvtot,Best.fb,Best.dv,eps,m,p] = DvTot(Dep,Flyby-Dep,Arr-Flyby,Date.LA.mjd2000,dP.ID,fbP.ID,aP.ID,muS,0,R);

%% Pork Chop plot
n = 100;
%
DEPa = Date.ED.mjd2000:(Dep-Date.ED.mjd2000)/30:Date.LD.mjd2000;
ndep = length(DEPa);
TOF1a= 0.2*ToF1:(Flyby-Dep)/10:2*ToF1;
ntof1 =  length(TOF1a);
%
DEP = linspace(Date.ED.mjd2000,Date.LD.mjd2000,n);
TOF1 = linspace(0.2*ToF1,2*ToF1,n);
TOF2 = linspace(0.2*ToF2,2*ToF2,n);
dv_lamb1 = zeros(ndep,ntof1);
dv_lamb2 = zeros(n,n);

for i = 1 : ndep
    
    dep = DEPa(i);
    [kep,~] = uplanet(dep, dP.ID);  % Mercury at departure time
    [r1, v1] = kep2car(kep, muS); % position and velocity
   
    for j = 1:ntof1
       
        tof1 = TOF1a(j);
        [kep,~] = uplanet(dep+tof1, fbP.ID);  % Mercury at departure time
        [r2, v2] = kep2car(kep, muS); % position and velocity
        [~,~,~,~,VI1,VF1,~,~]  = lambertMR(r1,r2,tof1*86400,muS,0,0,0,0);
        
        dv_lamb1(i,j) = norm(VI1'-v1);  %+norm(VF1'-v2);
     
        
            clear VI1 VF1 r2 v2
    end
        clear r1 v1 
end

 for j = 1:n
       
        tof1 = TOF1(j);
        [kep,~] = uplanet(Dep+tof1, fbP.ID);  % Mercury at departure time
        [r2, ~] = kep2car(kep, muS); % position and velocity
for k = 1:n
        
            tof2 = TOF2(k);
            [kep,~] = uplanet(Dep+tof1+tof2, aP.ID);  % with fixed Departure time
            [r3, v3] = kep2car(kep, muS); % position and velocity
            [~,~,~,~,VI2,VF2,~,~]  = lambertMR(r2,r3,tof2*86400,muS,0,0,0,0); 
            
            dv_lamb2(j,k) = norm(VF2'-v3);% + norm(VI2'-v2);
            clear VI2 VF2 r3 v3        
end
 end
clear i j k n
clear dep tof1 tof2 kep 

%% plots setup 1st pork-chop plot

nx1 = round(linspace(1,length(DEPa),8)); % number of dates seen on X axis
ny1 = round(linspace(1,length(TOF1a),8)); % number of dates seen on Y axis
 dep_dates1 = zeros(length(nx1),6);
 
for j = nx1 % just for the axis written values
    dep_dates1(j,:) = mjd20002date(DEPa(j));
    strX1(j,:) = sprintf('%2.0f/%2.0f/%4.0f',dep_dates1(j,3:-1:1));
end

%% plots setup 2nd pork-chop plot

DEP2 = Dep+TOF1; % [days]

nx2 = round(linspace(1,length(DEP2),8)); % number of dates seen on X axis
ny2 = round(linspace(1,length(TOF2),8)); % number of dates seen on Y axis
dep_dates2 = zeros(length(nx2),6);

for j = nx2 % just for the axis written values
    dep_dates2(j,:) = mjd20002date(DEP2(j));
    strX2(j,:) = sprintf('%2.0f/%2.0f/%4.0f',dep_dates2(j,3:-1:1));
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%            Plots             %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pork chp plot of the two transfer leg
figure()
hold on
contour(DEPa',(TOF1a)',dv_lamb1',0:1:50); %'
plot(Dep,(Flyby-Dep),'Marker','hexagram','MarkerSize',9,'MarkerFaceColor','k')
xlabel('Departure Date from Mercury [dd/mm/yyyy]');
ylabel('Time Of Flight Me-E [years]');
xtickangle(45);
grid on
set(gca,'XTick',DEPa(nx1))
set(gca,'xticklabel',strX1(nx1,:))
set(gca,'YTick',TOF1a(ny1))
set(gca,'yticklabel',TOF1a(ny1)/365)
hold on
title('First Lambert leg: Mercury-->Earth Pork Chop Plot','FontSize',15)

figure()
hold on
plot(Flyby,(Arr-Flyby),'Marker','hexagram','MarkerSize',9)
contour((Dep+TOF1),TOF2,dv_lamb2',0:0.5:10);

xlabel('Departure Date from Earth [dd/mm/yyyy]','FontSize',15);
ylabel('Time Of Flight E-S [years]','FontSize',15);
xtickangle(45);
set(gca,'XTick',DEP2(nx2))
set(gca,'xticklabel',strX2(nx2,:))
set(gca,'YTick',TOF2(ny2))
set(gca,'yticklabel',TOF2(ny2)/365)
hold on
grid on
title('Second Lambert leg: Earth-->Saturn Pork Chop Plot','FontSize',15) 

clear nx1 ny1 dep_dates1 strX1 DEP2 nx2 ny2 dep_dates2 strX2 j ndep ntof1
clear TOF2 TOF1 TOF1a DEP Tof2g Tof1g
clear dv_lamb1 dv_lamb2


%%
tvect = linspace(Dep,Arr,2000);
[r1,~,r2,~] = LambertArcs(Dep,Flyby-Dep,Arr-Flyby,muS,dP.ID,fbP.ID,aP.ID,0,tvect);
clear rM
rM(:,1) = [r1(1:end-1,1);r2(:,1)];
rM(:,2) = [r1(1:end-1,2);r2(:,2)];
rM(:,3) = [r1(1:end-1,3);r2(:,3)];

ManouverPlot(1,dP.ID,fbP.ID,aP.ID,muS,Dep,Flyby-Dep,Arr-Flyby,rM);

%% Hyperbolic flyby
HyperbolaPlot(Best.fb.rp,Best.fb.vp_m,Best.fb.vp_p,3,muE);

dt = Fb_tof(m,p,Flyby,fbP.ID);
clear m p
%% Numerical result
File = fopen("Result.txt",'w');
fprintf(File,'%.4f\n',Dep);
fprintf(File,'%.4f\n',Flyby);
fprintf(File,'%.4f\n',Arr);
fprintf(File,'%.3f\n',Best.Dvtot);
fprintf(File,'%.2f\n', Best.fb.rp);
fprintf(File,'%.3f\n', dt);
fprintf(File,'%f\n', eps);
fprintf(File,'37340.199');
fclose(File);
%%
Dep = mjd20002date(Dep);
Flyby = mjd20002date(Flyby);
Arr = mjd20002date(Arr);

%%
a = fopen("Assignment.txt",'w');
fprintf(a,'----------------------------------------------\n');
fprintf(a,'     INTERPLANETARY MISSION ASSIGNMENT        \n');
fprintf(a,'----------------------------------------------\n');
fprintf(a,'----------------------------------------------\n');
fprintf(a,'Genetic Alghoritm \n\n');
if flag == 1 
    fprintf(a,'Elapsed time: %.4f [s]   \n',Ga.elapsedtime);
end
fprintf(a,'dV        = %.4f [km/s]\n',Ga.Dvtot);
fprintf(a,'lambert1  = %.4f [km/s]\n',Ga.dv.lambert1);
fprintf(a,'lambert2  = %.4f [km/s]\n',Ga.dv.lambert2);
fprintf(a,'flyby     = %.4f [km/s]\t\n',Ga.dv.flyby);
fprintf(a,'rp        = %.4f [km]\n', Ga.fb.rp);
fprintf(a,'dt        = %.4f days\n', sum(tGa));
if flag ==1
end
fprintf(a,'----------------------------------------------\n');
fprintf(a,'----------------------------------------------\n');
fprintf(a,'MultiStep Grid search \n\n');
fprintf(a,'Elapsed time: %.4f [s]   \n',Gs.elapsedtime);
fprintf(a,'dV        = %.4f [km/s]\n',Gs.Dvtot);
fprintf(a,'lambert2  = %.4f [km/s]\n',Gs.dv.lambert1);
fprintf(a,'lambert1  = %.4f [km/s]\n',Gs.dv.lambert2);
fprintf(a,'flyby     = %.4f [km/s]\t\n',Gs.dv.flyby);
fprintf(a,'rp        = %.4f [km]\n', Ga.fb.rp);
fprintf(a,'dt        = %.4f days\n', sum(tGs));
fprintf(a,'----------------------------------------------\n');
fprintf(a,'Departure Date\t%.0f/%.0f/%.0f \t%.0f:%.0f:%.0f \t\n', Dep(3),Dep(2),Dep(1),Dep(4),Dep(5),Dep(6));
fprintf(a,'Flyby Date     \t%.0f/%.0f/%.0f \t%.0f:%.0f:%.0f     \n', Flyby(3),Flyby(2),Flyby(1),Flyby(4),Flyby(5),Flyby(6));
fprintf(a,'Arrival Date\t%.0f/%.0f/%.0f \t%.0f:%.0f:%.0f    \n', Arr(3),Arr(2),Arr(1),Arr(4),Arr(5),Arr(6));
fclose(a);

