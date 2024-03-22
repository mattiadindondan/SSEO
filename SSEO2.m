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

%%
% Departure, flyby, arrival identification number
dP.ID =  3;  % Departure
fbP.ID = 3 ; % Flyby
aP.ID = 5 ;  % Arrival
% Planetary constant
muS = astroConstants(4);
muE = astroConstants(13);

% Departure
Date.D.date = [2011, 08, 05, 0, 0, 0]; %[year, month, day, hour, minute, second]
Date.D.mjd2000 = date2mjd2000(Date.D.date);

% Fly-By
Date.FB.date = [2013, 10, 09, 0, 0, 0];
Date.FB.mjd2000 = date2mjd2000(Date.FB.date);

% Arrival
Date.A.date = [2016, 07, 05, 0, 0, 0];
Date.A.mjd2000 = date2mjd2000(Date.A.date);

%%
R = 5%astroConstants(23) + 100; 
t1 = Date.D.mjd2000;
t2 = Date.FB.mjd2000;
t3 = Date.A.mjd2000;
[Dvtot,fb,dv] = DvTot(t1,t2,t3,1e+10,dP.ID,fbP.ID,aP.ID,muS,0,R)



%%
%%
tvect = linspace(t1,t3,2000);
[r1,~,r2,~] = LambertArcs(t1,t2-t1,t3-t2,muS,dP.ID,fbP.ID,aP.ID,0,tvect);
clear rM
rM(:,1) = [r1(1:end-1,1);r2(:,1)];
rM(:,2) = [r1(1:end-1,2);r2(:,2)];
rM(:,3) = [r1(1:end-1,3);r2(:,3)];

ManouverPlot(0,dP.ID,fbP.ID,aP.ID,muS,t1,t2-t1,t3-t2,rM);

%%

[Gs,Time] = multiStepapproach(t1*[0.99,1.01],(t2-t1)*[0.99,1.01],(t3-t2)*[0.99,1.01],t3*1.1,dP.ID,fbP.ID,aP.ID,muS,10,R,1e-5,100,10,10)
[Dvtot,fb,dv,~,~] = DvTot(Time(1),Time(2),Time(3),1e+10,dP.ID,fbP.ID,aP.ID,muS,0,R)
%
%%
%%
tvect = linspace(Time(1),sum(Time),2000);
[r1,~,r2,~] = LambertArcs(Time(1),Time(2),Time(3)    ,muS,dP.ID,fbP.ID,aP.ID,0,tvect);
clear rM
rM(:,1) = [r1(1:end-1,1);r2(:,1)];
rM(:,2) = [r1(1:end-1,2);r2(:,2)];
rM(:,3) = [r1(1:end-1,3);r2(:,3)];

ManouverPlot(0,dP.ID,fbP.ID,aP.ID,muS,Time(1),Time(2),Time(3),rM);

mjd20002date(Time(1))
mjd20002date(Time(1)+Time(2))
mjd20002date(sum(Time))