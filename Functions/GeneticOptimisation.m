function [Ga,t] = GeneticOptimisation(Edep,Ldep,ToF1,ToF2,LA,muS,P1,P2,P3,R)
%
%
% PROTOTYPE:
% [Ga,t] = GeneticOptimisation(Edep,Ldep,ToF1,ToF2,LA,muS,P1,P2,P3,R)
%
%
% DESCRIPTION:
% The function find the optimun admisible value of Dvtot in a defined time windows 
% using a Genetic algorithm.
% (To find the best value, you need to calculate the genetic algorithm
% more than one times)
% 
% 
% INPUT:
% Edep                     [1x1]           Earliest Departure               [d]
% Ldep                     [1x1]           Latest Departure                 [d]
% ToF1                     [1x1]           First Time of Flight             [d]
% ToF2                     [1x1]           Second Time of Flight            [d]
% LArr                     [1x1]           Latest Arrival                   [d]
% muS                      [1x1]           Gravitational parameter          [km^3/s^2]
% P1                       [1x1]           Depature Planet                  [-]
% P2                       [1x1]           flyby Planet                     [-]
% P3                       [1x1]           Arrival Planet                   [-]
% R                        [1x1]           minimun Adimisible               [km] 
%                                          pericentre radius
% OUTPUT
% 
% Ga                       [struct]        Genetic algorithm result         [-] 
% t                        [3x1]           [Dep ToF1 ToF2]                  [d]
% 
% CONTRIBUTORS
%
% Monai Francesco
% Dora Campana
% Arda Varlı
% Marco Barbieri
% Versions: 2023-10-01 First version

tic
rng default  
lb = [Edep  0.3*ToF1  0.5*ToF2];              % Lower Boundary   
ub = [Ldep    1.5*ToF1    1.5*ToF2];          % Upper Boundary

options = optimoptions('ga','MaxGeneration',50,'PopulationSize',5000,...
    'FunctionTolerance',0.001,'PlotFcn',"gaplotbestf",'Display','off');
f = @(t) DVVV(t(1),t(2),t(3),LA,P1,P2,P3,muS,0,R);
n = 3;
t1 = zeros(n,3);
v = zeros(n,1);
for i = 1:n
[t1(i,:), v(i)] = ga(@(t) f(t),3,[],[],[],[],lb,ub,[],options);
end
clear i
[~,i] = min(v);
t = t1(i,:);
if(sum(t)<LA)
[Ga.Dvtot,Ga.fb,Ga.dv,~,~,~] = DvTot(t(1),t(2),t(3),LA,P1,P2,P3,muS,0,R);
Ga.elapsedtime = toc;

fprintf('Genetic Alghoritm \n\n');
fprintf('Elapsed time: %.4f s   \n',Ga.elapsedtime);
fprintf('dV          = %.4f     \t[km/s]\n',Ga.Dvtot);
fprintf('lambert1    = %.4f     \t[km/s]\n',Ga.dv.lambert1);
fprintf('lambert2    = %.4f     \t[km/s]\n',Ga.dv.lambert2);
fprintf('flyby       = %.5f     \t[km/s]\n',Ga.dv.flyby);
fprintf('rp          = %.4f  \t[km]\n', Ga.fb.rp);
fprintf('dt          = %.4f  \t[d]\n', sum(t));
else
    error('Outside the time window')
end
end
