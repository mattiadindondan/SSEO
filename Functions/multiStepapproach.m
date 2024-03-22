function [Gs,Time] = multiStepapproach(dep,TOF1g,TOF2g,LA,P1,P2,P3,muS,m,R,toll,ndep,nTOF1,nTOF2)
%
%
% PROTOTYPE:[Gs,Time] = multiStepapproach(dep,TOF1g,TOF2g,LA,dP,fbP,aP,muS,m,R,toll,ndep,nTOF1,nTOF2
%
%
% DESCRIPTION:
% Grid search algorithm iterated until the error is less a defined
% tollerance.
% Every iteretion the 'time windows' are scaled and defined around the
% local minimun found at the iteration before
% 
% INPUT:
% dep                      [1x2]           [Earliest Latest] Departure      [d]
% TOF1g                    [1x2]           [min max] ToF1                   [d]
% TOF2g                    [1x2]           [min max] ToF2                   [d]
% LA                       [1x1]           Latest Arrival                   [d]
% P1                       [1x1]           Depature Planet                  [-]
% P2                       [1x1]           flyby Planet                     [-]
% P3                       [1x1]           Arrival Planet                   [-]
% m                        [1x1]           Decreasing parameter             [-]
% muS                      [1x1]           Gravitational parameter          [km^3/s^2]
% R                        [1x1]           minimun Adimisible               [km] 
%                                          pericentre radius
% toll                     [1x1]           tolerance                        [km/s]
% ndep                     [1x1]           number of dep considered         [-]
% nTOF1                    [1x1]           number of TOF1 considered        [-]
% nTOF2                    [1x1]           number of TOF2 considered        [-]
%
%
% OUTPUT
% 
% Gs                       [struct]        Grid search result               [-] 
% t                        [3x1]           [Dep ToF1 ToF2]                  [d]
% 
% CONTRIBUTORS
%
% Monai Francesco
% Dora Campana
% Arda Varlı
% Marco Barbieri
% Versions: 2023-10-01 First version

fprintf('----------------------------------------------\n');
fprintf('\nGrid Search \n\n');

tic

DepStar = dep(1);
DepEnd  = dep(2);

TOF1Star = TOF1g(1) ;
TOF1End  = TOF1g(2);

TOF2Star = TOF2g(1);
TOF2End  = TOF2g(2);

err = 1;
K = 200;
a = 0;
  DEPm = (sum(dep))/2;
  TOF1m = TOF1g;
  TOF2m = TOF2g;

dVmin = 100;
figure()
% axis([0 5 15 50])
xlabel('iteration')
ylabel('function value')
grid on
hold on
while (err > toll)
a = a+1;
DEP = linspace(DepStar, DepEnd, ndep);
TOF1 = linspace(TOF1Star,TOF1End, nTOF1);
TOF2 = linspace(TOF2Star, TOF2End, nTOF2);

for i = 1 : ndep

    for j = 1 : nTOF1
      

        for k = 1 : nTOF2 % window for the flyby
   
               [dv,~,~,~,~,~] = DvTot(DEP(i),TOF1(j),TOF2(k),LA,P1,P2,P3,muS,0,R);

               if(dv <= dVmin && not(isnan(dv)))
                  K = dVmin;
                  dVmin = dv;
                   DEPm = DEP(i);
                   TOF1m = TOF1(j);
                   TOF2m = TOF2(k);
               end
        end
    end

end
        err = (K-dVmin);
        plot(a,dVmin,'Marker','o','Color','b');
        hold on
         fprintf('Iteration N° %d \n  error: %.5f\n',a,err)



dDEP  = DepEnd - DepStar;
dTOF1 = TOF1End - TOF1Star;
dTOF2 = TOF2End - TOF2Star;

if(DEPm-dDEP/m > DepStar)
    DepStar  =  DEPm-dDEP/m;
end
if(DEPm+dDEP/m < DepEnd)
    DepEnd = DEPm+dDEP/m;
end

if(TOF1m-dTOF1/m > TOF1Star)
TOF1Star = TOF1m-dTOF1/m;
end
if(TOF1m+dTOF1/m < TOF1End)
TOF1End =  TOF1m+dTOF1/m;
end
if(TOF2m-dTOF2/m > TOF2Star)
TOF2Star = TOF2m-dTOF2/m;
end
if(TOF2m+dTOF2/m < TOF2End)
TOF2End =  TOF2m+dTOF2/m;
end  



end

Time = [DEPm, TOF1m, TOF2m];

[Gs.Dvtot,Gs.fb,Gs.dv,~,~,~] = DvTot(Time(1),Time(2),Time(3),LA,P1,P2,P3,muS,0,R);

Gs.elapsedtime  = toc;

fprintf('Elapsed time: %.4f s   \n',Gs.elapsedtime);
fprintf('dV          = %.4f     \t[km/s]\n',Gs.Dvtot);
fprintf('lambert1    = %.4f     \t[km/s]\n',Gs.dv.lambert1);
fprintf('lambert2    = %.4f     \t[km/s]\n',Gs.dv.lambert2);
fprintf('flyby       = %.5f     \t[km/s]\n',Gs.dv.flyby);
fprintf('rp          = %.4f  \t[km]\n', Gs.fb.rp);
fprintf('dt          = %.4f  \t[d]\n', sum(Time));
end
     