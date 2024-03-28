function ManouverPlot(flag,n1,n2,n3,muS,DEP,TOF1,TOF2,rM)
%
%
% PROTOTYPE:
%  Plot(n1,n2,n3,muS,DEP,TOF1,ARR,rM)
%
% If flag ==1 FILM
% DESCRIPTION:
% 3D plot of the trajectory Manouver with defined TOFs. 
%
% INPUT:
% n1                       [1x1]          Depature ibody                   [-]
% n2                       [1x1]          flyby ibody                      [-]
% n3                       [1x1]          Arrival ibody                    [-]
% muS                      [1x1]          Gravitational parameter          [km^3/s^2]
%
% DEP                      [1x1]          Depature time                    [d]
% TOF                      [1x1]          Time of flight                   [d]
% ARR                      [1x1]          Arrival time                     [d]
% rM                       [nx3]         Manouver trajectory               [km]
% CONTRIBUTORS
%
% Monai Francesco
% Dora Campana
% Arda Varlı
% Marco Barbieri
% Versions: 2023-10-01 First version
[P1.Kep,~] = uplanet(DEP,n1);
kep2car(P1.Kep,muS);
[P1.X, P1.Y, P1.Z] = plotOrbit(P1.Kep,muS,2*pi,0.001);

[P2.Kep,~] = uplanet(DEP+TOF1,n2);
[P2.X, P2.Y, P2.Z] = plotOrbit(P2.Kep,muS,2*pi,0.001);

[P3.Kep,~] = uplanet(DEP+TOF1+TOF2,n3);
[P3.X, P3.Y, P3.Z] = plotOrbit(P3.Kep,muS,2*pi,0.001);


gcf=figure();
set(gcf, 'WindowState', 'maximized')
axis equal
hold on
% orbits
plot3(P1.X, P1.Y, P1.Z,'LineStyle','--','color',[0.4, 0.4, 0.4],'LineWidth',2)
% plot3(P2.X, P2.Y, P2.Z,'LineStyle','--','color',"#4DBEEE",'LineWidth',2)
plot3(P3.X, P3.Y, P3.Z,'LineStyle','--','color',"#D95319",'LineWidth',2)
plot3(rM(:,1), rM(:,2), rM(:,3),'color',"blue",'LineWidth',2)

% manoeuvre positions
plot3(0,0,0,'o','Color',"#EDB120",'MarkerSize',7,'MarkerFaceColor',"#EDB120")
plot3(P1.X(1), P1.Y(1), P1.Z(1),'Marker','o','color','m','MarkerSize',7)
plot3(P2.X(1), P2.Y(1), P2.Z(1),'Marker','o','color','blue','MarkerSize',7)
plot3(P3.X(1), P3.Y(1), P3.Z(1),'Marker','o','color','red','MarkerSize',7)
plot3(0,0,0,'Marker','o','color','yellow','MarkerSize',7)

xlabel('X [km]','FontSize',15)
ylabel('Y [km]','FontSize',15)
zlabel('Z [km]','FontSize',15)
grid on



% planets' positions
if flag ==1
 upd_i = 0;   % Time interval between updates in seconds

pointHandle = plot3(NaN, NaN,NaN, 'Marker','square', 'MarkerSize', 7, 'MarkerFaceColor', 'black');
pointHandle1 = plot3(NaN, NaN,NaN, 'o', 'MarkerSize', 7, 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor','k');
pointHandle2 = plot3(NaN, NaN,NaN, 'o', 'MarkerSize', 7, 'MarkerFaceColor', "#0072BD", 'MarkerEdgeColor',"#77AC30");
pointHandle3 = plot3(NaN, NaN,NaN, 'o', 'MarkerSize', 7, 'MarkerFaceColor', "#EDB120","MarkerEdgeColor","#A2142F");
% time vector
t = linspace(DEP,DEP+TOF1+TOF2,size(rM,1));
% legend('Mercury orbit','Earth orbit','Saturn Orbit','','','','','','','Spacecraft','Mercury','Earth','Saturn');

for i = 1:3:size(rM,1)
    [Kep1,~] = uplanet(t(i),n1);
    [r11,~] = kep2car(Kep1,muS);
      
    [Kep2,~] = uplanet(t(i),n2);
    [r22,~] = kep2car(Kep2,muS);
    
    [Kep3,~] = uplanet(t(i),n3);
    [r33,~] = kep2car(Kep3,muS);


   set(pointHandle,'XData',rM(i,1),'YData', rM(i,2),'ZData',rM(i,3));
   set(pointHandle1,'XData',r11(1),'YData', r11(2),'ZData',r11(3)); 
   set(pointHandle2,'XData',r22(1),'YData', r22(2),'ZData',r22(3));
   set(pointHandle3,'XData',r33(1),'YData', r33(2),'ZData',r33(3));

    grid on;
    
    % Pause for the specified interval to control the speed of the animation
    pause(upd_i);
    
    % Update the figure
    drawnow;
   
end
end

end