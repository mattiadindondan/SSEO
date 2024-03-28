close all
clear 
clc
%REVERSE SIZING OF PROPULSION SYSTEM

%%%%%% BIPROPELLANT REVERSE SIZING %%%%%%%%%

deltaV=1693;  % [m/s] deltaV cost
I_sp= 317; % [s] Specific Impulse
m_dry= 1593; % [kg] Dry Mass
OF_ratio= 0.85;  % O/F ratio
R=2077.3; % [J/(Kg*K)] specific constant gas for He
y=1.67; % heat ratio for He
g_0= 9.81;  % [m/s^2]  gravitational constant
rho_hydrazine= 1032; % [kg/m^3]  Fuel density
rho_MON=1443;  % [kg/m^3]  Oxidier density
rho_tank= 4430; % [kg/m^3] 2810 for Allumium (Al7075) or 4430 for Titanium (Ti6A14V)
sigma=9500*10^6; % [Pa] 503 for Allumium (Al7075) or 9500 for Titanium (Ti6A14V)
%rho_He=0.1784; % [kg/m^3] He density
P_chamber=18*10^5; % [Pa] Pressure of Combustion Chamber
v_feed=10; % [m/s] check the value of feed velocity
margin_deltaV= 1.1; % margin of deltaV
T_tank=300; %[K] the range for bi-propellant 290-300k Tank Temperature
m_engine=4.5; % [kg] engine mass


deltaV_marg=margin_deltaV*deltaV; % [m/s] deltaV with margin

MR=exp(deltaV_marg/I_sp/g_0);  % Mass Ratio Mox/Mfuel

m_fin= 1.2*m_dry;   % [kg] check 1.2 if it is the margin

m_0=MR*m_fin;   % [kg] Oxidier Mass

m_prop=(m_0-m_fin);   % [kg] Propellant Mass
m_prop_real=m_prop*1.055;  % [kg] check the margin

m_fuel=m_prop_real/(1+OF_ratio);  % [Kg] Fuel Mass computed with real propellant mass
m_ox=OF_ratio*m_fuel; % [kg] Oxidier Real Mass 

%%%%COMPUTE THE VOLUME OF OXIDIER AND FUEL WITH MARGIN%%%%%%%%%%
V_fuel=m_fuel/rho_hydrazine; % [m^3] Fuel Volume
V_ox=m_ox/rho_MON; % [m^3] Oxidier Volume
V_fuel_marg=V_fuel*1.1;  % [m^3] Fuel Volume + margin
V_ox_marg=V_ox*1.1; % [m^3] Oxidier Volume + margin
V_prop= V_ox_marg+V_fuel_marg; % [m^3] Propellant Volume



deltaP_inj=0.3*P_chamber; % [Pa] deltaP Injection
deltaP_feed=(1/2*rho_MON*v_feed^2); % [Pa] deltaP feed
P_tank=(P_chamber+deltaP_feed+deltaP_inj); % [Pa] Tank Pressure
P_press_i=10*P_tank; % [Pa] Pressurize Gas Pressure check the 10 times

m_press=((P_tank*V_prop)/(R*T_tank))*(y/(1-(P_tank/P_press_i))); % [kg] pressurize gas mass 
m_press_real=1.2*m_press; % [kg] pressurize gas mass + margin
V_press=m_press_real*R*T_tank/P_press_i; % [m^3] press Volume

%%% COMPUTE THA TANK MASSES %%%%%
r_tank_ox=(3/4*V_ox_marg/pi/2)^(1/3);  % [m] radius of spherycal tank of oxidier
r_tank_fuel=(3/4*V_fuel_marg/pi/4)^(1/3); % [m] radius of spherycal tank of fuel
r_tank_press=(3/4*V_press/pi/4)^(1/3); % [m] radius of spherycal tank of pressurize gas
t_tank_ox=P_tank*r_tank_ox/(2*sigma);  % [m] thickness of oxidier tank
t_tank_fuel=P_tank*r_tank_fuel/(2*sigma); % [m] thickness of fuel tank
t_tank_press=P_tank*r_tank_press/(2*sigma); % [m] thickness of pressurize gas tank

% compute mass of propulsion system considering 1 tank for each type of
% element
m_tank_ox=rho_tank*4/3*pi*((r_tank_ox+t_tank_ox)^3-r_tank_ox^3);  % [kg] oxidier tank mass
m_tank_fuel=rho_tank*4/3*pi*((r_tank_fuel+t_tank_fuel)^3-r_tank_fuel^3); % [kg] fuel tank mass
m_tank_press=rho_tank*4/3*pi*((r_tank_press+t_tank_press)^3-r_tank_press^3); % [kg] pressurize gas tank mass

rho_tank*((1+P_tank/(2*sigma))^3-1)


%%%%%%% Volume Sizing %%%%%%%
%n=[1:6];
%vfuel=V_fuel_marg./n;
%vox=V_ox_marg./n;

% vfuel = 2*vox %assumption
%rox=(3/4*vox/pi).^(1/3);
%rfuel=(3/4*vfuel/pi).^(1/3);
%figure()
%plot(n,rox)
%hold on
%plot(n,rfuel)
%grid on

%r_tank_ox=(3/4*V_ox_marg/pi)^(1/3);  % [m] radius of spherycal tank of oxidier
%r_tank_fuel=(3/4*V_fuel_marg/pi)^(1/3); % [m] radius of spherycal tank of fuel
%r_tank_press=(3/4*V_press/pi)^(1/3); % [m] radius of spherycal tank of pressurize gas
%tox=P_tank*rox/(2*sigma);  % [m] thickness of oxidier tank
%tfuel=P_tank*rfuel/(2*sigma); % [m] thickness of fuel tank
%t_tank_press=P_tank*r_tank_press/(2*sigma); % [m] thickness of pressurize gas tank

% compute mass of propulsion system
%mtankox=rho_tank*4/3*pi*((rox+tox).^3-rox.^3);  % [kg] oxidier tank mass
%mtankfuel=rho_tank*4/3*pi*((rfuel+tfuel).^3-rfuel.^3); % [kg] fuel tank mass
%m_tank_press=rho_tank*4/3*pi*((r_tank_press+t_tank_press)^3-r_tank_press^3); % [kg] pressurize gas tank mass

% tfuel
% n.*mtankfuel
%mox=vox*rho_MON.*n;
%mfuel=vfuel*rho_hydrazine;

%figure()

%plot(n,tox);
%grid on

r = max(r_tank_ox,r_tank_fuel);
V = 4/3*pi*r^3;
V_ox_final = 2*V;
V_fuel_final = 4*V;


Mass_fuel_final = rho_hydrazine*V_fuel_final;
Mass_ox_final = rho_MON*V_ox_final;


t_tank=P_tank*r/(2*sigma) 
m_tank=rho_tank*4/3*pi*((r+t_tank)^3-r^3)


%%%%%% MONOPROPELLANT REVERSE SIZING %%%%%%%%%
deltaVmono=180; % [m/s] deltaV 
m_dry=1593; % [kg] dry mass
I_sp=226.1; % [s] specific impulse of the thrusters
P_in=27.5*10^5; % [Pa] initial pressure of gas
P_fin=6.2*10^5; % [Pa] final pressure of gas
B=P_in/P_fin; % Blow-Down Ratio
rho_hydrazine=1032; % [kg/m^3] hydrazine density
T_tank=300; % [K] range 290-300 K
R=2077.3; % [J/(kgK)] specific constant gas for He
g_0=9.81; % [m/s^2] gravitational constant
sigma=9500*10^6; % [Pa] 503 for Allumium (Al7075) or 9500 for Titanium (Ti6A14V)
rho_tank=4430; % [kg/m^3] 2810 for Allumium (Al7075) or 4430 for Titanium (Ti6A14V)
m_thrusters=0.49; % [kg] on datasheet  

deltaV_margmono=2*deltaVmono; % deltaV + margin

MR=exp(deltaV_margmono/I_sp/g_0); % Mass Ratio
m_fin=1.2*m_dry; % [kg] final mass
m_0mono=m_fin*MR; % [kg] initial mass
m_propmon=m_0mono-m_fin; % [Kg] propellant mass
m_prop_realmono=m_propmon*1.055; % [kg] propellant mass + margin

V_propmono=m_prop_realmono/rho_hydrazine; % [m^3] propellant volume
V_prop_marg=1.1*V_propmono;  % [m^3] propellant volume + margin

%%%%% COMPUTE THE TOTAL PROPULSION SYSTEM MASS

M_propsyst=1.1*((6*m_tank)+(4*m_tank_press)+m_press_real+m_engine+(12*m_thrusters)) % [kg] total mass of primary propulsion system + margin







