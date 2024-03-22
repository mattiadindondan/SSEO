close all,clear all,clc;
%REVERSE SIZING OF PROPULSION SYSTEM

%%%%%% BIPROPELLANT REVERSE SYZING %%%%%%%%%

deltaV=9900;  % [m/s] deltaV cost
I_sp= 317; % [s] Specific Impulse
m_dry= 1593; % [kg] Dry Mass
OF_ratio= 0.85;  % O/F ratio
R=2077.3; % [J/(Kg*K)] specific constant gas for He
y=1.67; % heat ratio for He
g_0= 9.81;  % [m/s^2]  gravitational constant
rho_hydrazine= 1004.5; % [kg/m^3]  Fuel density
rho_MON=1433;  % [kg/m^3]  Oxidier density
rho_tank= 2810; % [kg/m^3] for Allumium (Al7075) or 2780 for Titanium (Ti6A14V)
sigma=503; % [MPa] for Allumium (Al7075) or 950 for Titanium (Ti6A14V)
%rho_He=0.1784; % [kg/m^3] He density
P_chamber=21; % [bar] Pressure of Combustion Chamber
v_feed=10; % [m/s] check the value of feed velocity
margin_deltaV= 1.1; % margin of deltaV
T_tank=290; %[K] the range for bi-propellant 290-300k Tank Temperature
m_thruster=NaN; % [kg] thruster mass


deltaV_marg=margin_deltaV*deltaV; % deltaV with margin

MR=exp(deltaV_marg/I_sp/g_0);  % Mass Ratio Mox/Mfuel

m_fuel= 1.2*MR;   % [kg] check 1.2 if it is the margin

m_ox=MR*m_fuel;   % [kg] Oxidier Mass

m_prop=m_ox-m_fuel;   % [kg] Propellant Mass
m_prop_real=m_prop*1.055;  % [kg] check the margin

m_fuel_real=m_prop_real/(1+OF_ratio);  % [Kg] Fuel Mass computed with real propellant mass
m_ox_real=OF_ratio*m_fuel_real; % [kg] Oxidier Real Mass 

%%%%COMPUTE THE VOLUME OF OXIDIER AND FUEL WITH MARGIN%%%%%%%%%%
V_fuel=m_fuel_real/rho_hydrazine; % [m^3] Fuel Volume
V_ox=m_ox_real/rho_MON; % [m^3] Oxidier Volume
V_fuel_marg=V_fuel+(0.1*V_fuel);  % [m^3] Fuel Volume + margin
V_ox_marg=V_ox+(0.1*V_ox); % [m^3] Oxidier Volume + margin
V_prop= V_ox_marg+V_fuel_marg; % [m^3] Propellant Volume
V_prop_marg=V_prop+(V_prop*0.1); % [m^3] Propellant Volume + margin


deltaP_inj=0.3*P_chamber; % [bar] deltaP Injection
deltaP_feed=1/2*rho_MON*v_feed^2; % [atm] deltaP feed
P_tank=P_chamber+deltaP_feed+deltaP_inj; % [bar] Tank Pressure
P_press_i=10*P_tank; % [MPa] Pressurize Gas Pressure check the 10 times

m_press=((P_tank*V_prop_marg)/(R*T_tank))*(y/(1-(P_tank/P_press_i))); % [kg] pressurize gas mass 
m_press_real=1.2*m_press; % [kg] pressurize gas mass + margin
V_press=m_press_real*R*T_tank/P_press_i; % [m^3] press Volume

%%% COMPUTE THA TANK MASSES %%%%%
r_tank_ox=(3/4*V_ox_marg/pi)^(1/3);  % [m] radius of spherycal tank of oxidier
r_tank_fuel=(3/4*V_fuel_marg/pi)^(1/3); % [m] radius of spherycal tank of fuel
r_tank_press=(3/4*V_press/pi)^(1/3); % [m] radius of spherycal tank of pressurize gas
t_tank_ox=P_tank*r_tank_ox/sigma;  % [m] thickness of oxidier tank
t_tank_fuel=P_tank*r_tank_fuel/sigma; % [m] thickness of fuel tank
t_tank_press=P_tank*r_tank_press/sigma; % [m] thickness of pressurize gas tank

% compute mass of propulsion system
m_tank_ox=rho_tank*4/3*pi*((r_tank_ox+t_tank_ox)^3-r_tank_ox^3);  % [kg] oxidier tank mass
m_tank_fuel=rho_tank*4/3*pi*((r_tank_fuel+t_tank_fuel)^3-r_tank_fuel^3); % [kg] fuel tank mass
m_tank_press=rho_tank*4/3*pi*((r_tank_press+t_tank_press)^3-r_tank_press^3); % [kg] pressurize gas tank mass


M_propsyst=(1.1*m_tank_fuel)+(1.1*m_tank_ox)+(1.1*m_tank_press)+m_press_real+(1.1*m_thruster); % [kg] total mass of primary propulsion system + margin
