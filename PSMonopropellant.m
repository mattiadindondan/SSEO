close all,clear all,clc;
%REVERSE SIZING OF PROPULSION SYSTEM

%%%%%% MONOPROPELLANT REVERSE SIZING %%%%%%%%%
deltaV=180; % [m/s] deltaV
m_dry=1593; % [kg] dry mass
I_sp=226.1; % [s] specific impulse of the thrusters
P_in=27.5*10^5; % [Pa] initial pressure of gas
P_fin=6.2*10^5; % [Pa] final pressure of gas
B=P_in/P_fin; % Blow-Down Ratio
rho_hydrazine=1008; % [kg/m^3] hydrazine density
T_tank=293; % [K] range 290-300 K
R=2077.3; % [J/(kgK)] specific constant gas for He
g_0=9.81; % [m/s^2] gravitational constant
sigma=950*10^6; % [Pa] 503 for Allumium (Al7075) or 950 for Titanium (Ti6A14V)
rho_tank=2780; % [kg/m^3] 2810 for Allumium (Al7075) or 2780 for Titanium (Ti6A14V)
m_thrusters=0.49; % [kg] on datasheet  

deltaV_marg=2*deltaV; % deltaV + margin

MR=exp(deltaV_marg/I_sp/g_0); % Mass Ratio
m_fin=1.2*m_dry; % [kg] final mass
m_0=m_fin*MR; % [kg] initial mass
m_prop=m_0-m_fin; % [Kg] propellant mass
m_prop_real=m_prop*1.055; % [kg] propellant mass + margin

V_prop=m_prop_real/rho_hydrazine; % [m^3] propellant volume
V_prop_marg=1.1*V_prop;  % [m^3] propellant volume + margin

V_in_gas=V_prop_marg/(B-1); % [m^3] initial gas volume

P_tank= P_in; % [Pa] tank pressure
m_press= P_tank*V_in_gas/(R*T_tank); % [Kg] pressurize gas mass
m_press_marg=m_press*1.2; % [kg] pressurize gas mass + margin
V_tank=(V_prop_marg+V_in_gas)*1.01;  % [m^3] volume of tank + 1% of margin
r_tank=(3/4*V_tank/pi)^1/3; % [m] radius of tank
t_tank=P_tank*r_tank/(2*sigma); % [m] thickness of tank
m_tank=rho_tank*4/3*pi*((r_tank+t_tank)^3-r_tank^3);  % [kg] tank's mass

Mtot=1.1*(m_tank+m_press_marg+(12*m_thrusters)); % total mass in kg + margin


