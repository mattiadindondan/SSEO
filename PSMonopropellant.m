close all,clear all,clc;
%REVERSE SIZING OF PROPULSION SYSTEM

%%%%%% MONOPROPELLANT REVERSE SIZING %%%%%%%%%
deltaV=NaN;
m_dry=NaN;
I_sp=NaN;
P_in=NaN;
P_fin=NaN;
B=P_in/P_fin; % Blow-Down Ratio
rho_hydrazine=1008; % [kg/m^3] hydrazine density
g_0=9.81; % [m/s^2] gravitational constant


deltaV_mol=NaN;

MR=exp(deltaV_mol/I_sp/g_0);
m_fin=1.2*m_dry;
m_0=m_fin*MR;
m_prop=m_0-m_fin;
m_prop_real=m_prop*1.055;

V_prop=m_prop_real/rho_hydrazine;
V_prop_marg=1.1*V_prop;

V_in_gas=V_prop_marg


