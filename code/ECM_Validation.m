%% Energy 295 Spring 2019 - Homework 3, Problem 2
%  Validation of an equivalent circuit battery model 
%  over HPPC profiles at T = 23degC and 45degC, 
%  and a US06 profile at T = 23degC.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc

%% Prepare lookup tables for ECM parameters

% Load R0, R1, C1 as a function of SOC for T = 23 and 45degC
%load OPT_Params_T23_GradientDescent
load OPT_Params_T23_GD
load OPT_Params_T45

R0_values = [R0_val_T23; R0_val_T45];
R1_values = [R1_val_T23; R1_val_T45];
C1_values = [C1_val_T23; C1_val_T45];
SOC_values = fliplr(SOC_pulse_range);
Temp_values = [23, 45];

figure(); surf(SOC_values,Temp_values,R0_values)
xlabel('SOC [-]','FontSize', 16); zlabel('R_{0} [\Omega]','FontSize', 16)
ylabel('Temperature [^{o}C]','FontSize', 16); set(gca,'FontSize', 16)
title('R0 as a function of SOC and Temperature','FontSize', 16)

save('ECM_Params_Opt.mat','R0_values','R1_values','C1_values',...
    'SOC_values','Temp_values')

%% Validate against HPPC Profile at T = 45degC

clear all; warning off
load OCV_SOC_curve % Load OCV vs SOC curves
load ECM_Params_Opt % ECM Parameters
load NMC_Cell_H4_T45_HPPC % Load experimental current and voltage data

I_data = -NMC_Cell_H4_T45_HPPC(:,1);          % Load current data
V_data = NMC_Cell_H4_T45_HPPC(:,2);           % Load voltage data
dt = 1;                                       % Sampling time
t_data = [0:dt:(length(I_data)-1)*dt]';       % Set timestamp data
initialSoC = 1;                               % Set Initial SOC
Tend = t_data(end);                           % Set simulation end time
Tstart = t_data(1);                           % Set simulation start time
temp_battery = 45;                            % Temperature

options = simset('SrcWorkspace','current'); 
sim('ECM_Final',[Tend],options);

% Compute RMSE for identification profile
RMSE_V_HPPC_T45  = rms(v_expt-vbatt)*100/mean(v_expt)

figure(); hold all; grid on
plot(clk,v_expt,'b','LineWidth',2); plot(clk,vbatt,'--r','LineWidth',2)
ylabel('Voltage [V]','FontSize', 16); xlabel('Time [s]','FontSize', 16)
set(gca,'FontSize', 16); legend('Experimental','Simulated','FontSize', 16)
title('HPPC (T=45^{o}C) Voltage Response','FontSize', 16)
xlim([0 t_data(end)]); ylim([3.4 4.25])

%% Validate against HPPC Profile at T = 23degC

clear all; warning off
load OCV_SOC_curve % Load OCV vs SOC curves
load ECM_Params_Opt % ECM Parameters
load NMC_Cell_H4_T23_HPPC % Load experimental current and voltage data

I_data = -NMC_Cell_H4_T23_HPPC(:,1);          % Load current data
V_data = NMC_Cell_H4_T23_HPPC(:,2);           % Load voltage data
dt = 1;                                       % Sampling time
t_data = [0:dt:(length(I_data)-1)*dt]';       % Set timestamp data
initialSoC = 1;                               % Set Initial SOC
Tend = t_data(end);                           % Set simulation end time
Tstart = t_data(1);                           % Set simulation start time
temp_battery = 23;                            % Temperature

options = simset('SrcWorkspace','current'); 
sim('ECM_Final',[Tend],options);

% Compute RMSE for identification profile
RMSE_V_HPPC_T23  = rms(v_expt-vbatt)*100/mean(v_expt)

figure(); hold all; grid on
plot(clk,v_expt,'b','LineWidth',2); plot(clk,vbatt,'--r','LineWidth',2)
ylabel('Voltage [V]','FontSize', 16); xlabel('Time [s]','FontSize', 16)
set(gca,'FontSize', 16); legend('Experimental','Simulated','FontSize', 16)
title('HPPC (T=23^{o}C) Voltage Response: Gradient Descent','FontSize', 16)
xlim([0 t_data(end)]); ylim([3.3 4.25])

%% Validate against US06 Profile at T = 23degC

clear all; warning off
load OCV_SOC_curve % Load OCV vs SOC curves
load ECM_Params_Opt % ECM Parameters
load NMC_Cell_H4_T23_US06 % Load experimental current and voltage data

I_data = -NMC_Cell_H4_T23_US06(:,1);          % Load current data
V_data = NMC_Cell_H4_T23_US06(:,2);           % Load voltage data
dt = 0.1;                                       % Sampling time
t_data = [0:dt:(length(I_data)-1)*dt]';       % Set timestamp data
initialSoC = interp1(ocv_points(1,:),soc_points,V_data(1));% Set Initial SOC
Tend = t_data(end);                           % Set simulation end time
Tstart = t_data(1);                           % Set simulation start time
temp_battery = 23;                            % Temperature

options = simset('SrcWorkspace','current'); 
sim('ECM_Final',[Tend],options);

% Compute RMSE for identification profile
RMSE_V_US06_T23  = rms(v_expt-vbatt)*100/mean(v_expt)

figure(); hold all; grid on
plot(clk,v_expt,'b','LineWidth',2); plot(clk,vbatt,'--r','LineWidth',2)
ylabel('Voltage [V]','FontSize', 16); xlabel('Time [s]','FontSize', 16)
set(gca,'FontSize', 16); legend('Experimental','Simulated','FontSize', 16)
title('US06 Voltage Response: Gradient Descent','FontSize', 16)
xlim([0 t_data(end)]); ylim([3.55 3.9])