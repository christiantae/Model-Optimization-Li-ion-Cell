%  Parameter identification of an equivalent circuit battery model 
%  over a HPPC profile at T = 23degC.
%
%  Parameters to be identified: R0, R1, C1.  
%  Parameters are identified as a function of SOC
%
%  Objective function is to minimize error between simulated and measured
%  voltage. Optimization algorithm: fminsearch.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc

tic % Monitor time taken to run the code (tic - toc)
load OCV_SOC_curve % Load OCV vs SOC curve for T=23degC
load NMC_Cell_H4_T23_HPPC % Load experimental current and voltage data

%% Global Parameters

global I_pulse_values V_pulse_values ...
    SOC_pulse_range temp_battery dt

%% Initialization of Variables

Current_data = -NMC_Cell_H4_T23_HPPC(:,1);          % Load current data
Voltage_data = NMC_Cell_H4_T23_HPPC(:,2);           % Load voltage data
dt = 1;                                                 % Sampling time
Time_data = [0:dt:(length(Current_data)-1)*dt]';    % Set timestamp data
N_opt = 3;                                              % Number of variables ...
                                                        % to be optimized
datapoints = 8;                                         % Number of pulses
temp_battery = 23;                                      % Select temperature

%% SOC Calculation

capacity = interp1(temperature_points,capacity_points,temp_battery); % Q(T)
SOC_cc(1) = 1; % Start from a fully charged battery
for k=2:length(Current_data)
    SOC_cc(k) = SOC_cc(k-1) - Current_data(k) * dt / (capacity*3600);
    if SOC_cc(k)<0
       SOC_cc(k)=0;
    end
end

%% Data Extraction 
%  For identification, select only charge-discharge pulses at every 10% SOC

Current_magnitude = 3.9;        % Magnitude of discharge pulse
start_point = 1;                % Temporary variable
data_padding = 20;              % Ensure each pulse is long enough for identification
pulse_duration = 100;           % Discharge pulse duration

for i = 1:datapoints
    
    % Select discharge part
    index_start_dischg(i) = find(Current_data(start_point:end)>=Current_magnitude,1)+start_point-data_padding;
    index_end_dischg(i) = index_start_dischg(i) + pulse_duration;
    
    % Update temporary variable with new value to select next pulse
    start_point = index_end_dischg(i);
     
    % Store Current and Voltage of discharge and charge pulses
    I_pulse_values(:,i) = Current_data(index_start_dischg(i):index_end_dischg(i),1);
    V_pulse_values(:,i) = Voltage_data(index_start_dischg(i):index_end_dischg(i),1);
     
    % Store SOC values for the pulses
    SOC_pulse_range(i) = SOC_cc(index_start_dischg(i));
     
%     figure();plot(V_pulse_values); % Plot the selected pulses
end

%% Optimization Routine for Parameter Identification

% Use Initial Guess for parameters from the graphical method
load ECM_graphical_T23.mat % Graphical method values
x0 = ECM_graphical_T23;

%options = optimset('PlotFcns', @optimplotfval) %For nelder mead 
% Run optimization routine in a loop for all pulses
for i = 1:datapoints %1:datapoints
    dataset = i;  
%     R0 = optimizableVariable('R0', [0,1],'Type','real'); 
%     R1 = optimizableVariable('R1', [0,1],'Type','real'); 
%     C1 = optimizableVariable('C1', [900, 3000],'Type','real'); 
%     fun = @(x)ECM_fit_BO(x, i);
%     results = bayesopt(fun, [R0, R1, C1], 'AcquisitionFunctionName', 'expected-improvement', 'IsObjectiveDeterministic',true) 
%     x_val(:,i) = table2array(results.XAtMinObjective)'
%     err_val(:,i) = 
    %x = bestPoint(results,
    %[x_val(:,i), err_val(:,i)] = fminsearch(@(x1)ECM_fit(x1,dataset),x0(:,i)', options);  
    [x_val(:,i), err_val(:,i)] = grad_descent(x0(:,i)', dataset); 
    %[x_val(:,i), err_val(:,i)] = graphicalmethod(x0(:,i)', dataset) 
end

% Parse the optimized parameter values
R0_val_T23 = fliplr(x_val(1,:));
R1_val_T23 = fliplr(x_val(2,:));
C1_val_T23 = fliplr(x_val(3,:));
save('OPT_Params_T23_GD_train_error', 'err_val')
save('OPT_Params_T23_GD.mat','SOC_pulse_range','R0_val_T23','R1_val_T23','C1_val_T23')
toc % Monitor time taken to run the code (tic - toc)