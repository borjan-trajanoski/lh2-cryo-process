% runLoadingSimulation.m
% Main script to simulate LH2 loading from storage tank to trailer
% at the liquid terminal of an 86 TPD liquefaction plant.
%
% Adapted from Petitpas (2018) runNominal.m
%
% Configuration:
%   (ST) or 1 = 350 m3 horizontal storage tank (FEEDING)
%   (ET) or 2 = 60 m3 horizontal LH2 trailer   (RECEIVING)

clc
close all
clear functions  % reset persistent variables in LH2Simulate
tic;

%% 1. Initialize parameters
fprintf('Loading input parameters...\n');
inputs_StorageTankToTrailer;

%% 2. Run simulation
fprintf('Starting ODE simulation...\n');
nominal = LH2Simulate;

%% 3. Extract data
fprintf('Extracting data...\n');
Data_extraction;

%% 4. Plot results
plotLH2Data(nominal);

%% 5. Compute and display transfer loss summary
fprintf('\n========== TRANSFER LOSS SUMMARY ==========\n');

% Final masses in trailer
mL2_final = nominal.mL2(end);
mv2_final = nominal.mv2(end);
m_total_ET_final = mL2_final + mv2_final;
m_total_ET_initial = LH2Model.mL20 + LH2Model.mv20;
m_transferred_to_ET = m_total_ET_final - m_total_ET_initial + nominal.Boiloff_ET(end);

% BOG vented from trailer during fill
BOG_trailer = nominal.Boiloff_ET(end);

% Mass change in storage tank
m_total_ST_initial = LH2Model.mL10 + LH2Model.mv10;
mL1_final = nominal.mL1(end);
mv1_final = nominal.mv1(end);
m_total_ST_final = mL1_final + mv1_final;
m_lost_from_ST = m_total_ST_initial - m_total_ST_final;

% BOG vented from storage tank (during pressurization + depressurization)
% This is tracked via the cumulative Jvvalve1 integral
dt = diff(nominal.t);
Boiloff_ST = sum(dt .* nominal.AAA(2:end)); % AAA = Jvvalve1

fprintf('Storage Tank (ST):\n');
fprintf('  Initial mass:    %.1f kg (%.1f liquid + %.1f vapor)\n', ...
    m_total_ST_initial, LH2Model.mL10, LH2Model.mv10);
fprintf('  Final mass:      %.1f kg (%.1f liquid + %.1f vapor)\n', ...
    m_total_ST_final, mL1_final, mv1_final);
fprintf('  Mass removed:    %.1f kg\n', m_lost_from_ST);
fprintf('  BOG vented (ST): %.1f kg\n', Boiloff_ST);

fprintf('\nTrailer (ET):\n');
fprintf('  Initial mass:    %.1f kg (%.1f liquid + %.1f vapor)\n', ...
    m_total_ET_initial, LH2Model.mL20, LH2Model.mv20);
fprintf('  Final mass:      %.1f kg (%.1f liquid + %.1f vapor)\n', ...
    m_total_ET_final, mL2_final, mv2_final);
fprintf('  LH2 delivered:   %.1f kg\n', m_transferred_to_ET);
fprintf('  BOG vented (ET): %.1f kg\n', BOG_trailer);
fprintf('  Fill fraction:   %.1f%%\n', ...
    nominal.VL2(end)/LH2Model.VTotal2*100);

fprintf('\nTransfer Losses:\n');
fprintf('  --------------------------------------------------\n');
fprintf('  Trailer BOG (pdV compression): %.1f kg (%.2f%% of delivered)\n', ...
    BOG_trailer, BOG_trailer/m_transferred_to_ET*100);
fprintf('  This is the actual transfer loss comparable to\n');
fprintf('  Petitpas (2018) who reported 3.3%% for terminal loading.\n');
fprintf('  --------------------------------------------------\n');
fprintf('  ST vaporizer overhead:         %.1f kg (%.2f%% of delivered)\n', ...
    Boiloff_ST, Boiloff_ST/m_transferred_to_ET*100);
fprintf('  Note: With steady-state manifold (p_ST_final = p_ST_fast),\n');
fprintf('  the ST never depressurizes. This value should be ~0.\n');
fprintf('  Any nonzero ST venting indicates a tuning issue.\n');
fprintf('  --------------------------------------------------\n');

fprintf('\nScaled to 86 TPD Plant (22 trucks/day):\n');
n_trucks = 22;

% Use ONLY trailer-side BOG for the per-truck re-liquefaction calculation
daily_BOG_trailer = BOG_trailer * n_trucks;
fprintf('  Trailer BOG per day:           %.0f kg/day\n', daily_BOG_trailer);
fprintf('  As %% of production:           %.2f%%\n', daily_BOG_trailer/86000*100);
SEC = 6.983; % kWh/kgLH2 from Moro thesis
reliq_power_trailer = daily_BOG_trailer * SEC / 24; % kW
plant_power = 25020; % kW
fprintf('  Re-liquefaction power:         %.1f kW\n', reliq_power_trailer);
fprintf('  As %% of plant power (25 MW):  %.2f%%\n', reliq_power_trailer/plant_power*100);

fprintf('\n  Note: Static storage BOG (502 kg/day from BoilFAST) is NOT\n');
fprintf('  added here. During continuous operation the heat leak into the\n');
fprintf('  manifold is already captured by the simulation (QdotEL1, QdotEV1).\n');
fprintf('  The tanks never reach relief pressure because liquid is withdrawn\n');
fprintf('  every ~65 min. Static BOG only applies during plant downtime.\n');
fprintf('===========================================\n');

%% 6. Save comprehensive output data for Python plotting
% Save as .mat file for easy Python loading via scipy.io
% Also save key time series as CSV for portability

% Collect all key variables into a struct for saving
out.t_s    = nominal.t;            % time in seconds
out.t_min  = nominal.t / 60;       % time in minutes

% Storage tank (ST)
out.mL1    = nominal.mL1;          % liquid mass in ST [kg]
out.mv1    = nominal.mv1;          % vapor mass in ST [kg]
out.TL1    = nominal.TL1(:,end);   % bulk liquid temperature ST [K]
out.Ts1    = nominal.Ts1;          % surface temperature ST [K]
out.Tv1    = nominal.Tv1(:,end);   % bulk vapor temperature ST [K]
out.pv1_Pa = nominal.pv1;          % vapor pressure ST [Pa]
out.pv1_bar = nominal.pv1 / 1e5;   % vapor pressure ST [bar]
out.rhoL1  = nominal.rho_L1;       % liquid density ST [kg/m3]
out.rhov1  = nominal.rhov1;        % vapor density ST [kg/m3]
out.hL1    = nominal.hL1';         % liquid height ST [m]

% Trailer (ET)
out.mL2    = nominal.mL2;          % liquid mass in ET [kg]
out.mv2    = nominal.mv2;          % vapor mass in ET [kg]
out.TL2    = nominal.TL2(:,end);   % bulk liquid temperature ET [K]
out.Ts2    = nominal.Ts2;          % surface temperature ET [K]
out.Tv2    = nominal.Tv2(:,end);   % bulk vapor temperature ET [K]
out.Tw2    = nominal.Tw2;          % wall temperature ET [K]
out.pv2_Pa = nominal.pv2;          % vapor pressure ET [Pa]
out.pv2_bar = nominal.pv2 / 1e5;   % vapor pressure ET [bar]
out.rhoL2  = nominal.rho_L2;       % liquid density ET [kg/m3]
out.rhov2  = nominal.rhov2;        % vapor density ET [kg/m3]
out.hL2    = nominal.hL2;          % liquid height ET [m]
out.VL2    = nominal.VL2;          % liquid volume ET [m3]
out.fillpct = nominal.VL2 / LH2Model.VTotal2 * 100; % fill percentage ET

% Transfer
out.Jtr    = nominal.Jtr;          % transfer mass flow [kg/s]
out.Boiloff_ET = nominal.Boiloff_ET; % cumulative BOG vented from ET [kg]
out.ventstate = nominal.ETTTVenstate; % ET vent valve state

% Vent flows
out.Jvvalve1 = nominal.AAA;       % ST vent mass flow [kg/s]
out.Jvvalve2 = nominal.Jvalve222;  % ET vent mass flow [kg/s]

% Parameters
out.p_ET_high_bar = LH2Model.p_ET_high / 1e5;
out.p_ET_low_bar  = LH2Model.p_ET_low / 1e5;
out.VTotal2 = LH2Model.VTotal2;
out.R2 = LH2Model.R2;

save('output_loading.mat', 'out');
fprintf('Output saved to output_loading.mat\n');
toc;
