% calibration_acquisition.m - collect data for calibration

% See spectrometer.m for details on how spectrometer data is collected.
% This file essentially just runs the acquisition routine on a few known
% sources (single laser lines across the range) and saves the data
% It also bundles acquisition parameters into one structure so they are
% passed around easier
% This must be run once for each fiber chuck orientation to get both
% polarization calibration data

% This file is also used for the two-laser resolution and polarization
% proximity experiments, as the logic is essentially identical (just taking
% measurements in a loop and changing laser wavelengths)

clear;
delete (instrfindall); % Delete all existing instruments
MATLAB_EQUIPMENT_CODE_PATH = 'matlab-equipment-code';
addpath(MATLAB_EQUIPMENT_CODE_PATH);

% Initialize and connect Agilent power meter
agi = agilent816x_start(Address = 'GPIB1::20::INSTR'); 
kes = kes_start(); % Connect to keysight power supply
santec = santec_start(); % connect to main tunable laser
anr = anritsu_start(); % connect to OSA
% uncomment this for resolution and pol proximity experiments:
ven = venturi_start(); % connect to secondary tunable laser
%% Heater contact test (check that our probes have good contact)
test_V_compliance = 10; % volts
test_I = 1e-3; % current, amps
kes_config_I_source(kes, test_V_compliance, 1);
kes_set_I(kes, test_I, 1, false);
R1_meas = kes_measure_resistance(kes, 1);
kes_config_I_source(kes, test_V_compliance, 2);
kes_set_I(kes, test_I, 2, false);
R2_meas = kes_measure_resistance(kes, 2);
fprintf("Channel 1 R = %1.1f ohms, Channel 2 R = %1.1f ohms \r\n", ...
    R1_meas, R2_meas);
%% Choose sweep params and send to equipment
% optical settings
detector_range = -20; % power meter range, dBm

% --- Normal calibration settings --- %
% lambda_array_nm = 1480:30:1630; % sweep forward
% do_update_venturi = false;

% --- 2-laser resolution measurement settings --- %
laser_spacing = 3;
lambda_array_nm = (1525:25:1625) - laser_spacing/2;
do_update_venturi = true;

% --- 2-laser polarization proximity settings --- %
% lambda_array_nm = 1540:2.5:1560;
% do_update_venturi = false;

% OSA optical settings are changed infrequently and can be set by hand on
% the unit

% for run_acquisition, we must generate a structure with the following
% fields - if all these are constant, the calibration should be consistent
% see run_acquisition for more details on what these are

acq_struct.sampling_interval = 2e-3; % optical sampling interval
% NOTE: the following resistance is not precise, but does not need to be.
% Errors in this parameter will yield errors in the absolute amount of
% power delivered to the heater, but if the same resistance is used here
% every time, our calibration is still valid.
acq_struct.heater_R = 950; % ohm
acq_struct.P_range = 6; % W
% only change this if we want to try narrower sweeps at elevated temp:
acq_struct.P_tot = acq_struct.P_range; 

acq_struct.start_time = 0.9; % seconds
acq_struct.ramp_time = 10; % seconds
acq_struct.end_time = 0.1; % seconds
acq_struct.P_num = 2001; % # of pts to split POWER ramp into
acq_struct.compliance_current = 85e-3; % A
acq_struct.compliance_voltage = 85; % V

acq_struct.pause_between_calibrations = 15;

% call function to set up all equipment
setup_acquisition(agi, kes, acq_struct, detector_range);
%% Run Calibration
% when we compare OSA data to our spectrometer data, we switch fibers. 
% therefore we run the code once with our spectrometer and the below
% 'false', then again with do_osa = true once we swap fibers
do_osa = false;

% Set Save directory and prefix
[file_prefix, save_dir] = uiputfile('*', ...
    'Select location and prefix where data will be saved:');
if(~file_prefix)
    error("You must select a save location");
end
num_lambda = length(lambda_array_nm);
lambda_array_m = lambda_array_nm*1e-9;

% Warm-up time at start, if desired (not used for any paper data)
% kes_output(kes, true);
% fprintf("Pausing for %d seconds to warm-up...", start_warmup_time);
% pause(start_warmup_time);

for lambda_index = 1:num_lambda
    clear channel1 channel2 lambda
    this_lambda_nm = lambda_array_nm(lambda_index);
    lambda = lambda_array_m(lambda_index);
    % update wavelength
    santec_set_wavelength(santec, this_lambda_nm);
    % second wavelength: only used for resolution measurements
    if(do_update_venturi)
        venturi_set_wavelength(ven, this_lambda_nm + laser_spacing);
    end
    start_timestamp = datetime;
    if(do_osa)
        fprintf("Taking OSA measurement for %1.0f nm...", this_lambda_nm);
        anr_single(anr); % trigger sweep
        anr_wait_for_operation(anr); % wait for sweep to finish
        fprintf("Done.\n")
        [osa_lambda,osa_power_dbm] = anr_get_trace(anr,"A");
        this_filename = sprintf('%s_%d.mat', file_prefix, this_lambda_nm);
        save(fullfile(save_dir, this_filename), ...
            'osa_lambda', 'osa_power_dbm', 'start_timestamp');
    else
        fprintf("Running sweep for %1.0f nm...", this_lambda_nm);
        % call function to run acquisition
        [time_array, channel1, channel2] = run_acquisition(agi, kes,acq_struct);
        this_filename = sprintf('%s_%d.mat', file_prefix, round(this_lambda_nm));
        % save files along the way so we don't lose everything if the
        % acquisition fails/crashes!
        save(fullfile(save_dir, this_filename), ...
            'channel1', 'channel2', 'lambda', 'acq_struct', 'start_timestamp');
    end    
    if lambda_index ~= num_lambda && ~do_osa
        pause(acq_struct.pause_between_calibrations);
    end
end
kes_output(kes, false);