function setup_acquisition(agi, kes, acq_struct, detector_range)
    % Set up hardware to perform acquisition using run_acquisition.
    % Acqusition method: "push-pull" heating on two arms of MZI, using
    % fast user sequence functionality on Agilent B2962a.
    % Hardware trigger used to start power readings on power meter in Agilent
    % 8164b system
    
    % See run_acquisition for parameters that must be in acq_struct

    % --- Setup Keysight power supply --- %
    % generate list of power points and time array
    [I_inc, I_dec, supply_interval, heater_time] = generate_heater_lists(acq_struct);
    % this line is un-commented for inbalanced supply tests
    %[I_inc, I_dec, supply_interval, heater_time] = generate_heater_lists_unbalanced(acq_struct);

    % "background" currents - what is on when the user sequence isn't
    % running? i.e., when the sweep finishes, the power supply will revert
    % to this
    background_index = 1;
    % uncomment this if you want both heaters to be equal
    %[~, background_index] = min(abs(I_inc-I_dec)); 
    background_current_1 = I_inc(background_index);
    background_current_2 = I_dec(background_index);
    
    % setup channel 1 as increasing
    kes_setup_user_sequence(kes, 1, 'current', ... 
                I_inc, supply_interval, acq_struct.compliance_voltage);
    % last argument here 'false' is a flag to use amps instead of mA - be
    % careful!
    kes_set_I(kes, background_current_1, 1, false);

    % setup channel 2 as decreasing
    kes_setup_user_sequence(kes, 2, 'current', ... 
                I_dec, supply_interval, acq_struct.compliance_voltage);
    kes_set_I(kes, background_current_2, 2, false);
    

    % --- Setup Agilent optical power meter --- %
    % parameters for Agilent power meter
    [total_time,~,sampling_num] = generate_detector_params(acq_struct);
    agi_setup_logging(agi, sampling_num, ...
        DetectorIntTime=acq_struct.sampling_interval);
    agi_set_range(agi, detector_range, DetectorChannel = 1);

    % --- Print message + figure for user with some sweep params --- %
    power_per_optical_sample = acq_struct.P_range/acq_struct.ramp_time * acq_struct.sampling_interval;
    figure; hold on;
    plot(heater_time, I_inc);
    plot(heater_time, I_dec);
    hold off; xlabel("Time (s)"); ylabel("Current (A");
    fprintf("Sweep time %1.1f s, %d optical samples at %.3e mW per optical sample\n", ...
        total_time, sampling_num, 1e3*power_per_optical_sample);
end

