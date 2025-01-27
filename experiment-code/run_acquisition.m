function [time_array, channel1, channel2] = run_acquisition(agi, kes, acq_struct)
    % RUN_ACQUISITION run two-heater push-pull sequence set up using setup_acquisition
    % Parameters that must be in acquisition_structure:
    %     sampling_interval - seconds, how often to sample optically
    %     heater_R - ohms, heater resistance, approx (used for P -> I conversion)
    %     P_range - W, the power range on each heater
    %     P_tot = P_range - W, the total power being put into the chip at any instant in time, usually equals P_range
    %     start_time - seconds, time spent on the start heater point
    %     ramp_time - seconds, time spent on ramp
    %     end_time - seconds, time spent on last heater point
    %     P_num - number of heater (not optical) pts to split ramp into
    %     compliance_current - A
    %     compliance_voltage - V

    % calculate some timing variables
    [total_time,time_array,~] = generate_detector_params(acq_struct);

    % --- Run sweep --- %
    kes_output(kes, true); % turn on power supply output
    agi_arm_logging(agi, TriggerType = "complete"); % arm detector logging
    fprintf("Agilent armed...");
    kes_trig_sweep(kes); % trigger sweep on power supply
    % wait for sweep to finish by checking for completion on detector
    loggingSuccessful = agi_wait_for_logging(agi, ...
        EstLoggingTime = total_time + 3); 
    % regardless of result, turn off outputs
    kes_output(kes, false);
    % check if logging finished and return result if so
    if(loggingSuccessful)
        [channel1, channel2] = agi_get_logging_result(agi);
        agi_reset_triggers(agi);
    else
        error("Logging did not finish in alloted time, and likely never started. Check trigger connection from Keysight to Agilent.");
    end
end