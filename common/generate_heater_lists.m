% Generates parameters to send to power supply given params in acq_struct

% Generate list of current values to send to power supply, along with
% supply interval to give to supplies
% also generated is a list of time points for reference
% also, check against compliance limits is performed here
function [I_inc_out, I_dec_out, supply_interval, heater_time] = ...
            generate_heater_lists(acq_struct)
    % see experiment-code/run_acquisition.m for fields expected in acquisition_structure

    supply_interval = acq_struct.ramp_time/(acq_struct.P_num-1); % how often power supply is updated, min 10 us
    start_num = acq_struct.start_time/supply_interval;
    end_num = acq_struct.end_time/supply_interval;

    if(acq_struct.P_range > acq_struct.P_tot)
        error("P_range (%f) greater than P_tot (%f)!", P_range, P_tot);
    end
    P_diff = linspace(-acq_struct.P_range,acq_struct.P_range,acq_struct.P_num);
    P_inc = acq_struct.P_tot/2 + P_diff/2;
    P_dec = acq_struct.P_tot/2 - P_diff/2;
    I_inc = sqrt(P_inc/acq_struct.heater_R);
    I_dec = sqrt(P_dec/acq_struct.heater_R);
    start_ones = []; end_ones = [];
    if(start_num > 0)
        start_ones = ones(1,start_num);
    end
    if(end_num > 0)
        end_ones = ones(1,end_num);
    end
    I_inc_out = [I_inc(1)*start_ones, I_inc, I_inc(end)*end_ones];
    I_dec_out = [I_dec(1)*start_ones, I_dec, I_dec(end)*end_ones];
    total_time = acq_struct.start_time + acq_struct.ramp_time + acq_struct.end_time;
    heater_time = 0:supply_interval:total_time;
    
    % check output arrays against compliance settings
    fprintf("---\n");
    if(any(I_inc_out > acq_struct.compliance_current) || any(I_dec_out > acq_struct.compliance_current))
        error("Power/resistance settings will exceed compliance current! (I_max = %f and %f)", max(I_inc_out), max(I_dec));
    else
        fprintf("I_max = %f and %f\n", max(I_inc_out), max(I_dec_out));
    end
    if(any(acq_struct.heater_R*I_inc_out > acq_struct.compliance_voltage) || any(acq_struct.heater_R*I_dec_out > acq_struct.compliance_voltage))
        error("Power/resistance settings will exceed compliance voltage! (R*I_max = %f and %f)", ...
            acq_struct.heater_R*max(I_inc_out), acq_struct.heater_R*max(I_dec_out));
    else
        fprintf("R*I_max = %f and %f\n", acq_struct.heater_R*max(I_inc_out), acq_struct.heater_R*max(I_dec_out));
    end
end