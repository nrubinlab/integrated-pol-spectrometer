% Generate parameters/timing that will set up Agilent detector logging
function [total_time,time_array,sampling_num] = generate_detector_params(acq_struct)
    total_time = acq_struct.start_time + ...
        acq_struct.ramp_time + acq_struct.end_time;
    time_array = (0:acq_struct.sampling_interval:total_time)'; % more sneaky transposes!
    % hacky fix of off-by-one
    time_array = time_array(1:end-1);
    sampling_num = length(time_array);
end

