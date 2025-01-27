% Reconstruct spectrum in measurement file using provided calibration file
function [crop_t, crop_p, filtered_p, final_spectrum] = ...
        reconstruct_spectrum(measure_file, calib_file, desired_nu, do_apodize)
    % Using a .mat measurement file and .mat calibration file, return
    % spectrum sampled at desiredNu
    % Because of how two-pol extraction works, all that is needed to do to
    % extract both pol spectra is to call this twice, once with each pol's
    % calibration file

    load(measure_file, 'channel1', 'acq_struct');
    [~,time_array,~] = generate_detector_params(acq_struct);
    measure_acq_struct = acq_struct; clear acq_struct;
    load(calib_file, ...
        'acq_struct', ... % power settings
        'opl_poly', 'start_time', 'end_time', ... % time domain scaling + range
        'g_poly', 'g_mu', 'nu_min', 'nu_max'); % freq domain scaling
    if(~isequal(measure_acq_struct, acq_struct))
        % silence this error if pause_between_calibrations is the only
        % thing different
        if(measure_acq_struct.pause_between_calibrations == acq_struct.pause_between_calibrations)
            warning("Measurement and calibration acquisition structures do not match! If the function succeeds, treat these results with the utmost skepticism...");
        end
    end
    % This warning is useful in general but is disabled here to avoid
    % excess errors in some scripts were we slightly exceed calibration
    % domain on purpose
%     if(max(desired_nu) > nu_max || min(desired_nu) < nu_min)
%         warning("Domain of desiredNu [%f, %f] exceeds calibration domain [%f, %f]. Results outside the calibration domain may be inaccurate.")
%     end
    % crop measurement to region with interferogram
    time_crop_idxs = time_array > start_time & time_array < end_time;
    crop_t = time_array(time_crop_idxs);
    crop_p = channel1(time_crop_idxs);
    [filtered_p, final_spectrum] = recover_spectrum(crop_t, crop_p, ...
        desired_nu, opl_poly, g_poly, g_mu, do_apodize);
end

function [filtered_data, final_spectrum] = recover_spectrum(...
        measured_t, measured_x, desired_nu, ...
        opl_poly, g_poly, g_mu, do_apodize)
    c = 299792458;
    % NUDFT reconstruction using calibration data
    % measuredT are the t-axis sample points
    % measuredX are the optical powers/measurements

    % --- Preprocess interferogram --- %
    % center it to have zero mean
    % deprecated old method of centering uses moving average
    %movmean_length = round(length(measured_x)/10);
    %this_centered = measured_x./ movmean(measured_x,movmean_length) - 1;
    % another method (that is a bit slower but seems better) simply takes a
    % high pass filter 
    this_centered = highpass(measured_x, 0.05);
    % apodize it so we don't get sidelobes in our reconstruction
    this_window = hann(length(this_centered))';
    this_apodized = this_window.*this_centered;

    if(do_apodize)
        filtered_data = this_apodized;
    else
        filtered_data = this_centered;
    end

    % --- Calculate OPL Rescaling --- %
    opl_scaling = polyval(opl_poly, measured_t);
    % --- Calculate Frequency Rescaling --- %
    freq_scaling = polyval(g_poly, desired_nu, [], g_mu).*desired_nu/c;
    % --- NUDFT reconstruction!!! --- %
    final_spectrum = abs(nufft(filtered_data, opl_scaling, freq_scaling));
end