function [peak_lambdas, peak_widths] = measure_reconstructed_fwhms(data_dir, calib_file)
    % Spectrum reconstructions
    % this is quite an important setting for resolution measurement!!!
    do_apodize = true; % hann window is used for all reconstructions
    
    % Spectrum reconstruction
    lambda_range = [1.46e-6, 1.65e-6];
    output_lambda = linspace(lambda_range(1),lambda_range(2), 1001);
    c0 = physconst('LightSpeed');
    desiredNu = c0./output_lambda;
    file_list = dir(fullfile(data_dir,'*.mat'));
    numFiles = length(file_list);
    peak_lambdas = zeros(1,numFiles);
    peak_widths = zeros(1,numFiles);
    for fileIdx = 1:numFiles
        measure_file = fullfile(data_dir, file_list(fileIdx).name);
        [~,~,~, this_data] = ...
                reconstruct_spectrum(measure_file, calib_file, desiredNu, do_apodize);
        % extract resolution using simple FWHM measurement
        this_data = this_data/max(this_data); % normalize
        % aggresive findpeaks to choose only highest peak
        [~, locs, w, ~] = findpeaks(this_data, output_lambda*1e9, MinPeakHeight=0.9);
        peak_lambdas(fileIdx) = locs;
        peak_widths(fileIdx) = w;
    end
end