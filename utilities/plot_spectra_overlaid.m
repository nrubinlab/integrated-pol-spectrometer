% Plot many reconstructed spectra on top of each other
% Just so I don't have to repeat code, I also use this file for empirical 
% resolution extraction
%% Load a calibration file
clear;
[~, TE_calib_file, TM_calib_file, amp_calib_file] = prompt_load_calibration();
%% Load measurement files
[file, location] = uigetfile('.mat', 'Select One or More Files', 'MultiSelect', 'on');
if(iscell(file))
    numFiles = length(file);
else
    numFiles = 1;
    file = {file};
end
% %% Load an OSA file, if you choose
% [osa_file, osa_location] = uigetfile('.mat', 'Select OSA measurement');
% load(fullfile(osa_location, osa_file), ...
%         'osa_lambda', 'osa_power_dbm');
% osa_power_lin = 10.^(osa_power_dbm/10);
% osa_power_lin = osa_factor*osa_power_lin/max(osa_power_lin);
%% Spectrum reconstructions
% User settings
lambda_range = [1.47e-6, 1.64e-6];
do_apodize = true;
show_osa = true;
do_amplitude_correction = false;

% Spectrum reconstruction
output_lambda = linspace(lambda_range(1),lambda_range(2), 1001);
plot_lambda = 1e9*output_lambda;
c0 = 299792458;
desiredNu = c0./output_lambda;
TE_spectra = zeros(numFiles, length(output_lambda));
TM_spectra = zeros(numFiles, length(output_lambda));
for fileIdx = 1:numFiles
    measure_file = fullfile(location, file{fileIdx});
    [interfT, interfP, filteredP, TE_reconstruction] = ...
            reconstruct_spectrum(measure_file, TE_calib_file, desiredNu, do_apodize);
    [~, ~, ~, TM_reconstruction] = ...
            reconstruct_spectrum(measure_file, TM_calib_file, desiredNu, do_apodize);
    if(do_amplitude_correction)
        [TE_reconstruction, TM_reconstruction] = correct_amplitude(...
        output_lambda,TE_reconstruction,TM_reconstruction,amp_calib_file);
    end

    TE_spectra(fileIdx,:) = TE_reconstruction;
    TM_spectra(fileIdx,:) = TM_reconstruction;
end
% Plotting
figure(Units = "inches", Position=[3 3 6 4]);
tiledlayout(2,1);
nexttile(1); hold on; title("TE");
nexttile(2); hold on; title("TM");
for fileIdx = 1:numFiles
    nexttile(1);
    plot(plot_lambda, TE_spectra(fileIdx,:));
    nexttile(2);
    plot(plot_lambda, TM_spectra(fileIdx,:));
end
if(exist('osa_power_dbm','var') && show_osa)
    nexttile(1);
    plot(osa_lambda, osa_power_lin, 'k--');
end