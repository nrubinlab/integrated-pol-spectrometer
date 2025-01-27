%% Plots showing spectral reconstructions of cross-pol lasers, while one is
% kept constant and the other is tuned.
% this also gives some indication of calib repeatability + accuracy
% moving-TE data
clear;
data_dir = 'data\20241014\pol_proximity\pol_prox_TE_try2';
TE_calib_file = "data\20241014\calibrations\calib3\TE.mat";
TM_calib_file = "data\20241014\calibrations\calib3\TM.mat";
% Load data from files
file_list = dir(fullfile(data_dir,'*.mat'));

num_files = length(file_list);
%% Spectrum reconstructions
% User settings
%lambda_range = [1.47e-6, 1.64e-6];
lambda_range = [1.53e-6, 1.57e-6];
lambda_ticks = [1530, 1550, 1570];
do_apodize = true;
show_osa = true;
do_amplitude_correction = false;

% Spectrum reconstruction
output_lambda = linspace(lambda_range(1),lambda_range(2), 1001);
plot_lambda = 1e9*output_lambda;
c0 = 299792458;
desiredNu = c0./output_lambda;
TE_spectra = zeros(num_files, length(output_lambda));
TM_spectra = zeros(num_files, length(output_lambda));
TE_reconst_center = zeros(num_files, 1);
TM_reconst_center = zeros(num_files, 1);
TE_reconst_fwhm = zeros(num_files, 1);
TM_reconst_fwhm = zeros(num_files, 1);

for fileIdx = 1:num_files
    thisFile = file_list(fileIdx);
    measure_file = fullfile(data_dir, thisFile.name);
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

    [TE_reconst_center(fileIdx), TE_reconst_fwhm(fileIdx)] ...
        = measure_peak(output_lambda, TE_reconstruction);
    [TM_reconst_center(fileIdx), TM_reconst_fwhm(fileIdx)] ...
        = measure_peak(output_lambda, TM_reconstruction);
end
%% Plotting
figure(Units = "centimeters", Position=[3 3 7 5.5]);
tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
nexttile(1); hold on; title("TE");
xticks(lambda_ticks);xticklabels([]);
xlim('tight');

nexttile(3); hold on; title("TM");
xticks(lambda_ticks); 
xlim('tight');


te_colors = crameri('devon', num_files+2); % +2 to avoid too-white trace
tm_colors = flipud(crameri('lajolla', num_files+2));
vert_offset = 0.2;
for fileIdx = 1:num_files
    nexttile(1);
    set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
    plot(plot_lambda, fileIdx*vert_offset + TE_spectra(fileIdx,:)./max(TE_spectra(fileIdx,:)), Color = te_colors(fileIdx,:));
    
    nexttile(3);
    set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
    plot(plot_lambda, fileIdx*vert_offset + TM_spectra(fileIdx,:)./max(TM_spectra(fileIdx,:)), Color = tm_colors(fileIdx,:));
end

nexttile([2 1]);
hold on;
for fileIdx = 1:num_files
    mkr_size = 4;
    % put a 0.5 factor on error bars as they're the fwhm, not hwhm
    errorbar(fileIdx, 1e9*TE_reconst_center(fileIdx), ...
        0.5*1e9*TE_reconst_fwhm(fileIdx), ...
        'o', ...
        Color = te_colors(fileIdx,:), ...
        MarkerSize=mkr_size);
    errorbar(fileIdx, 1e9*TM_reconst_center(fileIdx), ...
        0.5*1e9*TM_reconst_fwhm(fileIdx), ...
        'o', ...
        Color = tm_colors(fileIdx,:), ...
        MarkerSize=mkr_size);
end

ylim([1535 1565]);
xlabel("Test no.");
ylabel("Peak λ");
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);

function [lambda, fwhm] = measure_peak(lambda_in, power_in)
    [~, max_idx] =  max(power_in);
    lambda = lambda_in(max_idx);
    threshold = (min(power_in) + max(power_in)) / 2;
    index1 = find(power_in >= threshold, 1, 'first');
    index2 = find(power_in >= threshold, 1, 'last');
    fwhm = lambda_in(index2) - lambda_in(index1);
end