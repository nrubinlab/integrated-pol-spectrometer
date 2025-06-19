% Plot polarization ratio versus angle of linear polarization incident on
% chip
%% Load a calibration file
% for these plots, use calibration files in data/20250525/calibs_1
clear;
[~, TE_calib_file, TM_calib_file, amp_calib_file] = prompt_load_calibration();

%% Load measurement files
% data is in data/20250525/pol_sweep - select all files in that dir
[file, location] = uigetfile('.mat', 'Select One or More Files', 'MultiSelect', 'on');
if(iscell(file))
    numFiles = length(file);
else
    numFiles = 1;
    file = {file};
end
%% Spectrum reconstructions
% User settings
lambda_range = [1.47e-6, 1.64e-6];
do_apodize = true;
show_osa = true;

% Spectrum reconstruction
output_lambda = linspace(lambda_range(1),lambda_range(2), 1001);
plot_lambda = 1e9*output_lambda;
c0 = 299792458;
desiredNu = c0./output_lambda;
TE_spectra = zeros(numFiles, length(output_lambda));
TM_spectra = zeros(numFiles, length(output_lambda));
angle_list = zeros(numFiles, 1);
for fileIdx = 1:numFiles
    measure_file = fullfile(location, file{fileIdx});
    [interfT, interfP, filteredP, TE_reconstruction] = ...
            reconstruct_spectrum(measure_file, TE_calib_file, desiredNu, ...
            do_apodize);
    [~, ~, ~, TM_reconstruction] = ...
            reconstruct_spectrum(measure_file, TM_calib_file, desiredNu, ...
            do_apodize);
    TE_spectra(fileIdx,:) = TE_reconstruction;
    TM_spectra(fileIdx,:) = TM_reconstruction;
    [~,num_str,~] = fileparts(file{fileIdx});
    num_str_split = split(num_str, '_');
    num_str = num_str_split(1);
    angle_list(fileIdx) = str2double(num_str);
end
%% Plotting
% Plot spectra just for reference - slight calibration drift is visible in
% these plots
colors = cool(numFiles);
figure(Units = "inches", Position=[3 3 6 4]);
tiledlayout(2,1);
nexttile(1); title("TE");
nexttile(2); title("TM"); ylim([1e-6,1e-3]);

colororder(colors);
for fileIdx = 1:numFiles
    nexttile(1);
    semilogy(plot_lambda, abs(TE_spectra(fileIdx,:))); hold on;
    nexttile(2);
    semilogy(plot_lambda, abs(TM_spectra(fileIdx,:))); hold on;
end
nexttile(1);
ylim([1e-6,1e-3]); ylabel("Reconstructed TE (a.u.)");
nexttile(2);
ylim([1e-6,1e-3]); ylabel("Reconstructed TM (a.u.)");
%% Power ratio vs angle
% integrate power over some window and take ratio between polarizations
% ratio method is used so misalignment errors don't show up in the
% measurement
% note that due to poor calibration, the laser shows up at 1580 on the TE
% channel but 1575 on the TM. Fortunately this is irrelevant for the
% experiment here (and if necessary, the experiment can be re-done with
% proper calibration to show this)
TE_center = 1580;
TM_center = 1576;
int_radius = 1;
TE_powers = zeros(numFiles, 1);
TM_powers = zeros(numFiles, 1);

for fileIdx = 1:numFiles
    TE_powers(fileIdx) = powerInRegion(plot_lambda, ...
        TE_spectra(fileIdx,:), TE_center - int_radius, TE_center + int_radius);
    TM_powers(fileIdx) = powerInRegion(plot_lambda, ...
        TM_spectra(fileIdx,:), TM_center - int_radius, TM_center + int_radius);
end
% shift and wrap angle to make plot x axis nice
zero_angle = 17;
shift_angle = 0; % at what angle do we put the zero
shift_under_mod = -45;
power_ratio = 0.36;
angle_wrapped = shift_under_mod + mod(angle_list - zero_angle + shift_angle - shift_under_mod, 180);
figure("Units", "centimeters", "Position", [5 1 8 6]); hold on;

% theoretical curve
angle_theory = linspace(shift_under_mod,shift_under_mod+180,1000);
plot(angle_theory,10*log10((tand(angle_theory - shift_angle)).^2), 'r');
% experimental curve
plot(angle_wrapped, 10*log10(TE_powers ./ TM_powers), 'ko', MarkerSize = 4);
yline(-30:10:30, 'k:');
ylim([-40,40]); ylabel("TE / TM ratio (dB)");
xticks(shift_under_mod + 0:45:180); xlabel("Fiber angle (deg)");
xlim('tight');
%plot(angle_wrapped, 10*log10(TM_powers), 'o');
legend("Ideal polarizer", "Experiment", Location = "northwest");
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);

function out = powerInRegion(lambda, psd, low, high)
    crop_idx = (lambda > low) & (lambda < high);
    lambda_crop = lambda(crop_idx);
    psd_crop = psd(crop_idx);
    out = trapz(lambda_crop, psd_crop);
end