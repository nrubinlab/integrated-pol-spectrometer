%% Load a calibration file
clear;

% --- TE-aligned data --- %
% data_dir = 'data\20241014\resolution\TE_3nm'; 
% which_polarization = 'TE'; % which polarization is this data?
% ----------------------- %

% --- TM-aligned data --- %
data_dir = 'data\20241014\resolution\TM_3nm'; 
which_polarization = 'TM';
% ----------------------- %

% Load data from files
file_list = dir(fullfile(data_dir,'*.mat'));
num_files = length(file_list);

osa_dir = 'data\20241014\resolution\OSA_3nm';
% this relies that osa files just happen to be sorted in the same order,
% which works
osa_list = dir(fullfile(osa_dir,'*.mat'));

TE_calib_file = "data\20241014\calibrations\calib\TE.mat";
TM_calib_file = "data\20241014\calibrations\calib\TM.mat";

%% Spectrum reconstructions
% User settings
lambda_range = [1.47e-6, 1.64e-6];
do_apodize = true;
show_osa = true;
do_amplitude_correction = false;

% approximate amplitude ratio TM/TE, calibrated using unpolarized source
pol_factor = 0.36; 

% Spectrum reconstruction
output_lambda = linspace(lambda_range(1),lambda_range(2), 1001);
plot_lambda = 1e9*output_lambda;
c0 = 299792458;
desiredNu = c0./output_lambda;
TE_spectra = zeros(num_files, length(output_lambda));
TM_spectra = zeros(num_files, length(output_lambda));
for fileIdx = 1:num_files
    thisFile = file_list(fileIdx);
    measure_file = fullfile(data_dir, thisFile.name);
    [interfT, interfP, filteredP, TE_reconstruction] = ...
            reconstruct_spectrum(measure_file, TE_calib_file, desiredNu, do_apodize);
    [~, ~, ~, TM_reconstruction] = ...
            reconstruct_spectrum(measure_file, TM_calib_file, desiredNu, do_apodize);
    TE_spectra(fileIdx,:) = TE_reconstruction;
    TM_spectra(fileIdx,:) = TM_reconstruction;
end
% Plotting
te_colors = crameri('devon', num_files+2); % +2 to avoid too-white trace
tm_colors = flipud(crameri('lajolla', num_files+2));
%figure(Units = "centimeters", Position=[3 3 3.8 4]);
figure(Units = "centimeters", Position=[3 3 6 4]);
%figure(Units = "centimeters", Position=[3 3 4 2.3]);
tiledlayout(2,1, Padding="tight", TileSpacing="tight");
nexttile(1); 
hold on;  ylim([0,1]);
nexttile(2);
hold on; ylim([0,1]);
for fileIdx = 1:num_files
    load(fullfile(osa_dir, osa_list(fileIdx).name), 'osa_lambda', 'osa_power_dbm');
    osa_power_lin = 10.^(osa_power_dbm/10);
    this_te = TE_spectra(fileIdx,:);
    this_tm = TM_spectra(fileIdx,:);
    if(strcmp(which_polarization,'TE'))
        normalize_TE = max(this_te);
        normalize_TM = max(this_te)*pol_factor;
    else
        normalize_TE = max(this_tm)/pol_factor;
        normalize_TM = max(this_tm);
    end
    nexttile(1);
    plot(plot_lambda, (this_te/normalize_TE), Color = te_colors(fileIdx,:));
    %title("TE");
    nexttile(2);
    plot(plot_lambda, (this_tm/normalize_TM), Color = tm_colors(fileIdx,:));
    %title("TM");
    
    if(strcmp(which_polarization,'TE'))
            % no point in plotting OSA when we plot cross pol, as there is
        % better than 20 dB extinction ratio and it doesn't show up
        nexttile(1);
    else
        nexttile(2);
    end
    plot(osa_lambda, osa_power_lin/max(osa_power_lin), 'k--');

end
nexttile(1);
xlim([1520,1630]);
xticks([1520,1575,1630]); xticklabels([]);
yticks([0,1]);
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
nexttile(2);
xlim([1520,1630]);
xticks([1520,1575,1630]);
yticks([0,1]);
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);