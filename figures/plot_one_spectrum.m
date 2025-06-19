%% Plot a single spectrum extraction for 
% Load a calibration directory
clear;
% this is the easiest way way to pick calibration files
[~, TE_calib_file, TM_calib_file] = prompt_load_calibration(); 
%% Load measurement folder
% for simplicity, data consists of files named 'measurement', 'OSA_TE', and
% 'OSA_TM' in a single folder
% usually, uigetdir() is used so the user can pick files with a UI.
% However, to make paper figure plots easy to reproduce, there is a toggle
% to turn this on/off

do_ui = false;

if(do_ui)
    measure_dirname = uigetdir(); % this is the easier way to pick datasets
end
% Datasets shown in figure (SET do_ui = false and uncomment!!!)

% a, 
% measure_dirname ='data\20241014\spectra_laser\1480TE_1630TM';
% normalize_type = "same";

% b
% measure_dirname ='data\20241014\spectra_ASE\ASE_TM';
% normalize_type = "TM";

% c
% measure_dirname ='data\20241014\spectra_ASE\ASE_unpol';
% normalize_type = "TM";

% d
% measure_dirname ='data\20241014\spectra_SLED\SLED_45_beating_1';
% normalize_type = "TE";

% e, normalize_type = "same"
measure_dirname ='data\20241014\spectra_laser\1480TM_1630TE';
normalize_type = "same";

% f
% measure_dirname ='data\20241014\spectra_SLED\SLED_TE_1480_TM';
% normalize_type = "same";

% g
% measure_dirname ='data\20241014\spectra_ASE\ASE_unpol_1630TM';
% normalize_type = "TE";

% h
% measure_dirname ='data\20241014\spectra_SLED\SLED_45_beating_2';
% normalize_type = "TE";


measure_file = fullfile(measure_dirname, 'measurement.mat'); % our spectrometer data
OSA_TE_file = fullfile(measure_dirname, 'OSA_TE.mat'); % OSA data with TE-oriented analyzer
OSA_TM_file = fullfile(measure_dirname, 'OSA_TM.mat'); % OSA data with TM-oriented analyzer
if(~isfile(measure_file))
    error("No measure.mat in chosen dir!");
end
osa_te_exists = false;
osa_tm_exists = false;
if(isfile(OSA_TE_file))
    load(OSA_TE_file, 'osa_lambda', 'osa_power_dbm');
    osa_te_lambda = osa_lambda;
    osa_te_dbm = osa_power_dbm;
    osa_te_lin = 10.^(osa_te_dbm/10);
    clear osa_lambda osa_power_dbm
    osa_te_exists = true;
end

if(isfile(OSA_TM_file))
    load(OSA_TM_file, 'osa_lambda', 'osa_power_dbm');
    osa_tm_lambda = osa_lambda;
    osa_tm_dbm = osa_power_dbm;
    osa_tm_lin = 10.^(osa_tm_dbm/10);
    clear osa_lambda osa_power_dbm
    osa_tm_exists = true;
end
%% Plot extracted spectrum - simple two-pane interferogram/reconstruction
% settings for figure in paper:
% to what polarization do we normalize the intensity? for accurate pol
% intensity ratios, do not pick "same"
if(do_ui)
    % only change normalize type here if we're in ui-selected data mode
    normalize_type = "TM"; % TE, TM, or same
end
pol_factor = 0.36; % calibration factor from looking at unpolarized ASE

show_axis_labels = false; % we add these manually to first plot, it's easiest
lin_range = [-0.1, 1.2];
lambda_range = [1.47e-6, 1.64e-6];
%figsize = [4.7 5]; % cm, for the first pane
figsize = [4.7 5.2]; % cm

do_apodize = true;
do_log = false;
osa_te = true;
osa_factor = 1;
show_interferogram = true;
interferogram_decimate = 1;
show_osa = true;
these_x_ticks = 1480:40:1640;
these_x_labels = ["1480", "", "", "", "1640"];

TEcolor = '#3081D0'; TMcolor = '#B31312';
% prepare x-axis
output_lambda = linspace(lambda_range(1),lambda_range(2), 1001);
plot_lambda = 1e9*output_lambda;
c0 = 299792458;
desiredNu = c0./output_lambda;
tic
[interfT, interfP, filteredP, TE_reconstruction] = ...
        reconstruct_spectrum(measure_file, TE_calib_file, desiredNu, do_apodize);
[~, ~, ~, TM_reconstruction] = ...
        reconstruct_spectrum(measure_file, TM_calib_file, desiredNu, do_apodize);
% plotting
figure("Units", "centimeters", "Position", [5 1 figsize]);
if(show_interferogram), num_tiles = 3; else, num_tiles = 2; end
t = tiledlayout(num_tiles,1,'TileSpacing','Compact','Padding','Compact');   
% interferogram plot
if(show_interferogram)
    nexttile;
    plot_t = decimate(interfT - min(interfT),interferogram_decimate);
    plot_p = decimate(interfP/max(interfP),interferogram_decimate);
    plot(plot_t, plot_p, 'k');
    % xlabel("Time (s)");
    xticks(0:2:10); xticklabels([]);
    yticklabels([]);
    % ylabel("Normalized Interferogram (a.u.)");
    % title("Preprocessed interferogram");
    set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
end

osa_norm_factor = 1;
switch(normalize_type)
    case 'TE'
        te_norm_factor = max(TE_reconstruction);
        tm_norm_factor = pol_factor*max(TE_reconstruction);
        if(osa_te_exists), te_osa_norm_factor = max(osa_te_lin); end
        if(osa_tm_exists), tm_osa_norm_factor = max(osa_te_lin); end
    case 'TM'
        te_norm_factor = max(TM_reconstruction)/pol_factor;
        tm_norm_factor = max(TM_reconstruction);
        if(osa_te_exists), te_osa_norm_factor = max(osa_tm_lin); end
        if(osa_tm_exists), tm_osa_norm_factor = max(osa_tm_lin); end
    case 'same'
        te_norm_factor = max(TE_reconstruction);
        tm_norm_factor = max(TM_reconstruction);
        if(osa_te_exists), te_osa_norm_factor = max(osa_te_lin); end
        if(osa_tm_exists), tm_osa_norm_factor = max(osa_tm_lin); end
end

if(do_log)
    TE_plot = 10*log10(TE_reconstruction/te_norm_factor);
    TM_plot = 10*log10(TM_reconstruction/tm_norm_factor);
    if(osa_te_exists), osa_te_plot = 10*log10(osa_te_lin/te_osa_norm_factor); end
    if(osa_tm_exists), osa_tm_plot = 10*log10(osa_tm_lin/tm_osa_norm_factor); end
else
    TE_plot = TE_reconstruction/te_norm_factor;
    TM_plot = TM_reconstruction/tm_norm_factor;
    if(osa_te_exists), osa_te_plot = osa_te_lin/te_osa_norm_factor; end
    if(osa_tm_exists), osa_tm_plot = osa_tm_lin/tm_osa_norm_factor; end
end

nexttile;

plot(plot_lambda, TE_plot, Color = TEcolor);
if(osa_te_exists && show_osa)
    hold on;
    plot(osa_te_lambda, osa_te_plot, 'k--');
    hold off;
end
if(do_log)
    ylim([-40,0]);
    yline(-40:10:0, ':', Color = '#888888');
else
    ylim(lin_range);
end

xlim(lambda_range*1e9);
xticks(these_x_ticks); xticklabels([]);

set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);

if(show_axis_labels)
    xlabel("Wavelength (nm)");
    ylabel("Norm. TE");
    
else
    ylabel(" ");
    yticklabels([]);
end

nexttile;
plot(plot_lambda, TM_plot, Color = TMcolor);
if(osa_tm_exists && show_osa)
    hold on;
    plot(osa_tm_lambda, osa_tm_plot, 'k--');
    hold off;
end
if(do_log)
    ylim([-40,0]);
    yline(-40:10:0, ':', Color = '#888888');
else
    ylim(lin_range);
end
 xlim(lambda_range*1e9);
xticks(these_x_ticks); 
xtickangle(0);
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);

if(show_axis_labels)
    ylabel("Norm. TM");
    xlabel("Wavelength (nm)");
    xticklabels(these_x_labels);
else
    ylabel(" "); % necessary to keep scaling the same
    xlabel(" ");
    xticklabels([]);
    yticklabels([]);
end

%export_fig -clipboard -m2