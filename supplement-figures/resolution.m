%% Load data
clear;
% calibration files
TE_calib_file = "data\20241014\calibrations\calib\TE.mat";
TM_calib_file = "data\20241014\calibrations\calib\TM.mat";
% directories of calibration source data, i.e. interferograms of various
% single lasers across BW
TE_data_dir = "data\20241014\calibrations\TE_calib";
TM_data_dir = "data\20241014\calibrations\TM_calib";
%% Measure from reconstructions
% experimental width measurement
[TE_lambdas, TE_widths] = measure_reconstructed_fwhms(TE_data_dir, TE_calib_file);
[TM_lambdas, TM_widths] = measure_reconstructed_fwhms(TM_data_dir, TM_calib_file);
c0 = physconst('LightSpeed');
lambda_eval = linspace(1.48e-6, 1.63e-6);
nu_eval = c0./lambda_eval;
%% Plot
TE_theory = calculate_theoretical_resolution(TE_calib_file, nu_eval);
TM_theory = calculate_theoretical_resolution(TM_calib_file, nu_eval);
TEcolor = '#3081D0'; TMcolor = '#B31312';
figsize = [12 8];
figure(Units = "centimeters", Position=[3 3 figsize], ...
    PaperUnits="centimeters", PaperSize=figsize, PaperPositionMode="auto"); 
tiledlayout(2,1);
nexttile;
plot(TE_lambdas, TE_widths, 'o', Color = TEcolor); hold on;
plot(1e9*lambda_eval, 1e9*TE_theory, Color = TEcolor); hold off;
legend('Experiment', 'Theory', 'Location','southeast');
ylabel("TE FWHM (nm)");
nexttile;
plot(TM_lambdas, TM_widths, 'o', Color = TMcolor); hold on;
plot(1e9*lambda_eval, 1e9*TM_theory, Color = TMcolor); hold off;
legend('Experiment', 'Theory', 'Location','southeast');
xlabel("Wavelength (nm)");
ylabel("TM FWHM (nm)");



