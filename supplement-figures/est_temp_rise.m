% Load all those FILES
% FTS calibration data
TE_calib_file = 'data\20240922\Q1-0p4-0p4\Q1-0p4-0p4_calib\TE.mat';
TM_calib_file = 'data\20240922\Q1-0p4-0p4\Q1-0p4-0p4_calib\TM.mat';
load(TE_calib_file, 'opl_poly', 'acq_struct', 'start_time', 'end_time');
TE_opl_poly = opl_poly;
[I_inc_out, I_dec_out, ~, heater_time] = ...
            generate_heater_lists(acq_struct);
crop_idxs = (heater_time > start_time) & (heater_time < end_time);
heater_time = heater_time(crop_idxs);
P_inc_out = I_inc_out(crop_idxs).^2*acq_struct.heater_R;
P_dec_out = I_dec_out(crop_idxs).^2*acq_struct.heater_R;
P_difference = P_inc_out - P_dec_out;

clear opl_poly
load(TM_calib_file, 'opl_poly');
TM_opl_poly = opl_poly;

% Ring hot plate data (processed using ring_toc)
TE_ring_file = 'data\20241010\ring_TE_processed.mat';
load(TE_ring_file, 'delta_neffs', 'temps');
TE_delta_neffs = delta_neffs; TE_temps = temps;

TM_ring_file = 'data\20241010\ring_TM_processed.mat';
load(TM_ring_file, 'delta_neffs', 'temps');
TM_delta_neffs = delta_neffs; TM_temps = temps;
clear delta_neffs temps
%% Plot avg neff vs temp for TE vs TM
TEcolor = '#3081D0'; TMcolor = '#B31312';
figure(Units="centimeters",Position=[6 6 6 4.5]); hold on;
plot(TE_temps-min(TE_temps), mean(TE_delta_neffs, 2), Color = TEcolor);
plot(TM_temps-min(TM_temps), mean(TM_delta_neffs, 2), Color = TMcolor);
hold off;
xlabel("Temperature change (K)");
ylabel("n_{eff} change"); legend("TE", "TM", Location = 'northwest');
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
%% plot OPL difference versus power
% Extract approx TOC -> est temperature rise
TE_poly = polyfit(TE_temps, mean(TE_delta_neffs, 2),1);
TM_poly = polyfit(TM_temps, mean(TM_delta_neffs, 2),1);
TE_opls = polyval(TE_opl_poly, heater_time);
TM_opls = polyval(TM_opl_poly, heater_time);
TE_opls = TE_opls-mean(TE_opls);
TM_opls = TM_opls-mean(TM_opls);

figure(Units="centimeters",Position=[6 6 6 4.5]); hold on;
plot(P_difference, TE_opls, Color = TEcolor);
plot(P_difference, TM_opls, Color = TMcolor);
hold off; legend("TE", "TM", Location = 'northwest');
xlabel("Power difference (W)"); ylabel("OPL difference (m)");
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
%% plot approx temperature difference versus power
spec_length = 0.0295; % physical length of spectrometer (m)
TE_temps = TE_opls/(spec_length*TE_poly(1));
TM_temps = TM_opls/(spec_length*TM_poly(1));
figure(Units="centimeters",Position=[6 6 6 4.5]); hold on;
plot(P_difference, TE_temps, Color = TEcolor);
plot(P_difference, TM_temps, Color = TMcolor);
hold off; legend("TE", "TM", Location = 'northwest');
xlabel("Power difference (W)"); ylabel("Temperature difference (K)");
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
te_efficiency = mean(diff(TE_temps)./diff(P_difference))
tm_efficiency = mean(diff(TM_temps)./diff(P_difference))