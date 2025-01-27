%% Plots showing interferogram and FFT for 45 deg input laser
clear;
data_45 = 'data\20240920\Q1_0p4_0p4\2000mW_45degs\2000mW_45degs_1540.mat';

load(data_45, "acq_struct", "channel1");
channel1_uw = 1e6*channel1;
[~,detector_time,~] = generate_detector_params(acq_struct);
[I_inc_out, I_dec_out, ~, heater_time] = ...
            generate_heater_lists(acq_struct);
figure(Units = "centimeters", Position = [5 5 6.4 3.5]);
hold on;
yyaxis("left");
ylabel("Optical Power (μW)", Color = 'black');
ylim([-0.1, 1.6]*max(channel1_uw));

plot(detector_time, channel1_uw, 'k');

yyaxis("right");
ylabel("Heater Power (W)");
inc_color = '#e96900';
dec_color = '#429130';
P_inc_out = I_inc_out.^2*acq_struct.heater_R;
P_dec_out = I_dec_out.^2*acq_struct.heater_R;
ylim([-0.1, 1.6]*max(P_inc_out));

plot(heater_time, P_inc_out, '-', LineWidth=1.5);%, Color = inc_color);
plot(heater_time, P_dec_out, '--', LineWidth=1.5);%, Color = dec_color);
hold off;
xlabel("Time (s)");
xlim([0.4,1.6]);
ax = gca;
set(ax, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
ax.YAxis(1).Color = 'k';
%ax.YAxis(2).Color = 'k';
%% zoom-in of interferogram
figure(Units = "centimeters", Position = [5 5 2 2]);
plot(detector_time, channel1_uw, 'k');
x_span = 0.03; x_center = 1; x_range = x_center + 0.5*[-x_span, x_span];
xlim(x_range); xticks([]); yticks([]);
%% fft of interferogram
N = length(channel1_uw);
preprocessed_interferogram = channel1_uw.*hamming(N)';

this_fft = abs(fft(channel1_uw.*hamming(N)'));
this_fft = this_fft/max(this_fft);
plot_idx = 1:round(length(this_fft)/2);
dt = detector_time(2) - detector_time(1);

fft_freq = (0:(N-1))/(N*dt);
figure(Units = "centimeters", Position = [5 5 6.4 3.5]);
semilogy(fft_freq(plot_idx), this_fft(plot_idx), 'k');
xlim([0,400]);
ylim([1e-3,1]);
yticks([1e-3,1e-2,1e-1,1]);
xlabel("Freqency (Hz)");
ylabel("Noramlized FFT");
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);