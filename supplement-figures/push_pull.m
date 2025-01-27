%% Plots showing impact of push-pull tuning on temp-domain nonlinearity
clear;
% --- TRUE for first column, FALSE for second --- %
balanced = true;
% ----------------------------------------------- %

if(balanced)
    % with push/pull:
    te_dirname = 'data\20241010\TE_6W_balanced';
    tm_dirname = 'data\20241010\TM_6W_balanced';
    processed_calib_dirname = 'data\20241010\bal_calibs';
else
    % without push/pull
    te_dirname = 'data\20241010\TE_6W_unbalanced';
    tm_dirname = 'data\20241010\TM_6W_unbalanced';
    processed_calib_dirname = 'data\20241010\unbal_calibs';
end

% Load data from files
[te_results, te_names, te_lambdas, time_array, num_tests] = read_data_from_dir(te_dirname);
[tm_results, tm_names, tm_lambdas] = read_data_from_dir(tm_dirname);

te_colors = crameri('devon', num_tests+2); % +2 to avoid too-white trace
tm_colors = flipud(crameri('lajolla', num_tests+2));
%% Plot interferogram and heater tuning

selected_index = 3;
TE_calib_file = fullfile(processed_calib_dirname, 'TE.mat');
load(TE_calib_file, 'acq_struct', 'start_time', 'end_time', 'opl_poly'); 
% BE CAREFUL!!!! must un-comment the correct function here to generate
% correct heater vs time curve
if(balanced)
    [I_inc_out, I_dec_out, supply_interval, heater_time] = ...
                generate_heater_lists(acq_struct);
else
    [I_inc_out, I_dec_out, supply_interval, heater_time] = ...
        generate_heater_lists_unbalanced(acq_struct);
end

P_inc = acq_struct.heater_R*I_inc_out.^2;
P_dec = acq_struct.heater_R*I_dec_out.^2;
figure(Units="centimeters", Position=[3 3 10 3]);
yyaxis left
plot(time_array-start_time, 1e6*te_results(:,selected_index), 'k');
ylabel("Optical Power (uW)");
yyaxis right
hold on;
plot(heater_time-start_time, P_inc, LineWidth=2);
plot(heater_time-start_time, P_dec, LineWidth=2);
ylabel("Heater Power (W)");
hold off;
xlabel("Time (s)");
xlim([0, end_time-start_time])
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
%% Plot OPL versus time, and also scanning speed for that interferogram
% some steps in here are special because unbalanced tuning is not fit well
% by a polynomial fit

time_crop_idxs = time_array > start_time & time_array < end_time;
crop_time = time_array(time_crop_idxs);

this_data = te_results(time_crop_idxs,selected_index);
interp_factor = 16;
this_lambda = te_lambdas(selected_index);

[zero_crossing_times, zero_crossing_opls] = opls_from_interferogram_interp( ...
    crop_time, this_data, this_lambda, interp_factor);

% numerical derivative
opls_derivative = gradient(zero_crossing_opls)./gradient(zero_crossing_times);
% smooth data
%opls_derivative = lowpass(opls_derivative, 0.1);

opls_lin_poly = polyfit(zero_crossing_times, zero_crossing_opls,1);
opls_lin = polyval(opls_lin_poly, zero_crossing_times);

% reduce # of data points for more lightweight vector plot
dec_factor = 10;
centered_opl = zero_crossing_opls - mean(zero_crossing_opls);
plot_opl = downsample(centered_opl, dec_factor);
zeroed_time = zero_crossing_times - min(zero_crossing_times);
plot_time = downsample(zeroed_time, dec_factor);
dev_from_lin = zero_crossing_opls - opls_lin';
plot_dev = downsample(dev_from_lin, dec_factor);
plot_deriv = downsample(opls_derivative, dec_factor);

figure(Units="centimeters", Position=[5 5 10 4.5]); 
tiledlayout(2, 1,'TileSpacing','Compact','Padding','Compact');
nexttile(1);
plot(plot_time, 1e3*plot_opl, 'k.');
xticklabels([]); xlim('tight');
ylabel("ΔOPL (mm)");
set(gca, 'fontsize', 6, 'ticklength', [0.03, 0.03]);

nexttile(2);
%plot(plot_time, 1e6*plot_dev, 'k.');
plot(plot_time, 1e3*plot_deriv, 'k.');
xlabel("Time (s)"); xlim('tight');
ylabel("Sweep speed (mm/s))");
set(gca, 'fontsize', 6, 'ticklength', [0.03, 0.03]);
ylim([1,2]);
%% Calculate OPL derivatives without polyfit (plot not used)
num_wavs = length(te_lambdas);
te_derivs = cell(num_wavs,1);
te_derivs_interp = cell(num_wavs,1);
tm_derivs = cell(num_wavs,1);
tm_derivs_interp = cell(num_wavs,1);
figure; hold on;
interp_factor = 50;

for i = 1:num_wavs
    % TE
    [te_times, te_opls] = opls_from_interferogram_interp( ...
        crop_time, te_results(time_crop_idxs,i), te_lambdas(i), interp_factor);
    te_derivs{i} = [te_times', ...
        gradient(te_opls')./gradient(te_times')];
    te_derivs_interp{i} = interp1(te_derivs{i}(:,1), te_derivs{i}(:,2), crop_time, 'linear', 'extrap');
    plot(te_derivs{i}(:,1), ...
        te_derivs{i}(:,2), ...
        '.-', ...
        "Color", te_colors(i,:))
    % TM
    [tm_times, tm_opls] = opls_from_interferogram_interp( ...
        crop_time, tm_results(time_crop_idxs,i), tm_lambdas(i), interp_factor);
    tm_derivs{i} = [tm_times', ...
        gradient(tm_opls')./gradient(tm_times')];
    tm_derivs_interp{i} = interp1(tm_derivs{i}(:,1), tm_derivs{i}(:,2), crop_time, 'linear', 'extrap');

    plot(tm_derivs{i}(:,1), ...
        tm_derivs{i}(:,2), ...
        '.-', ...
        "Color", tm_colors(i,:))
end

%% g(nu)
center_index = 2;
selected_index = 3;
movmean_length = 200; % use a rolling average to suppress noise
figure(Units="centimeters", Position=[5 5 10 3]); 
tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
nexttile(1); hold on; xlabel("Time"); ylabel(sprintf("g_{TE}, λ = %.0f nm", 1e9*te_lambdas(selected_index))); %ylim([0.97,1.03]);

this_g_te = movmean(te_derivs_interp{selected_index},movmean_length, 'Endpoints','fill')...
        ./movmean(te_derivs_interp{center_index},movmean_length, 'Endpoints','fill');
nexttile(1);
plot(crop_time - min(crop_time), this_g_te, ...
    "Color", te_colors(selected_index,:));
ylim('tight');
set(gca, 'fontsize', 6, 'ticklength', [0.03, 0.03]);

nexttile(2); hold on; xlabel("Time"); ylabel(sprintf("g_{TM}, λ = %.0f nm", 1e9*te_lambdas(selected_index)))
this_g_tm = movmean(tm_derivs_interp{selected_index},movmean_length, 'Endpoints','fill')...
    ./movmean(tm_derivs_interp{center_index},movmean_length, 'Endpoints','fill');
nexttile(2);
plot(crop_time - min(crop_time), this_g_tm, ...
    "Color", tm_colors(selected_index,:));
ylim('tight');
set(gca, 'fontsize', 6, 'ticklength', [0.03, 0.03]);
%% h(nu)
center_index = 2;
figure(Units="centimeters", Position=[5 5 10 3]); 
tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
h = movmean(te_derivs_interp{center_index},movmean_length, 'Endpoints','fill')...
    ./movmean(tm_derivs_interp{center_index},movmean_length, 'Endpoints','fill');
nexttile(1);
plot(crop_time - min(crop_time), h, 'k-');
xlim('tight'); xlabel("Time (s)");
ylabel(sprintf("h_{TM}, λ = %.0f nm", 1e9*te_lambdas(center_index)));
set(gca, 'fontsize', 6, 'ticklength', [0.03, 0.03]);


%% Polynomial fit to OPL data (plot not used)
% Note that unbalanced sweep is not fit well by a polynomial because of the
% sharp edge at the middle of the sweep. Fixing this by using some sort of
% interpolation instead of a polynomial is something that we could do, but
% it's not worth spending time on for this supplemental plot

time_crop_idxs = time_array > start_time & time_array < end_time;
crop_time = time_array(time_crop_idxs);

fit_order = 15; % 7 for balanced, 15 for unbalanced
num_wavs = length(te_lambdas);
te_polys = zeros(num_wavs, fit_order+1); te_mus = zeros(num_wavs, 2);
tm_polys = zeros(num_wavs, fit_order+1); tm_mus = zeros(num_wavs, 2);

figure; hold on;
for i = 1:num_wavs
    % TE
    this_data = te_results(time_crop_idxs,i);
    this_lambda = te_lambdas(i);
    this_normalized = highpass(this_data, 0.05);
    [~,~,zero_crossing_indices] = zerocrossrate(this_normalized);
    zero_crossing_times = crop_time(zero_crossing_indices);
    zero_crossing_opls = this_lambda/2*(1:length(zero_crossing_times));
    [te_polys(i,:), ~, te_mus(i,:)] = polyfit(zero_crossing_times, zero_crossing_opls,fit_order);

    plot(zero_crossing_times, ...
        zero_crossing_opls' - polyval(te_polys(i,:), zero_crossing_times, [], te_mus(i,:)), ...
        'o')

    % TM
    this_data = tm_results(time_crop_idxs,i);
    this_lambda = tm_lambdas(i);
    this_normalized = highpass(this_data, 0.05);
    [~,~,zero_crossing_indices] = zerocrossrate(this_normalized);
    zero_crossing_times = crop_time(zero_crossing_indices);
    zero_crossing_opls = this_lambda/2*(1:length(zero_crossing_times));
    [tm_polys(i,:), ~, tm_mus(i,:)] = polyfit(zero_crossing_times, zero_crossing_opls,fit_order);

%     plot(zero_crossing_times, ...
%         zero_crossing_opls' - polyval(tm_polys(i,:), zero_crossing_times, [], tm_mus(i,:)), ...
%         'o')
end
%% FFT plots, corrected & uncorrected
plot_corrected = true;

fmin = 200;
fmax = 1000;
ymax = 3;
ymin = 1e-3;
log_ticks = [1e-2,1e-1,1];

figure(Units="centimeters", Position=[5 5 10 4.8]); 
%figure(Units="centimeters", Position=[5 5 7 7]); 
t = tiledlayout(2, 1,'TileSpacing','Compact','Padding','Compact');
% te_colors = crameri('batlow', num_tests);
% tm_colors = crameri('batlow', num_tests);

% settings for uncorrected FFT
zeropadfactor = 1;
N = length(crop_time)*zeropadfactor;
this_dt = (time_array(2)-time_array(1));
this_freq_axis = (0:N-1)/(N*this_dt); this_max_n = round(N/2);

% settings for corrected FFT
%t_corrected = polyval(opl_poly, crop_time);
t_corrected = polyval(te_polys(center_index,:), crop_time, [], te_mus(center_index,:));
t_center = median(crop_time);
%t_slope_at_zero = polyval(polyder(opl_poly), t_center);
t_slope_at_zero = polyder_with_mu(te_polys(center_index,:), ...
                                t_center, te_mus(center_index,:));
t_corrected = t_corrected/t_slope_at_zero;
% since we aren't doing freq scaling here, freq units here are ill-defined
% after nonlinearity correction. However, since the nonlinearity is small,
% scaling using dt like this gets it close enough for a visualization
f_sample_corrected = linspace(0,1,5000)/this_dt;

for i = 1:num_tests
    te_transmission = te_results(time_crop_idxs,i);
    te_detrended = te_transmission./movmean(te_transmission, 100) - 1;
    te_detrended = hann(length(te_detrended)).*te_detrended;

    tm_transmission = tm_results(time_crop_idxs,i);
    tm_detrended = tm_transmission./movmean(tm_transmission, 100) - 1;
    tm_detrended = hann(length(tm_detrended)).*tm_detrended;

    if(~plot_corrected)
        te_fft = abs(fft(te_detrended,N));
        tm_fft = abs(fft(tm_detrended,N));
    else
        te_fft = abs(nufft(te_detrended, t_corrected, f_sample_corrected));
        tm_fft = abs(nufft(tm_detrended, t_corrected, f_sample_corrected));
    end

    te_fft = te_fft/max(te_fft);
    te_label = sprintf("%1.0f nm", 1e9*te_lambdas(i));
    
    tm_fft = tm_fft/max(tm_fft);
    tm_label = sprintf("%1.0f nm", 1e9*tm_lambdas(i));

    plot_n = 1:this_max_n;
    if(~plot_corrected)
        plot_f = this_freq_axis(plot_n);
        plot_te = te_fft(plot_n);
        plot_tm = tm_fft(plot_n);
    else
        plot_f = f_sample_corrected;
        plot_te = te_fft;
        plot_tm = tm_fft;
    end

    % plot FFT

    nexttile(1);
    
    semilogy(plot_f, plot_te, ...
        "DisplayName", te_label, ...
        "Color", te_colors(i,:));
    nexttile(2);
    semilogy(plot_f, plot_tm, ...
        "DisplayName", tm_label, ...
        "Color", tm_colors(i,:));
    if(i == 1)

        nexttile(1);
        hold on;
        set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
        xlim([fmin,fmax]);
        ylim([ymin, ymax]);
        yticks(log_ticks);
        ylabel("FFT (a.u)");
        xticklabels([]);
        nexttile(2);
        hold on;
        xlabel("Frequency (Hz)");
        set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
        xlim([fmin,fmax]);
        ylim([ymin, ymax]);
        yticks(log_ticks);
        ylabel("FFT (a.u)");
    end
end
nexttile(2);
hold off;
%legend("Interpreter", "none", "Location", "best");



function out = polyder_with_mu(p, x, mu)
    out = polyval(polyder(p), x, [], mu)/mu(2);
end

function [test_results, test_names, test_lambdas, time_array, num_tests, common_acq_struct] = read_data_from_dir(dirname)
    file_list = dir(fullfile(dirname,'*.mat'));
    test_lambdas = [];
    test_results = [];
    test_names = strings();
    for thisIdx = 1:length(file_list)
        thisFile = file_list(thisIdx);
        load(fullfile(dirname, thisFile.name), 'acq_struct', 'channel1', 'lambda'); 
        if(thisIdx == 1)
            common_acq_struct = acq_struct;
        else
            if(~isequal(common_acq_struct, acq_struct))
                error("Time arrays do not match for file %s", thisFile.name);
            end
        end
        % TODO check channel1 shape before doing this
        test_results = [test_results, channel1'];
        %testPowers = [testPowers, thesePowers];
        test_names = [test_names, thisFile.name];
        test_lambdas = [test_lambdas, lambda];
        clear acq_struct channel1 lambda
    end
    [~,time_array,~] = generate_detector_params(common_acq_struct);
    test_names = test_names(2:end); % remove empty string in first spot
    num_tests = size(test_results,2);
end

function [times, opls] = opls_from_interferogram_interp(this_time, this_data, this_lambda, interp_factor)
    % we know we sampled at nyquist rate, so let's upsample to get more
    % clean derivatives in our plots (this is not necessary for actual spectrum
    % extractions, as taking the derivative of a polyfit yields good results)
    this_data_interp = interpft(this_data, interp_factor*length(this_data));
    crop_time_interp = linspace(this_time(1), this_time(end), interp_factor*length(this_time));
    %this_normalized = highpass(this_data_interp, 0.05/interp_factor);
    movmean_window_relative = 0.1;
    movmean_window = round(movmean_window_relative*length(this_data_interp));
    this_normalized = this_data_interp./movmean(this_data_interp, movmean_window) - 1;
    [~,~,zero_crossing_indices] = zerocrossrate(this_normalized);
    times = crop_time_interp(zero_crossing_indices);
    opls = this_lambda/2*(1:length(times));
    % delete the first/last few as they're usually bad
    remove_num = 1;
    times = times(remove_num+1:end-remove_num);
    opls = opls(remove_num+1:end-remove_num);
end