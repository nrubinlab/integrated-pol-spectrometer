%% Plots showing FFTs of calibration data (handful of lasers across BW), and nonlinearity correction
clear;

te_dirname = 'data\20241014\calibrations\TE_calib';
tm_dirname = 'data\20241014\calibrations\TM_calib';
processed_calib_dirname = 'data\20241014\calibrations\calib';

% Load data from files
[te_results, te_names, te_lambdas, time_array, num_tests] = read_data_from_dir(te_dirname);
[tm_results, tm_names, tm_lambdas, ~, ~] = read_data_from_dir(tm_dirname);
%% Plot OPL versus time, and again with linear term subtracted off
TE_calib_file = fullfile(processed_calib_dirname, 'TE.mat');
load(TE_calib_file, ...
        'acq_struct', ... % power settings
        'opl_poly', 'start_time', 'end_time'); % time domain scaling + range
time_crop_idxs = time_array > start_time & time_array < end_time;
this_time = time_array(time_crop_idxs);


nonlinear_opl_index = 3; % choose nth index of TE data for nonlinear OPL vis
this_data = te_results(time_crop_idxs,nonlinear_opl_index);
this_lambda = te_lambdas(nonlinear_opl_index);
this_normalized = highpass(this_data, 0.05);
[~,~,zero_crossing_indices] = zerocrossrate(this_normalized);
zero_crossing_times = this_time(zero_crossing_indices);
zero_crossing_opls = this_lambda/2*(1:length(zero_crossing_times));

opls_lin_poly = polyfit(zero_crossing_times, zero_crossing_opls,1);
opls_lin = polyval(opls_lin_poly, zero_crossing_times);

% reduce # of data points for more lightweight vector plot
dec_factor = 10;
centered_opl = zero_crossing_opls - mean(zero_crossing_opls);
plot_opl = decimate(centered_opl, dec_factor);
zeroed_time = zero_crossing_times - min(zero_crossing_times);
plot_time = decimate(zeroed_time, dec_factor);
dev_from_lin = zero_crossing_opls - opls_lin';
plot_dev = decimate(dev_from_lin, dec_factor);

figure(Units="centimeters", Position=[5 5 5.5 4.5]); 
tiledlayout(2, 1,'TileSpacing','Compact','Padding','Compact');
nexttile(1);
plot(plot_time, 1e3*plot_opl, 'k.');
xticklabels([]);
ylabel("ΔOPL (mm)");
set(gca, 'fontsize', 6, 'ticklength', [0.03, 0.03]);

nexttile(2);
plot(plot_time, 1e6*plot_dev, 'k.');
xlabel("Time (s)");
ylabel("ΔOPL nonlin. (um)");
set(gca, 'fontsize', 6, 'ticklength', [0.03, 0.03]);

%% FFT plots, corrected & uncorrected
% --- FALSE for panes d/e, TRUE for panes g/h  --- %
plot_corrected = false;
% ------------------------------------------------ %

fmin = 25;
fmax = 105;
%fmin = 300;
%fmax = 1000;
ymax = 3;
ymin = 1e-3;
log_ticks = [1e-2,1e-1,1];

figure(Units="centimeters", Position=[5 5 7 4.8]); 
t = tiledlayout(2, 1,'TileSpacing','Compact','Padding','Compact');
te_colors = crameri('devon', num_tests+2); % +2 to avoid too-white trace
tm_colors = flipud(crameri('lajolla', num_tests+2));

% settings for uncorrected FFT
crop_time = time_array(time_crop_idxs);
zeropadfactor = 1;
N = length(crop_time)*zeropadfactor;
this_dt = (time_array(2)-time_array(1));
this_freq_axis = (0:N-1)/(N*this_dt); this_max_n = round(N/2);

% settings for corrected FFT
t_corrected = polyval(opl_poly, crop_time);
t_center = median(crop_time);
t_slope_at_zero = polyval(polyder(opl_poly), t_center);
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

function [test_results, test_names, test_lambdas, time_array, num_tests] = read_data_from_dir(dirname)
    file_list = dir(fullfile(dirname,'*.mat'));
    test_lambdas = [];
    test_results = [];
    test_names = strings();
    if(isempty(file_list))
        error("File list length == 0! Either %s does not exist or it contains no files.", dirname);
    end
    for thisIdx = 1:length(file_list)
        thisFile = file_list(thisIdx);
        thisCompleteFile = fullfile(dirname, thisFile.name);
        if(~isfile(thisCompleteFile))
            error("File %s does not exist!", thisCompleteFile)
        end
        load(thisCompleteFile, 'acq_struct', 'channel1', 'lambda'); 
        if(thisIdx == 1)
            common_acq_struct = acq_struct;
        else
            if(~isequal(common_acq_struct, acq_struct))
                error("Time arrays do not match for file %s", thisFile.name);
            end
        end
        % TODO check channel1 shape before doing this
        test_results = [test_results, channel1'];
        test_names = [test_names, thisFile.name];
        test_lambdas = [test_lambdas, lambda];
        clear acq_struct channel1 lambda
    end
    [~,time_array,~] = generate_detector_params(common_acq_struct);
    test_names = test_names(2:end); % remove empty string in first spot
    num_tests = size(test_results,2);
end