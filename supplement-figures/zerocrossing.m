clear;
% this file is just a modified version of calibration_processing.m so there
% are some extraneous sections left over (but all still have to be run)

% load in a folder of .mat files 
%dirname = uigetdir();

%  --- FIG S4 --- %
dirname = 'data\20241014\calibrations\TE_calib'; % TE data
which_polarization = 'TE';
% --------------- %

%  --- FIG S3 --- %
% dirname = 'data\20241014\calibrations\TM_calib'; % TM data
% which_polarization = 'TM';
% --------------- %

% Load data from files
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
[~,common_time_array,~] = generate_detector_params(common_acq_struct);
test_names = test_names(2:end); % remove empty string in first spot
num_tests = size(test_results,2);

te_colors = crameri('devon', num_tests+2); % +2 to avoid too-white trace
tm_colors = flipud(crameri('lajolla', num_tests+2));
%% Plot raw data and crop region
if(strcmp(which_polarization, 'TE'))
    colors = te_colors;
else
    colors = tm_colors;
end

% --- FALSE for pane a, TRUE for pane b --- %
show_cropped = true;
% ----------------------------------------- %

% we have to crop data to the region over which the sweep actually occurs
start_time = 0.94; end_time = 10.805; % this is chosen by user as needed
time_crop_idxs = common_time_array > start_time & common_time_array < end_time;
crop_time = common_time_array(time_crop_idxs);
figure("Units", "centimeters", "Position", [6,6,4,10]); 
hold on;
for this_idx = 1:num_tests
    this_data = test_results(:,this_idx);
    this_data = this_data/max(this_data) + this_idx - 1.5;
    this_time = common_time_array;
    if(~show_cropped)
        % if we're not cropping, reduce # of data points to reasonable
        % number so we can still do vector export
        % purposefully DON'T use decimate as this lpf's data
        factor = 1;
        this_time = downsample(this_time,factor);
        this_data = downsample(this_data,factor);
    end
    this_color = colors(this_idx,:);
    plot(this_time, this_data, "DisplayName", "Original", Color = this_color);
    
end
%xline([start_time, end_time]);


xlabel("Time (s)");
ylabel("Norm + shifted power (a.u.)");
if(show_cropped)
    xlim([5,5.2]);
else
    xlim([start_time,end_time]);
end

%title(test_names(test_idx), "Interpreter", "none");
%legend("Location","north");

%% Compile & preprocess data (plot not used)
figure(Units="inches", Position=[3 1 10 5.5]); 
t = tiledlayout(2, 1,'TileSpacing','Compact','Padding','Compact');
%colors = cool(num_tests);
preprocessed_interferograms = cell(0);
N = length(crop_time);
for i = 1:num_tests
    %thisHeater = testResults{i}(:,1);
    this_transmission = test_results(time_crop_idxs,i);
    this_color = colors(i,:);
    this_detrended = this_transmission./movmean(this_transmission, 100) - 1;
    % plot raw detrended interferogram
    nexttile(1);
    plot(crop_time, 2*i+this_detrended, ...
        "DisplayName", test_names(i), ...
        "Color", this_color);
    if(i == 1)
        hold on; xlabel("Time (s)"); ylabel("Detrended power");
    end

    this_fft = abs(fft(this_detrended)).^2;
    % our power vs. data point should be super linear
    this_dt = (max(crop_time) - min(crop_time))./(N-1);
    this_freq_axis = (0:N-1)/(N*this_dt); this_max_n = round(N/2);
    % legacy functionality, only necessary on poor-performing devices
    % find absolute max within this range
    peak_search_range = [0, 1e5]; 
    freqCropIdxs = this_freq_axis > peak_search_range(1) & this_freq_axis < peak_search_range(2);
    [this_max, max_idx] = max(this_fft(freqCropIdxs));
    max_idx = max_idx + find(freqCropIdxs,1) - 1;
    this_peak_freq = this_freq_axis(max_idx);

    % plot FFT
    nexttile(2);
    plot_n = 1:this_max_n;
    semilogy(this_freq_axis(plot_n), this_fft(plot_n), ...
        "DisplayName", test_names(i), ...
        "Color", this_color);
    if(i == 1)
        hold on;
        xlabel("Frequency (Hz)");
        ylabel("Original FFT");
        %xlim([0,0.08]);
        %ylim([1, 1e5]);
    end
    plot(this_freq_axis(max_idx), this_max, '*', ...
        "Color", this_color, ...
        HandleVisibility="off");
    this_preprocessed_interferogram = {[crop_time, this_detrended]};
    preprocessed_interferograms = [preprocessed_interferograms, this_preprocessed_interferogram];
end
nexttile(2);
%xline(peak_search_range, 'k--', HandleVisibility="off");
hold off; %legend("Interpreter", "none", "Location", "best");
%% Zero-crossing OPL extraction
% --- FALSE for pane c, TRUE for pane d --- %
dev_from_poly = false;
% ----------------------------------------- %

figure("Units", "centimeters", "Position", [6,6,6,5]); 
hold on;
offset = 1;
clear result_polyfits
polyfit_order = 7;

for i = 1:num_tests
    this_time = preprocessed_interferograms{i}(:,1);
    this_normalized = preprocessed_interferograms{i}(:,2);
    this_lambda = test_lambdas(i);
    this_name = sprintf("%.0f nm", this_lambda*1e9);

    [~,~,zero_crossing_indices] = zerocrossrate(this_normalized);
    zero_crossing_times = this_time(zero_crossing_indices);
    zero_crossing_opls = this_lambda/2*(1:length(zero_crossing_times));
    [this_polyfit, S] = polyfit(zero_crossing_times, zero_crossing_opls, polyfit_order);
    % R_squared = 1 - (S.normr/norm(zeroCrossingOPLs - mean(zeroCrossingOPLs)))^2
    % subtract mean for plotting to center at zero
    time_plot = zero_crossing_opls;
    
    if(dev_from_poly)
        y_plot = 1e6*(polyval(this_polyfit, zero_crossing_times) - zero_crossing_opls');
        ylabel("Deviation (um)");
        %title("Deviation of OPL vs. T from polynomial fit")
    else
        y_plot = 1e6*(zero_crossing_opls - mean(zero_crossing_opls));
        ylabel("OPL difference (um)");
        %title("OPL Difference vs. Power")
    end
    dec_factor = 5; % gotta do this to make vector export reasonable
    plot(downsample(zero_crossing_times, dec_factor), ...
        downsample(y_plot, dec_factor), '.', ...
            "Color", colors(i,:), "DisplayName", this_name);
    result_polyfits(i,:) = this_polyfit;
end
hold off;
%legend("Interpreter","none", Location="best");
xlabel("Time (s)"); xlim('tight')
%% plot g(t) at different lambdas (plot not used)
center_index = 3;
% TODO using P_start/P_end here is a little sus
%P_start = -500; P_end = 500;
%powerArray = linspace(P_start, P_end);
center_deriv = polyval(polyder(result_polyfits(center_index,:)),crop_time);
figure; hold on;
mean_g = zeros(1,length(test_lambdas));
for i = 1:length(test_lambdas)
    this_name = sprintf("%1.3e m", test_lambdas(i));
    this_g = polyval(polyder(result_polyfits(i,:)),crop_time)./center_deriv;
    plot(crop_time, this_g, "Color", colors(i,:), "DisplayName", this_name);
    mean_g(i) = mean(this_g); % for now, use mean g instead of center G
end
hold off;
xlabel("Power (mW)"); ylabel("g factor");
legend;
%% polyfits
c = 299792458; % assuming we're using lambda in free space!
test_freqs = c./test_lambdas;
nu_min = min(test_freqs);
nu_max = max(test_freqs);
nu_plot = linspace(nu_min, nu_max);
% OPL fit at nu0
fit_index = 2;
nu0 = test_freqs(fit_index);
opl_poly = result_polyfits(center_index,:);
% do a polyfit on g, in frequency units as that makes most sense during
% extraction
g_poly_order = 3;
[g_poly, ~, g_mu] = polyfit(test_freqs, mean_g, g_poly_order);
figure("Units", "centimeters", "Position", [6,6,6,5]);  hold on;
plot(1e-12*test_freqs, mean_g, 'ro');
plot(1e-12*nu_plot, polyval(g_poly, nu_plot, [], g_mu), 'b');
plot(1e-12*test_freqs, ones(size(test_freqs)), 'k:');
hold off;
xlabel("Frequency (THz)"); ylabel("g(ν) factor");
xlim('tight');
legend('Data', 'Polynomial Fit', location = 'best');
%% save calibration data
acq_struct = common_acq_struct;
[output_filename, output_path] = uiputfile('*', 'Select location to save data:');
if(output_filename)
    save(strcat(output_path,output_filename), ...
        'acq_struct', ... % acquisition settings
        'opl_poly', 'start_time', 'end_time', ... % time domain scaling + range
        'g_poly', 'g_mu', 'nu_min', 'nu_max'); % freq domain scaling
else
    disp("File save cancelled");
end