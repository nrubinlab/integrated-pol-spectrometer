% Perform analysis of calibration datasets to generate calibration files
% this must be run once on each polarization's dataset
% the contents of the generated calibration file are simply polynomial fits
% of the time and frequency rescaling needed in the nufft spectrum
% extraction technique

clear;
% load in a folder of .mat files 
dirname = uigetdir();
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
    test_names = [test_names, thisFile.name];
    test_lambdas = [test_lambdas, lambda];
    clear acq_struct channel1 lambda
end
[~,common_time_array,~] = generate_detector_params(common_acq_struct);
test_names = test_names(2:end); % remove empty string in first spot
num_tests = size(test_results,2);
%% plot heater powers and recorded power for one data point
test_idx = 3;
% the start/end times need to be adjusted here to select only the portion
% of the linear ramp that is far enough away from start/end to be quite
% linear

%start_time = 0.94; end_time = 10.805; 
start_time = 0.94; end_time = 1.9;
time_crop_idxs = common_time_array > start_time & common_time_array < end_time;
crop_time = common_time_array(time_crop_idxs);
figure("Units", "inches", "Position", [3,2,4,5]); 
tiledlayout(2,1);
nexttile;
hold on;
plot(common_time_array, 1000*test_results(:,test_idx), "DisplayName", "Original");
plot(crop_time, 1000*test_results(time_crop_idxs,test_idx), "DisplayName", "Cropped");
xlabel("Time (s)");
ylabel("Optical Power (mW)");
title(test_names(test_idx), "Interpreter", "none");
legend("Location","north");
nexttile;
semilogy(abs(fft(1000*test_results(time_crop_idxs,test_idx))));
xlabel("FFT index"); ylabel("FFT magnitude (log)");
%% Compile data
% This section must be run to compile the data into a convenient data
% structure. Also in this section is a bandpass filter functionality.

% this is a deprecated functionality not used in the final results in the
% paper. In early spectrometer prototypes, the sidewall roughness was
% extremely high and the interferograms were very messy in the time domain,
% but a clear spike was visible in the frequency domain. As such, a
% narrow bandpass filter was used to select only the sinusoidal component
% (with slight nonlinearities) before performing zero crossing OPL
% extraction
peak_search_range = [0, 1e5]; % units of Hz, only necessary on poor-performing devices
do_bandpass = false;
bandpass_width = 200;
figure(Units="inches", Position=[3 1 10 5.5]); 
if do_bandpass, num_plot_tiles = 4; else, num_plot_tiles = 2; end
t = tiledlayout(num_plot_tiles, 1,'TileSpacing','Compact','Padding','Compact');
colors = cool(num_tests);
preprocessed_interferograms = cell(0);
N = length(crop_time);
for i = 1:num_tests
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
    % find absolute max within this range
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
    end
    plot(this_freq_axis(max_idx), this_max, '*', ...
        "Color", this_color, ...
        HandleVisibility="off");
    if(do_bandpass)
        % perform filtering
        this_filtered = bandpass(this_detrended, ...
            [this_peak_freq-bandpass_width/2, this_peak_freq+bandpass_width/2], ...
            1/this_dt);
        
        % plot filtered FFT
        nexttile(3);
        filtered_fft = abs(fft(this_filtered)).^2;
        semilogy(this_freq_axis(plot_n), filtered_fft(plot_n), ...
            "DisplayName", test_names(i), ...
            "Color", this_color);
        if(i == 1)
            hold on;
            xlabel("Frequency (Hz)");
            ylabel("Filtered FFT");
            
        end
    
        % plot filtered interferograms
        nexttile(4);
        plot(crop_time, this_filtered./(2*max(this_filtered)) + i - 1, ...
            "DisplayName", test_names(i), ...
            "Color", this_color);
        if(i == 1)
            hold on;
            xlabel("Time (s)");
            ylabel("Filtered interferogram");
        end
        this_preprocessed_interferogram = {[crop_time, this_filtered]};
    else
        this_preprocessed_interferogram = {[crop_time, this_detrended]};
        
    end
    preprocessed_interferograms = [preprocessed_interferograms, this_preprocessed_interferogram];
end
nexttile(2);
%xline(peak_search_range, 'k--', HandleVisibility="off");
hold off; %legend("Interpreter", "none", "Location", "best");
%% Zero-crossing OPL extraction
figure("Units", "inches", "Position", [3,3,5,2.5]); 
hold on;
offset = 1;
colors = cool(num_tests);
clear result_polyfits

polyfit_order = 7;

% all these are just plot parameters and don't impact calibration
% typically various plots are viewed here before proceeding to check that
% the calibration result is robust
plot_interferogram = false;
plot_envelope = false;
plot_normalized = true;
dev_from_linear = true;
dev_from_poly = true;
for i = 1:num_tests
    this_time = preprocessed_interferograms{i}(:,1);
    this_data = preprocessed_interferograms{i}(:,2);
    this_lambda = test_lambdas(i);
    this_name = test_names(i);
    this_normalized = this_data; % legacy variable name
    if(plot_interferogram)
        if(plot_normalized)
            this_plot_data = this_normalized + offset*(i-1);
        else
            this_plot_data = 1000*this_data;
        end
            
        plot(this_time, ...
             this_plot_data,...
             "Color", colors(i,:), "DisplayName", this_name, "LineWidth",0.25)
        if(plot_envelope)
            plot(this_time, ...
             1000*this_upper,...
             "Color", colors(i,:), "HandleVisibility", "off");
            plot(this_time, ...
             1000*this_lower,...
             "Color", colors(i,:), "HandleVisibility", "off");
        end
    end
    [~,~,zero_crossing_indices] = zerocrossrate(this_normalized);
    zero_crossing_times = this_time(zero_crossing_indices);
    zero_crossing_opls = this_lambda/2*(1:length(zero_crossing_times));
    [this_polyfit, S] = polyfit(zero_crossing_times, zero_crossing_opls, polyfit_order);
    % uncomment this line to see comically high R values
    % R_squared = 1 - (S.normr/norm(zeroCrossingOPLs - mean(zeroCrossingOPLs)))^2
    
    if(~plot_interferogram)
        if(dev_from_poly)
            plot(zero_crossing_times, ...
                1e6*(polyval(this_polyfit, zero_crossing_times) - zero_crossing_opls'), '*', ...
                "Color", colors(i,:), "DisplayName", this_name);
            ylabel("Deviation (um)");
            title("Deviation of OPL vs. T from polynomial fit")
        elseif(dev_from_linear)
            linPolyfit = polyfit(zero_crossing_times, zero_crossing_opls, 1);
            plot(zero_crossing_times, ...
                1e6*(polyval(linPolyfit, zero_crossing_times) - zero_crossing_opls'), '*', ...
                "Color", colors(i,:), "DisplayName", this_name);
            ylabel("Deviation (um)");
            title("Deviation of OPL vs. T from linear");
        else
            plot(zero_crossing_times, ...
                1e6*(zero_crossing_opls - this_polyfit(end)), '*', ...
                "Color", colors(i,:), "DisplayName", this_name);
            plot(this_time, 1e6*(polyval(this_polyfit, this_time) - this_polyfit(end)), ...
                        Color = colors(i,:), DisplayName = this_name + " fit")
            ylabel("OPL difference (um)");
            title("OPL Difference vs. Power")
        end
    end
    result_polyfits(i,:) = this_polyfit;
end
hold off;
legend("Interpreter","none");
xlabel("Time (s)");
if(plot_interferogram)
    if(plot_normalized)
        ylabel("Normalized Power");
    else
        ylabel("Optical Power (mW)");
    end
end
%% plot g(t) at different lambdas
% the polyfit of this g(t) vs lambda is the frequency rescaling that will
% be used in final calibration file
center_index = 3;
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
%% polyfits for final export 
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
g_poly_order = 3; % this is all that's needed in nearly all cases!
[g_poly, ~, g_mu] = polyfit(test_freqs, mean_g, g_poly_order);
figure; hold on;
plot(test_freqs, mean_g, 'o');
plot(nu_plot, polyval(g_poly, nu_plot, [], g_mu));
hold off;
xlabel("Frequency"); ylabel("g(nu)");
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