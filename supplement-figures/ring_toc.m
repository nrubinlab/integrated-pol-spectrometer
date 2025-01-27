% this file tracks resonances and saves a result used for temp plots
% in est_temp_rise.m
% data is in data/20241010/Q2_ring_TE and Q2_ring_TM
[file, location] = uigetfile('.mat', 'Shift-click all files', 'MultiSelect', 'on');
if(iscell(file))
    num_files = length(file);
else
    num_files = 1;
    file = {file};
end
%% Import all spectra
found_idxs = cell(num_files,1);

TE_colors = crameri('devon', num_files+2); % +2 to avoid too-white trace
TM_colors = flipud(crameri('lajolla', num_files+2));
colors = TM_colors;
temps = zeros(1,num_files);
labels = [];

figure; hold on;
plot_shift = 1;
for file_idx = 1:num_files
    load(fullfile(location, file{file_idx}), 'lambdaArrayLooped', 'channel1Looped');
    % get temperature from filename
    split_name = split(file{file_idx},'.');
    temps(file_idx) = str2double(replace(split_name{1}, 'p', '.'));
    this_label = replace(split_name{1}, 'p', '.') + "C";
    labels = [labels; this_label];
    this_lambda = lambdaArrayLooped{1};
    this_t = channel1Looped{1};
    this_norm = this_t/max(this_t);
    [~, peak_idxs] = findpeaks(-this_norm, MinPeakProminence=0.5);
    plot(this_lambda,this_norm+plot_shift*file_idx, Color = colors(file_idx,:), DisplayName=this_label);
    plot(this_lambda(peak_idxs),this_norm(peak_idxs)+plot_shift*file_idx, '*', Color = colors(file_idx,:), HandleVisibility='off');
    found_idxs{file_idx} = peak_idxs;
    if(file_idx == 1)
        % assumption that all lambda are the same!
        imported_spectra = zeros(num_files, length(this_norm));
    end
    imported_spectra(file_idx,:) = this_norm+file_idx;
end
hold off; legend;
xlabel("Wavelength (nm)");
ylabel("Transmission, norm + shift");
%% Resonance tracking
% only consider resonances that appear in 1st spectrum
last_peaks = found_idxs{1};
approx_fsr_idx = mean(diff(last_peaks));
% each column of this should be the same resonance
sorted_peak_idxs = NaN(num_files, length(last_peaks));
sorted_peak_idxs(1,:) = last_peaks;
for file_idx = 2:num_files
    this_peaks = found_idxs{file_idx};
    for peak_idx = 1:length(this_peaks)
        this_peak = this_peaks(peak_idx);
        [nearest_dist, nearest_idx] = min(abs(last_peaks - this_peak));
        if(nearest_dist < approx_fsr_idx / 2) % criterion for this actually being the closest
            sorted_peak_idxs(file_idx,nearest_idx) = this_peak;
        end
    end
    last_peaks = sorted_peak_idxs(file_idx,:);
end
% trim incomplete resonances
sorted_peak_idxs(:,any(isnan(sorted_peak_idxs),1)) = [];
%% Plot tracked resonances (supp figure plot)
figure(Units="centimeters",Position=[6 6 9 4.5]);  hold on; 
hold on;
factor = 30; % downsample traces by this factor for lightweight vector plot
plot_lambda = downsample(this_lambda, factor);
for file_idx = num_files:-1:1 % flip order so legend is in same order
    this_norm = imported_spectra(file_idx,:);
    this_plot = downsample(this_norm, factor);
    plot(plot_lambda,this_plot, Color = colors(file_idx,:), DisplayName=labels(file_idx));
end
for peak_idx = 1:size(sorted_peak_idxs,2)
    this_idxs = sorted_peak_idxs(:,peak_idx);
    this_plot_idxs = sub2ind(size(imported_spectra), (1:num_files)', this_idxs);
    plot(this_lambda(this_idxs), imported_spectra(this_plot_idxs), 'k.-',HandleVisibility='off');
end
legend(Location = "westoutside");
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
xlabel("Wavelength (nm)"); ylabel("Transmission (norm. + shift)")
%% Calc effective index change
ring_length = 450e-6; % meters, used for group index calc
sorted_peak_lambda = 1e-9*this_lambda(sorted_peak_idxs);
[all_fsrs, ~] = gradient(sorted_peak_lambda);
avg_fsr = mean(all_fsrs,1);
avg_lambda = mean(sorted_peak_lambda,1);
c0 = 2.998e8;
ng = avg_lambda.^2./(avg_fsr*ring_length);
delta_lambdas = sorted_peak_lambda - sorted_peak_lambda(1,:);
delta_neffs = ng/avg_lambda * delta_lambdas;
avg_delta_neff = mean(delta_neffs,2);
figure; hold on;
plot(temps,avg_delta_neff, 'ko-');
%% save data
[output_filename, output_path] = uiputfile('*', 'Select location to save data:');
if(output_filename)
    save(strcat(output_path,output_filename), ...
        'temps', 'delta_neffs'); % freq domain scaling
else
    disp("File save cancelled");
end
