[file, location] = uigetfile('.mat', 'Select One or More Files', 'MultiSelect', 'on');
if(iscell(file))
    numFiles = length(file);
else
    numFiles = 1;
    file = {file};
end
%%
start_time = 0.91; end_offset = 0.11; % this is chosen by user as needed
%start_time = .91; end_time = 11; % this is chosen by user as needed
%time_crop_idxs = [2550:7500];
plot_cropped = true;
figure;
for i = 1:numFiles
    load(fullfile(location, file{i}), ...
        'acq_struct', 'channel1', 'channel2');
    [~,time_array,~] = generate_detector_params(acq_struct);
    end_time = time_array(end) - end_offset;
    time_crop_idxs = time_array > start_time & time_array < end_time;
    %thisDetrended = channel1./movmean(channel1,50) - 1;
    thisData = channel1;
    %thisDetrended = channel1 - movmean(channel1,100);
    if(plot_cropped)
        plot(time_array(time_crop_idxs), thisData(time_crop_idxs), 'DisplayName', file{i});
    else
        plot(time_array, channel1, 'DisplayName', file{i});
    end
    
    if(i == 1)
        hold on;
    end
end
hold off; legend();
xlabel("Time (s)"); ylabel("Optical power (W)");
%% Naive FFTs
zeropadfactor = 8;
doNormalize = true;
maxPower = 0;
lambda_name = false;
figure(Units = "inches", Position=[3 3 5 2.5]); 
colors = crameri('devon', numFiles+2); % +2 to avoid too-white trace
for i = 1:numFiles
   load(fullfile(location, file{i}), ...
        'acq_struct', 'lambda', 'channel1', 'channel2');
    [~,time_array,~] = generate_detector_params(acq_struct);
    end_time = time_array(end) - end_offset;
    time_crop_idxs = time_array > start_time & time_array < end_time;
    thisMax = max(thisData);
    if(thisMax > maxPower)
        maxPower = thisMax;
    end
    %thisData = channel1./movmean(channel1,75) - 1;
    %thisDetrended = channel1 - movmean(channel1,100);
    thisData = channel1(time_crop_idxs); % don't normalize when we're worried about amplitude
    thisData = highpass(thisData,0.05);
    fftLength = length(thisData)*zeropadfactor;
    thisFFT = abs(fft(thisData.*hamming(length(thisData))', fftLength));
    thisIdx = 1:round(length(thisFFT)/2);
    N = length(thisFFT);
    this_f = (0:(N-1))/(N*acq_struct.sampling_interval);
    %this_f = (0:(N-1))/(N);
    if(exist('lambda', 'var') && lambda_name)
        this_name = sprintf("%1.0fnm", 1e9*lambda);
    else
        this_name = file{i};
    end
    this_color = colors(i,:);
    
    if(doNormalize)
        this_data = thisFFT(thisIdx)/max(thisFFT(thisIdx));
    else
        this_data = thisFFT(thisIdx);
    end
    semilogy(this_f(thisIdx),this_data, 'DisplayName', this_name,...
        'Color',this_color);
    if(i == 1)
        hold on;
    end
end
hold off; legend("Interpreter","none");
%xlabel("Normalized frequency (a.u.)");
xlabel("Frequency (Hz)");
ylabel("Power spectral density");
t = title(sprintf("First file %s, max power = %1.2f uW", file{1}, 1e6*maxPower), "Interpreter","none"); 