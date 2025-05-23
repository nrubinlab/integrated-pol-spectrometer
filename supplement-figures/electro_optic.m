clear;
location = 'data\20240813\';
%% Time domain - use data with narrow time range
file = {'narrow_pol1.CSV', 'narrow_pol2.CSV', 'narrow_saw.CSV'};
numFiles = length(file);
yIncrement = -1;
figure("Units", "centimeters", Position=[7.5 7.5 9 4.5]); hold on;
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
reduc_factor = 10;
for i = 1:numFiles
    trace = importTektronix(fullfile(location, file{i}));
    thisOffset = (i-1)*yIncrement;

    thisPlot = downsample(normalize(trace.V)/5, reduc_factor);
    thisTime = downsample(trace.t, reduc_factor);
    %thisPlot = trace.V;
    plot(1e6*thisTime, thisPlot + thisOffset, 'DisplayName', file{i});
    yline(thisOffset, 'k:', 'HandleVisibility','off');
    clear trace
end
legend('interpreter', 'none', Location = "Southeast");
xlabel("Time (us)");
ylabel("Measured signal (a.u. + offset)");
legend("Polarization 1", "Polarization 2", "Sawtooth Input");
%% FFTs - use data with wide time range
file = {'wide_pol1.CSV', 'wide_pol2.CSV'};
numFiles = length(file);
figure("Units", "centimeters", Position=[7.5 7.5 15 4.5]);
tiledlayout(1,2, "TileSpacing","tight","Padding","tight");
nexttile(1); 
nexttile(2); 
horiz_offset = 1e-2; % MHz
for i = 1:numFiles
    trace = importTektronix(fullfile(location, file{i}));
    %trace_zeroed = trace.V - mean(trace.V);
    trace_zeroed = normalize(trace.V);
    thisFFT = abs(fft(trace_zeroed)).^2;
    thisDt = trace.t(2)-trace.t(1); % assumes uniform
    N = length(thisFFT);
    thisFreqAxis = (0:N-1)/(N*thisDt); thisMaxN = round(N/2);
    %thisPlot = thisFFT(1:thisMaxN).*thisFreqAxis(1:thisMaxN);
    nexttile(1);
    semilogy(1e-6*thisFreqAxis(1:thisMaxN) + i*horiz_offset, thisFFT(1:thisMaxN), 'DisplayName', file{i});
    nexttile(2);
    semilogy(1e-6*thisFreqAxis(1:thisMaxN) + i*horiz_offset, thisFFT(1:thisMaxN), 'DisplayName', file{i});
    if(i == 1)
        nexttile(1); hold on; xlim('tight'); xlabel("Frequency (MHz)");
        set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
        ylim([1e2, 1e8]); ylabel("|FFT|^2 (a.u.)")
        nexttile(2); hold on; xlim([0,1.5]); xlabel("Frequency (MHz)");
        set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
        ylim([1e2, 1e8]); 
    end
    clear trace
end
nexttile(1); legend("Polarization 1", "Polarization 2", Location = "westoutside");

function trace = importTektronix(filename, dataLines)
    %IMPORTFILE Import data from a text file
    %  TRACE = IMPORTFILE(FILENAME) reads data from text file FILENAME for
    %  the default selection.  Returns the data as a table.
    %
    %  TRACE = IMPORTFILE(FILE, DATALINES) reads data for the specified row
    %  interval(s) of text file FILENAME. Specify DATALINES as a positive
    %  scalar integer or a N-by-2 array of positive scalar integers for
    %  dis-contiguous row intervals.
    %
    %  Example:
    %  trace = importfile("C:\Users\Kyle\Documents\Polarization_Spectrometer\electro-optic\20240813\TEK00003.CSV", [1, Inf]);
    %
    %  See also READTABLE.
    %
    % Auto-generated by MATLAB on 13-Aug-2024 20:32:36
    
    %% Input handling
    
    % If dataLines is not specified, define defaults
    if nargin < 2
        dataLines = [1, Inf];
    end
    
    %% Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 2);
    
    % Specify range and delimiter
    opts.DataLines = dataLines;
    opts.Delimiter = ",";
    
    % Specify column names and types
    opts.VariableNames = ["t", "V"];
    opts.VariableTypes = ["double", "double"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % Import the data
    trace = readtable(filename, opts);

end