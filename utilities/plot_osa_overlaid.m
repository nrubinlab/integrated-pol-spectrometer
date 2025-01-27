[file, location] = uigetfile('.mat', 'Select One or More Files', 'MultiSelect', 'on');
if(iscell(file))
    numFiles = length(file);
else
    numFiles = 1;
    file = {file};
end
%%
colors = cool(numFiles);
plot_cropped = false;
plot_lin = true;
figure;
for i = 1:numFiles
    load(fullfile(location, file{i}), ...
        'osa_lambda', 'osa_power_dbm');
    osa_power_lin = 10.^(osa_power_dbm/10);
    if(plot_lin), plot_data = osa_power_lin; else, plot_data = osa_power_dbm; end
    plot(osa_lambda,plot_data, 'DisplayName', file{i}, ...
        Color = colors(i,:));
    if(i == 1)
        hold on;
    end
    %pause
end
hold off; legend(Interpreter = "none");
xlabel("Wavelength (nm)");
%ylim([-50,-5]);