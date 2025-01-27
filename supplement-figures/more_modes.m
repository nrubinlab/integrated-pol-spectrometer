%% Load data exported from FEMwell
clear;
sim_file = "data/simulation/800x220_Si_4modes_2024-12-03-10-35-15.mat";
load(sim_file);
% --- Use this line for fig S11 --- %
% mode_names = ["TE1", "TE2", "TM1"]; plot_lambda = linspace(1.4e-6, 1.7e-6,300);
% --------------------------------- %
% --- Use this line for fig S12 --- %
mode_names = ["TE1", "TE2", "TM1", "TM2"]; plot_lambda = linspace(1.4e-6, 1.65e-6,300);
% --------------------------------- %
colors = ["#006eae", "#00ae5b", "#df6464", "#c952c2"];
%% Plot loaded data
figure; 
for mode_idx = 1:length(mode_names)
    surf(indexTemps, indexLambda, indexSimResult.(mode_names(mode_idx))', ...
        'FaceColor', colors(mode_idx));
    if(mode_idx == 1)
        hold on;
    end
end
ylabel("Wavelength"); xlabel("Temperature"); zlabel("Effective Index");
zlim([0,3]);
%% neff plot
figure("Units", "centimeters", Position=[7.5 7.5 9 4.5]); hold on;
tempIdx = 1;
for mode_idx = 1:length(mode_names)
    this_data = indexSimResult.(mode_names(mode_idx));
    plot(1e9*indexLambda, this_data(tempIdx,:), ...
        'Color', colors(mode_idx), ...
        DisplayName=mode_names(mode_idx));
end
yline(1.44, 'k--', HandleVisibility='off'); % cladding index
% hold off;
xlim('tight');
ylim([1,3]);
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
xlabel("Wavelength (nm)");
ylabel("Effective Index");
legend(Location = 'westoutside');
%% dneff/dT plot

figure("Units", "centimeters", Position=[7.5 7.5 6 4.5]); hold on;
for mode_idx = 1:length(mode_names)
    this_mode = mode_names(mode_idx);
    this_neff_toc = neff_toc_from_file(sim_file, this_mode, plot_lambda);
    plot(1e9*plot_lambda, this_neff_toc, ...
        'Color', colors(mode_idx), ...
        DisplayName=this_mode);
end
hold off;
ylim([0,2.5e-4]);
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
xlabel("Wavelength (nm)");
ylabel("dn_{eff}/dT (K^{-1})")
%% dng/dT plot
figure("Units", "centimeters", Position=[7.5 7.5 6 4.5]); hold on;
for mode_idx = 1:length(mode_names)
    this_mode = mode_names(mode_idx);
    this_ng_toc = ng_toc_from_file(sim_file, this_mode, plot_lambda);
    plot(1e9*plot_lambda, this_ng_toc, ...
        'Color', colors(mode_idx), ...
        DisplayName=this_mode);
end
hold off;
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
xlabel("Wavelength (nm)");
ylabel("dn_{g}/dT (K^{-1})")
%% resolution, normalized for temperature difference and spectrometer length
figure("Units", "centimeters", Position=[7.5 7.5 6 4.5]); hold on;
for mode_idx = 1:length(mode_names)
    this_mode = mode_names(mode_idx);
    this_ng_toc = ng_toc_from_file(sim_file, this_mode, plot_lambda);
    this_res = res_from_ng_toc(this_ng_toc, plot_lambda);
    % 1e2 converts m to cm, 1e9 converts m to nm
    plot(1e9*plot_lambda, 1e2*1e9*this_res, ...
        'Color', colors(mode_idx), ...
        DisplayName=this_mode);
end
hold off;
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
xlabel("Wavelength (nm)");
ylabel("Resolution FOM (nm·K·cm)");
%%
% get all unique combinations of two modes
%combinations = nchoosek(["TE1", "TE2", "TM1", "TM2"],2);
combinations = nchoosek(mode_names,2);
figure("Units", "centimeters", Position=[7.5 7.5 9 4.5]); hold on;
min_upper_lambda = [];
for i = 1:length(combinations)
    upper_lambda = BW_for_mode_pair(sim_file, combinations(i,1), combinations(i,2), plot_lambda);
    plot(1e9*plot_lambda, 1e9*upper_lambda, DisplayName=sprintf("%s-%s", combinations(i,1), combinations(i,2)));
    if(i == 1)
        min_upper_lambda = upper_lambda;
    else
        min_upper_lambda = min(min_upper_lambda, upper_lambda);
    end
end
plot(1e9*plot_lambda, 1e9*min_upper_lambda, 'k--', LineWidth = 1.5, DisplayName="All");
plot(1e9*plot_lambda,1e9*plot_lambda,'k:', HandleVisibility='off');
hold off; legend(Location = "westoutside");
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
xlabel("Lower Wavelength (nm)");
ylabel("Upper wavelength (nm)");
%% plot 'all 4' separable BW using center vs. BW type plot
figure("Units", "centimeters", Position=[7.5 7.5 5 4.5]); hold on;
plot(1e9*(plot_lambda+min_upper_lambda)/2, 1e9*(min_upper_lambda - plot_lambda), 'k');
xlabel("Center wavelength (nm)"); ylabel("Bandwidth (nm)");
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
xlim('tight');
function upper_lambda = BW_for_mode_pair(sim_file, mode1, mode2, lower_lambda)
    % for the set of provided lower lambdas, returns the max lambda of the
    % seprable BW (or the max lambda of the simulation, whichever is lower)
    c0 = physconst("LightSpeed");

    neff_TOC_1 = neff_toc_from_file(sim_file, mode1, lower_lambda);
    neff_TOC_2 = neff_toc_from_file(sim_file, mode2, lower_lambda);

    f_1 = (c0./lower_lambda').*neff_TOC_1;
    f_2 = (c0./lower_lambda').*neff_TOC_2;

    upper_lambda = zeros(size(lower_lambda));
    for i = 1:length(lower_lambda)
        if(f_1(i) > f_2(i)) % check which is higher
            % solve inequality by finding crossing
            [~, idx] = min(abs(f_2(i) - f_1));
        else
            [~, idx] = min(abs(f_1(i) - f_2));
        end
        upper_lambda(i) = lower_lambda(idx);
    end
end

function out = res_from_ng_toc(ng_toc, lambda)
    % this gives the wavelength resolution without the L or delta T terms
    % to get the actual resolution, you take this number and divide by the
    % spectrometer length times max temp difference
    c0 = physconst("lightspeed");
    a = 2; % FWHM constant for Hann window
    dv = a*c0./(2*ng_toc);
    out = dv'.*lambda.^2/c0;
end

function out = ng_toc_from_file(sim_file, mode_string, lambda)
    thisPoly = polyfitFromFile(sim_file, mode_string);
    derivPlotT = linspace(300,600,25); % hardcoding out of lazyness
    numT = length(derivPlotT); numL = length(lambda);
    [meshT, meshL] = meshgrid(derivPlotT, lambda);
    derivSampleVec = [reshape(meshT,[],1), reshape(meshL,[],1)];
    lambdaDeriv = reshape(polyvaln(polydern(thisPoly, 2), derivSampleVec), numL, numT);
    thisNeff = reshape(polyvaln(thisPoly, derivSampleVec), numL, numT);
    ngroup = thisNeff - meshL.*lambdaDeriv;
    dngdT_TE = diff(ngroup, 1, 2)./diff(meshT, 1, 2);
    out = dngdT_TE(:,1);
end

function out = neff_toc_from_file(sim_file, mode_string, lambda)
    derivSampleVec = [300*ones(size(lambda')), lambda'];
    thisPoly = polyfitFromFile(sim_file, mode_string);
    out = polyvaln(polydern(thisPoly, 1), derivSampleVec);
end