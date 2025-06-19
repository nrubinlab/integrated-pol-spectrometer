%% Plots of mode parameters
% Load data exported from FEMwell
clear;
simFile = "data/simulation/700x250_Si3N4_SU8_2025-05-28-20-28-20.mat";
load(simFile);
% Color palette
TEcolor = '#006eae'; TMcolor = '#df6464'; bothColor = '#12b31f';
%% -- CONVENIENCE PLOTS, THESE FIRST ONES ARE NOT USED IN PAPER -- %%
%% Plot loaded data
figure; 
surf(indexTemps, indexLambda, indexSimResult.TE', 'FaceColor', TEcolor);
hold on;
surf(indexTemps, indexLambda, indexSimResult.TM', 'FaceColor', TMcolor);
hold off;
legend("TE", "TM");
ylabel("Wavelength"); xlabel("Temperature"); zlabel("Effective Index");
zlim([0,3]);
% Do nth order 2-variable polynomial fit of data
polyOrder = 3;
%% Polynomial fit of index data
[meshTemps, meshLambda] = meshgrid(indexTemps, indexLambda);
vectorTemps = reshape(meshTemps,[],1);
vectorLambda = reshape(meshLambda,[],1);
sampleVector = [vectorTemps, vectorLambda];
teVector = reshape(indexSimResult.TE',[],1); % these transposes are really sneaky!
tmVector = reshape(indexSimResult.TM',[],1);
polyTE = polyfitFromFile(simFile, "TE");
polyTM = polyfitFromFile(simFile, "TM");
% residuals
figure;
plot3(vectorTemps, vectorLambda, (teVector-polyvaln(polyTE, sampleVector))./teVector, 'o', Color = TEcolor); hold on;
plot3(vectorTemps, vectorLambda, (tmVector-polyvaln(polyTM, sampleVector))./tmVector, 'o', Color = TMcolor); hold off;
xlabel("Temperature"); ylabel("Wavelength"); zlabel("Relative fit error");
%% Plot temp derivative of loaded data
derivPlotT = linspace(300,600,25); 
%derivPlotL = 1e-6*linspace(1.3,2,51);
derivPlotL = 1e-6*linspace(.6,1.1,51);
numT = length(derivPlotT); numL = length(derivPlotL);
[meshT, meshL] = meshgrid(derivPlotT, derivPlotL);
derivSampleVec = [reshape(meshT,[],1), reshape(meshL,[],1)];
teTempDeriv = polyvaln(polydern(polyTE, 1), derivSampleVec);
tmTempDeriv = polyvaln(polydern(polyTM, 1), derivSampleVec);
figure;
surf(meshT, meshL, reshape(teTempDeriv, numL, numT), 'FaceColor', TEcolor); hold on;
surf(meshT, meshL, reshape(tmTempDeriv, numL, numT), 'FaceColor', TMcolor); hold off;
legend("TE", "TM");
ylabel("Wavelength"); xlabel("Temperature"); zlabel("Effective Index TOC");
%zlim([0, 3e-4]);
%% Plot group index of loaded data
% group index
teLambdaDeriv = reshape(polyvaln(polydern(polyTE, 2), derivSampleVec), numL, numT);
tmLambdaDeriv = reshape(polyvaln(polydern(polyTM, 2), derivSampleVec), numL, numT);
teNeff = reshape(polyvaln(polyTE, derivSampleVec), numL, numT);
tmNeff = reshape(polyvaln(polyTM, derivSampleVec), numL, numT);
ngroupTE = teNeff - meshL.*teLambdaDeriv;
ngroupTM = tmNeff - meshL.*tmLambdaDeriv;
figure;
surf(meshT, meshL, ngroupTE, 'FaceColor', TEcolor); hold on;
surf(meshT, meshL, ngroupTM, 'FaceColor', TMcolor); hold off;
legend("TE", "TM");
ylabel("Wavelength"); xlabel("Temperature"); zlabel("Group Index");
%% Group index TOC
% doing a rough sample based derivative for now
dngdT_TE = diff(ngroupTE, 1, 2)./diff(meshT, 1, 2);
dngdT_TM = diff(ngroupTM, 1, 2)./diff(meshT, 1, 2);
figure;
surf(derivPlotT(1:end-1), derivPlotL, dngdT_TE, 'FaceColor', TEcolor); hold on;
surf(derivPlotT(1:end-1), derivPlotL, dngdT_TM, 'FaceColor', TMcolor); hold off;
legend("TE", "TM");
ylabel("Wavelength"); xlabel("Temperature"); zlabel("Group Index TOC");
%% neff vs lambda
figure("Units", "centimeters", Position=[7.5 7.5 3 2]);
hold on;
tempIdx = 1;
plot(1e9*indexLambda, indexSimResult.TE(tempIdx,:), 'Color', TEcolor);
plot(1e9*indexLambda, indexSimResult.TM(tempIdx,:), 'Color', TMcolor);
hold off;
ylim([0,3]);
%xlim([1300, 1900]);
xticks(1300:300:1900);
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
%% -- PLOTS USED IN PAPER -- %%
% neff TOC vs lambda
plotLambda = linspace(min(indexLambda), 1.1e-6);
derivSampleVec = [300*ones(size(plotLambda')), plotLambda'];
teNeffTOC = polyvaln(polydern(polyTE, 1), derivSampleVec);
tmNeffTOC = polyvaln(polydern(polyTM, 1), derivSampleVec);
figure("Units", "centimeters", Position=[7.5 7.5 6 6]);
hold on;
plot(1e9*plotLambda, teNeffTOC, 'Color', TEcolor, DisplayName = 'TE');
plot(1e9*plotLambda, tmNeffTOC, 'Color', TMcolor, DisplayName = 'TM');
yline(0, 'k:', HandleVisibility = 'off');
hold off;
%ylim([0,2.5e-4]);
%xlim([1300, 1900]);
%xticks(1300:300:1900);
set(gca, 'fontsize', 10, 'ticklength', [0.03, 0.03]);
xlabel("Wavelength (nm)");
ylabel("Effective index TOC");
xlim([780, 1100]);
legend(Location = 'southwest');
%% ng TOC vs lambda
% THIS SECTION WILL NOT WORK UNLESS YOU RUN *ALL* THE 'CONVENIENCE
% PLOT' SECTIONS!
dngdT_TE = diff(ngroupTE, 1, 2)./diff(meshT, 1, 2);
dngdT_TM = diff(ngroupTM, 1, 2)./diff(meshT, 1, 2);

% show alongside ng the resolution assuming these spectrometer params
L = 10 * 1e-2; % 
Tmax = 100;
% the resolution expression in the paper is frequency but this is the
% wavenumber version
TE_res = 1 ./ (2*L*Tmax*dngdT_TE(:,1)');
TM_res = 1 ./ (2*L*Tmax*dngdT_TM(:,1)');

figure("Units", "centimeters", Position=[7.5 7.5 6 6]);
hold on;
tempIdx = 1;
plot(1e9*derivPlotL, dngdT_TE(:,1), 'Color', TEcolor, DisplayName = 'TE');
plot(1e9*derivPlotL, dngdT_TM(:,1), 'Color', TMcolor, DisplayName = 'TM');
ylabel("Group index TOC");
hold off;
%ylim([0,1e-3]);
xlim([780, 1050]);
%xticks(1300:300:1900);
set(gca, 'fontsize', 10, 'ticklength', [0.03, 0.03]);
xlabel("Wavelength (nm)");
legend(Location = 'southwest');
%% Resolution calculation
figure("Units", "centimeters", Position=[7.5 7.5 6 6]);
hold on;
plot(1e9*derivPlotL, 1e-2*TE_res, 'Color', TEcolor, DisplayName = 'TE');
plot(1e9*derivPlotL, 1e-2*TM_res, 'Color', TMcolor, DisplayName = 'TM');
xlabel("Wavelength (nm)");
ylabel("Resolution (cm^{-1})");
xlim([780, 1050]);
legend(Location = 'northwest');
set(gca, 'fontsize', 10, 'ticklength', [0.03, 0.03]);
%% These plots are not used but are helpful to visualize the separable BW
polyTE = polyfitFromFile(simFile, "TE");
polyTM = polyfitFromFile(simFile, "TM");
lambdaMin = 0.6e-6; lambdaMax = 0.8e-6;
lowLambdas = linspace(lambdaMin, lambdaMax);
lambda0 = 1.555e-6; t0 = 350;
thisQ_TE = calcQ(polyTE, t0, lambda0, lowLambdas);
thisH = calcH(polyTE, polyTM, t0, lowLambdas);

upperLambdas = zeros(size(lowLambdas));
for i = 1:length(lowLambdas)
    [~, idx] = min(abs(thisH(i)*thisQ_TE(i) - thisQ_TE));
    upperLambdas(i) = lowLambdas(idx);
end
centerLambdas = (lowLambdas + upperLambdas)/2;
bandwidth = upperLambdas - lowLambdas;
thisBandwidth = interp1(centerLambdas, bandwidth, lambda0);

figure("Units", "centimeters", Position=[7.5 7.5 3 2]); hold on;
plot(1e9*centerLambdas, 1e9*bandwidth, 'k', LineWidth = 1);
hold off;
ylim([0, 400]);
%xlim([1400, 1700]);
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
%%
figure;
plot(lowLambdas, lowLambdas);   
plot(lowLambdas, upperLambdas);
function q = calcQ(polyFit, t0, lambda0, lambda)
    % nu can be 1d array but nu0 and t0 cannot be
    % nu must be column vector
    lambda = reshape(lambda, [], 1);
    numerator = polyvaln(polydern(polyFit, 1), [t0*ones(size(lambda)), lambda]);
    denominator = polyvaln(polydern(polyFit, 1), [t0, lambda0]);
    g = numerator/denominator;
    q = (2.998e8./lambda).*g;
end

function h = calcH(polyFitTE, polyFitTM, t0, lambda)
    % nu can be 1d array but nu0 and t0 cannot be
    % nu must be column vector
    lambda = reshape(lambda, [], 1);
    numerator = polyvaln(polydern(polyFitTM, 1), [t0*ones(size(lambda)), lambda]);
    denominator = polyvaln(polydern(polyFitTE, 1), [t0*ones(size(lambda)), lambda]);
    h = numerator./denominator;
end