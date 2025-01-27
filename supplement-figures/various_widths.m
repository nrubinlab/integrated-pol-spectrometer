%% Load data
clear;
data_dir = 'data\20240922\organized_for_supplement';
dir_list = dir(data_dir);
% remove '.' and '..'
dir_list = dir_list(3:end);
%% Iterate through and plot h for each waveguide geometry
c0 = physconst("LightSpeed");
lambdaExp = 1e-9*(1480:30:1630);
nuExp = c0./lambdaExp;

lambdaSim = linspace( 1.3e-6, 1.9e-6);
lambda0 = 1.55e-6;
t0 = 350; % this is just an approximation
% to find factor to adjust 0.5/1.0 spectrometer,
% using area ratio needed to predict area of spiral from 1.0 and 0.5
% spirals
% mixed: 25861 um
% 0.5:14741 -> 2.948 cm length
% 1.0 would be 2x the area, or 29482
%
area_mixed = 25861; % um^2
area_0p5 = 14741; area_1p0 = 2*area_0p5;
tapered_wg_ratio = (area_mixed + area_0p5)/(area_1p0+area_0p5);
sim_0p5 = 'data\simulation\500x220_Si_2024-11-18-20-15-00.mat';
hSim_0p5 = calcH_from_sim_file(sim_0p5, t0, lambdaSim);

figure(Units="centimeters",Position=[6 6 5.5 4.5]);  hold on;
for dir_idx = 1:length(dir_list)
    this_dir = fullfile(dir_list(dir_idx).folder, dir_list(dir_idx).name);
    % extract waveguide width from naming convention
    
    split_str = split(dir_list(dir_idx).name, {'-','_'});
    width_string = split_str{2};
    this_width = str2double(replace(width_string,'p','.'));

    TE_calib_file = fullfile(this_dir, 'TE.mat');
    TM_calib_file = fullfile(this_dir, 'TM.mat');
    sim_file = fullfile(this_dir, 'simulation.mat');

    hExp = calcH_from_calib(TE_calib_file, TM_calib_file, nuExp);
    
    % H from simulation
    hSim = calcH_from_sim_file(sim_file, t0, lambdaSim);
    if(this_width == 1)
        % apply correction on wide waveguides, which taper in and out of
        % this width (standard width is 0.5)
        hSim = tapered_wg_ratio*hSim + (1-tapered_wg_ratio)*hSim_0p5;
    end
    %
    p1 = plot(1e9*lambdaExp, hExp, '.', DisplayName=sprintf('w= %1.1f um exp.', this_width)); 
    plot(1e9*lambdaSim, hSim, Color = p1.Color, DisplayName=sprintf('w = %1.1f sim.', this_width));
end
hold off;
xlim([1400,1700]);
xlabel("Wavelength (nm)");
ylabel("h Parameter")
TEcolor = '#006eae'; TMcolor = '#df6464'; bothColor = '#12b31f';
%legend('location', 'southwest');
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
%% g(nu)
lambda0 = 1.54e-6;
t0 = 350; % this is just an approximation

sim_0p5 = 'data\simulation\500x220_Si_2024-11-18-20-15-00.mat';
[g_TE_sim_0p5, g_TM_sim_0p5] = calcG_from_sim(sim_0p5, t0, lambda0, lambdaSim);

figure(Units="centimeters",Position=[6 6 14 4.5]); 
tiledlayout(1,2,"TileSpacing","tight","Padding","tight");
nexttile(1); hold on;
xlim([1400,1700]); xlabel("Wavelength (nm)");
ylabel("g_{TE}"); legend('location', 'westoutside');
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);

nexttile(2); hold on;
xlim([1400,1700]); xlabel("Wavelength (nm)");
ylabel("g_{TM}"); %legend('location', 'southwest');
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
for dir_idx = 1:length(dir_list)
    this_dir = fullfile(dir_list(dir_idx).folder, dir_list(dir_idx).name);
    % extract waveguide width from naming convention
    
    split_str = split(dir_list(dir_idx).name, {'-','_'});
    width_string = split_str{2};
    this_width = str2double(replace(width_string,'p','.'));

    TE_calib_file = fullfile(this_dir, 'TE.mat');
    TM_calib_file = fullfile(this_dir, 'TM.mat');
    sim_file = fullfile(this_dir, 'simulation.mat');

    [g_TE_exp, g_TM_exp] = calcG_from_calib(TE_calib_file, TM_calib_file, nuExp);
    
    % H from simulation
    [g_TE_sim, g_TM_sim] = calcG_from_sim(sim_file, t0, lambda0, lambdaSim);
    if(this_width == 1.0)
        % apply correction on wide waveguides, which taper in and out of
        % this width (standard width is 0.5)
        g_TE_sim = tapered_wg_ratio*g_TE_sim + (1-tapered_wg_ratio)*g_TE_sim_0p5;
        g_TM_sim = tapered_wg_ratio*g_TM_sim + (1-tapered_wg_ratio)*g_TM_sim_0p5;
    end
    %
    nexttile(1);
    p1 = plot(1e9*lambdaExp, g_TE_exp, '.', DisplayName=sprintf('w= %1.2f um exp.', this_width)); 
    plot(1e9*lambdaSim, g_TE_sim, Color = p1.Color, DisplayName=sprintf('w = %1.2f um sim.', this_width));
    nexttile(2);
    p1 = plot(1e9*lambdaExp, g_TM_exp, '.', DisplayName=sprintf('w= %1.2f um exp.', this_width)); 
    plot(1e9*lambdaSim, g_TM_sim, Color = p1.Color, DisplayName=sprintf('w = %1.2f um sim.', this_width));
end
%% Resolution, experimental
extrap_lambda = linspace(1.48e-6, 1.63e-6);

figure(Units="centimeters",Position=[6 6 12 4.5]); 
tiledlayout(1,2,"TileSpacing","tight","Padding","tight");
nexttile(1); hold on;
xlim([1475,1635]); xlabel("Wavelength (nm)");
ylabel("TE resolution (nm)"); %legend('location', 'westoutside');
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);

nexttile(2); hold on;
xlim([1475,1635]); xlabel("Wavelength (nm)");
ylabel("TM resolution (nm)"); %legend('location', 'southwest');
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);
for dir_idx = 1:length(dir_list)
    this_dir = fullfile(dir_list(dir_idx).folder, dir_list(dir_idx).name);
    % extract waveguide width from naming convention
    
    split_str = split(dir_list(dir_idx).name, {'-','_'});
    width_string = split_str{2};
    this_width = str2double(replace(width_string,'p','.'));

    TE_calib_file = fullfile(this_dir, 'TE.mat');
    TM_calib_file = fullfile(this_dir, 'TM.mat');
    TE_data_dir = fullfile(this_dir, 'sweeps', 'TE');
    TM_data_dir = fullfile(this_dir, 'sweeps', 'TM');
    sim_file = fullfile(this_dir, 'simulation.mat');

    TE_theory = calculate_theoretical_resolution(TE_calib_file, c0./extrap_lambda);
    TM_theory = calculate_theoretical_resolution(TM_calib_file, c0./extrap_lambda);

    %[TE_sim, TM_sim] = calculate_sim_resolution(sim_file, TE_calib_file, TM_calib_file, t0, lambda0, lambdaSim);

    [TE_lambdas, TE_widths] = measure_reconstructed_fwhms(TE_data_dir, TE_calib_file);
    [TM_lambdas, TM_widths] = measure_reconstructed_fwhms(TM_data_dir, TM_calib_file);

    %
    nexttile(1);
    %p1 = plot(1e9*lambdaExp, g_TE_exp, '.', DisplayName=sprintf('w= %1.1f um exp.', this_width)); 
    %plot(1e9*lambdaSim, g_TE_sim, Color = p1.Color, DisplayName=sprintf('w = %1.1f sim.', this_width));
    p1 = plot(TE_lambdas, TE_widths, '.', DisplayName=sprintf('w= %1.1f um exp.', this_width));
    % dummy plot for colors
    plot(0,0);
    %plot(1e9*extrap_lambda, 1e9*TE_theory, DisplayName=sprintf('w= %1.1f um theor.', this_width), Color = p1.Color);
    %plot(1e9*lambdaSim, 1e9*TE_sim, DisplayName=sprintf('w= %1.1f um sim.', this_width), Color = p1.Color);
    
    nexttile(2);
    %p1 = plot(1e9*lambdaExp, g_TM_exp, '.', DisplayName=sprintf('w= %1.1f um exp.', this_width)); 
    %plot(1e9*lambdaSim, g_TM_sim, Color = p1.Color, DisplayName=sprintf('w = %1.1f sim.', this_width));
    p2 = plot(TM_lambdas, TM_widths, '.', DisplayName=sprintf('w= %1.1f um exp.', this_width));
    plot(0,0);
    %plot(1e9*extrap_lambda, 1e9*TM_theory, DisplayName=sprintf('w= %1.1f um theor.', this_width), Color = p1.Color);
    %plot(1e9*lambdaSim, 1e9*TM_sim, DisplayName=sprintf('w= %1.1f um sim.', this_width), Color = p2.Color);
end

%% Separable bandwidth
figure(Units="centimeters",Position=[6 6 6.8 4.5]);  hold on;
BW_0p5 = calcBW_from_sim(sim_0p5, t0, lambda0, lambdaSim);
for dir_idx = 1:length(dir_list)
    this_dir = fullfile(dir_list(dir_idx).folder, dir_list(dir_idx).name);
    % extract waveguide width from naming convention
    split_str = split(dir_list(dir_idx).name, {'-','_'});
    width_string = split_str{2};
    this_width = str2double(replace(width_string,'p','.'));

    sim_file = fullfile(this_dir, 'simulation.mat');
    [this_center, this_BW] = calcBW_from_sim(sim_file, t0, lambda0, lambdaSim);

%     if(this_width == 1)
%         % apply correction on wide waveguides, which taper in and out of
%         % this width (standard width is 0.5) - this is just a rough
%         % approximation
%         this_BW = tapered_wg_ratio*this_BW + (1-tapered_wg_ratio)*BW_0p5;
%     end
    %
    p1 = plot(1e9*this_center, 1e9*this_BW, DisplayName=sprintf('w= %1.1f um exp.', this_width)); 
    % dummy plot to get same colors as all the others :)
    plot(0,0);
    
end
plot(1555, 150, '.'); % demonstrated BW in this work
hold off;
xlim([1400,1700]);
xlabel("Center Wavelength (nm)");
ylabel("Separable BW (nm)")
TEcolor = '#006eae'; TMcolor = '#df6464'; bothColor = '#12b31f';
%legend('location', 'southwest');
set(gca, 'fontsize', 7, 'ticklength', [0.03, 0.03]);

function h = calcH_from_sim_file(sim_file, t0, lambda)

    polyTE = polyfitFromFile(sim_file, "TE");
    polyTM = polyfitFromFile(sim_file, "TM");

    % nu can be 1d array but nu0 and t0 cannot be
    % nu must be column vector
    lambda = reshape(lambda, [], 1);
    % polyTM and polyTM are just effective index vs temp and frequency, but
    % in the simplest case of a balanced spectrometer, h reduces to a ratio
    % of the TOC's
    numerator = polyvaln(polydern(polyTM, 1), [t0*ones(size(lambda)), lambda]);
    denominator = polyvaln(polydern(polyTE, 1), [t0*ones(size(lambda)), lambda]);

    h = numerator./denominator;
end

function [g_TE, g_TM] = calcG_from_sim(sim_file, t0, lambda0, lambda)
    polyTE = polyfitFromFile(sim_file, "TE");
    polyTM = polyfitFromFile(sim_file, "TM");

    % lambda can be 1d array but nu0 and t0 cannot be
    % nu must be column vector
    lambda = reshape(lambda, [], 1);
    % polyTM and polyTM are just effective index vs temp and frequency, but
    % in the simplest case of a balanced spectrometer, h reduces to a ratio
    % of the TOC's
    g_TE = polyvaln(polydern(polyTE, 1), [t0*ones(size(lambda)), lambda])...
        /polyvaln(polydern(polyTE, 1), [t0, lambda0]);
    g_TM = polyvaln(polydern(polyTM, 1), [t0*ones(size(lambda)), lambda])...
        /polyvaln(polydern(polyTM, 1), [t0, lambda0]);
end

function [res_TE, res_TM] = calculate_sim_resolution(sim_file, calib_TE, calib_TM, t0, lambda0, lambda)
    % calculate simulated resolution using dispersion information from
    % simulation, but max OPL data from experiment - there is a bit too
    % much uncertainty on the actual temperature difference achieved to
    % determine OPL only from simulation, and that is not the purpose of
    % these plots
    c0 = physconst("LightSpeed");
    load(calib_TE, 'opl_poly', 'start_time', 'end_time');
    u_range_te = abs(polyval(opl_poly, start_time) - polyval(opl_poly, end_time));
    load(calib_TM, 'opl_poly', 'start_time', 'end_time');
    u_range_tm = abs(polyval(opl_poly, start_time) - polyval(opl_poly, end_time));

    [g_TE, g_TM] = calcG_from_sim(sim_file, t0, lambda0, lambda);
    % lazy derivative of g, I really don't want to deal with more polyders
    % note we have to do derivative w/ respect to frequencys
    nu = c0./lambda';
    g_TE_deriv = gradient(g_TE)./gradient(nu);
    g_TM_deriv = gradient(g_TM)./gradient(nu);

    a = 2; % a value for the Hann window
    resolution_freq_TE = a*c0/u_range_te * (g_TE + nu.*g_TE_deriv).^-1;
    res_TE = c0*resolution_freq_TE./(nu.^2);

    resolution_freq_TM = a*c0/u_range_tm * (g_TM + nu.*g_TM_deriv).^-1;
    res_TM = c0*resolution_freq_TM./(nu.^2);
end

function [centerLambdas, bandwidth] = calcBW_from_sim(sim_file, t0, lambda0, lambda)
    c0 = physconst("LightSpeed");

    [g_TE, ~] = calcG_from_sim(sim_file, t0, lambda0, lambda);
    q_TE = (c0./lambda').*g_TE;
    thisH = calcH_from_sim_file(sim_file, t0, lambda);
    upperLambdas = zeros(size(lambda));
    for i = 1:length(lambda)
        [~, idx] = min(abs(thisH(i)*q_TE(i) - q_TE));
        upperLambdas(i) = lambda(idx);
    end
    centerLambdas = (lambda + upperLambdas)/2;
    bandwidth = upperLambdas - lambda;
end

function [g_TE, g_TM] = calcG_from_calib(calibTE, calibTM, nu_eval)
    TE = load(calibTE, 'opl_poly', 'start_time', 'end_time',...
        'g_poly', 'g_mu', 'nu_min', 'nu_max');
    TM = load(calibTM, 'opl_poly', 'start_time', 'end_time',...
        'g_poly', 'g_mu', 'nu_min', 'nu_max');
    %
    g_TE = polyval(TE.g_poly, nu_eval, [], TE.g_mu);
    g_TM = polyval(TM.g_poly, nu_eval, [], TM.g_mu);
end

function h = calcH_from_calib(calibTE, calibTM, nu_eval)
    TE = load(calibTE, 'opl_poly', 'start_time', 'end_time',...
        'g_poly', 'g_mu', 'nu_min', 'nu_max');
    TM = load(calibTM, 'opl_poly', 'start_time', 'end_time',...
        'g_poly', 'g_mu', 'nu_min', 'nu_max');
    %
    g_TE = polyval(TE.g_poly, nu_eval, [], TE.g_mu);
    g_TM = polyval(TM.g_poly, nu_eval, [], TM.g_mu);
    % we have to recover the numerator of g, not just g. Fortunately we can
    % easily find the denominator.
    center_slope_TE = polyval(polyder(TE.opl_poly),(TE.start_time + TE.end_time)/2);
    center_slope_TM = polyval(polyder(TM.opl_poly),(TM.start_time + TM.end_time)/2);
    %center_deriv = polyval(polyder(result_polyfits(center_index,:)),crop_time);
    h = g_TM./g_TE * center_slope_TM/center_slope_TE;
end
