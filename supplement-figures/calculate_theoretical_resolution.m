function resolution = calculate_theoretical_resolution(calib_file, nu_eval)
    load(calib_file, ...
        'opl_poly', 'start_time', 'end_time', ... % time domain scaling + range
        'g_poly', 'g_mu', 'nu_min', 'nu_max');
    %
    g = polyval(g_poly, nu_eval, [], g_mu);
    g_deriv = polyval(polyder(g_poly), nu_eval, [], g_mu)/g_mu(2);
    u_range = abs(polyval(opl_poly, start_time) - polyval(opl_poly, end_time));
    a = 2; % a value for the Hann window
    c0 = physconst('LightSpeed');
    resolution_freq = a*c0/u_range * (g + nu_eval.*g_deriv).^-1;
    resolution = c0*resolution_freq./(nu_eval.^2);
end
