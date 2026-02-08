clear; clc; close all;

components = {'N2', 'CO2', 'C1', 'C2', 'C3', 'iC4', 'nC4', 'iC5', 'nC5', ...
              'C6', 'C7', 'C8', 'C9', 'C10-C11', 'C12-C13', 'C14-C16', ...
              'C17-C18', 'C19-C21', 'C22-C24', 'C25-C29', 'C30-C37', 'C38-C80'};

comp_ref = [0.42; 0.69; 50.03; 7.848; 6.769; 1.04; 3.199; 1.16; 1.55; 1.88; ...
            3.499; 3.749; 2.28; 3.259; 2.589; 2.929; 1.46; 1.65; 1.17; 1.24; ...
            0.96; 0.63];
comp_ref = comp_ref / sum(comp_ref);

M_gmol = [28.0; 44.0; 16.0; 30.1; 44.1; 58.1; 58.1; 72.2; 72.2; 86.2; ...
          96.0; 107.0; 121.0; 140.1; 167.6; 204.7; 243.6; 275.3; 317.0; ...
          370.4; 456.8; 640.9];

Tc = [-147.0; 31.1; -82.5; 32.2; 96.6; 134.9; 152.1; 187.2; 196.4; ...
      234.2; 269.9; 290.9; 315.3; 345.9; 385.0; 432.7; 476.8; 511.1; ...
      553.3; 604.9; 683.3; 847.8] + 273.15;

Pc = [33.94; 73.76; 46.00; 48.84; 42.46; 36.48; 38.00; 33.84; 33.74; ...
      29.69; 29.60; 27.86; 25.83; 23.76; 21.58; 19.59; 18.24; 17.54; ...
      16.84; 16.23; 15.62; 14.72] * 1e5;

acentric = [0.040; 0.225; 0.008; 0.098; 0.152; 0.176; 0.193; 0.227; 0.251; ...
            0.296; 0.338; 0.374; 0.420; 0.483; 0.570; 0.685; 0.795; 0.881; ...
            0.984; 1.098; 1.229; 1.159];

c_JY = [-4.23; -1.91; -5.20; -5.79; -6.35; -7.18; -6.49; -6.20; -5.12; ...
        1.42; 8.45; 9.20; 10.32; 11.12; 10.48; 6.57; -1.30; -9.89; ...
        -23.43; -43.30; -79.96; -151.44];

R = 8.3144598;
n = length(comp_ref);
BIP = zeros(n, n);
BIP(1,2) = -0.0170; BIP(2,1) = -0.0170;
BIP(1,3) = 0.0311;  BIP(3,1) = 0.0311;
BIP(1,4) = 0.0515;  BIP(4,1) = 0.0515;
BIP(1,5) = 0.0852;  BIP(5,1) = 0.0852;
BIP(1,6) = 0.1033;  BIP(6,1) = 0.1033;
BIP(1,7) = 0.0800;  BIP(7,1) = 0.0800;
BIP(1,8) = 0.0922;  BIP(8,1) = 0.0922;
BIP(1,9) = 0.1000;  BIP(9,1) = 0.1000;
for i = 10:n, BIP(1,i) = 0.0800; BIP(i,1) = 0.0800; end
for i = 3:10, BIP(2,i) = 0.1200; BIP(i,2) = 0.1200; end
for i = 11:n, BIP(2,i) = 0.100;  BIP(i,2) = 0.100;  end

Vc = [89.8; 94.0; 99.0; 148.0; 203.0; 263.0; 255.0; 306.0; 304.0; 370.0; ...
      432; 492; 548; 600.22; 697.08; 810.90; 910.91; 980.42; 1061.51; ...
      1152.43; 1287.26; 1648.29];
Zc = Pc .* (Vc * 1e-6) ./ (R .* Tc);

h_ref = 175;
press_ref = 284e5;
temp_ref = 93 + 273.15;
dTdh = 0.025;
tol = 1e-10;
maxiter = 500;

Cp_coeffs = [
    31.15,     -0.014,    2.68e-5,   -1.17e-8;
    19.79,      0.073,   -5.60e-5,    1.72e-8;
    19.25,      0.052,    1.20e-5,   -1.13e-8;
     5.41,      0.178,   -6.94e-5,    8.71e-9;
    -4.22,      0.306,   -1.59e-4,    3.21e-8;
    -1.39,      0.385,   -1.83e-4,    2.90e-8;
     9.49,      0.331,   -1.11e-4,   -2.82e-9;
    -9.52,      0.507,   -2.73e-4,    5.72e-8;
    -3.63,      0.487,   -2.58e-4,    5.30e-8;
    -4.41,      0.582,   -3.12e-4,    6.49e-8;
    -5.15,      0.676,   -3.65e-4,    7.66e-8;
    -6.10,      0.771,   -4.20e-4,    8.85e-8;
     3.14,      0.677,   -1.93e-4,   -2.98e-8;
    25.20,      0.830,   -3.23e-4,    4.06e-8;
    30.14,      0.993,   -3.87e-4,    4.86e-8;
    36.83,      1.213,   -4.72e-4,    5.93e-8;
    43.82,      1.443,   -5.62e-4,    7.06e-8;
    49.52,      1.630,   -6.35e-4,    7.98e-8;
    57.03,      1.878,   -7.31e-4,    9.19e-8;
    66.63,      2.194,   -8.55e-4,    1.07e-7;
    82.18,      2.706,   -1.05e-3,    1.32e-7;
   115.26,      3.795,   -1.48e-3,    1.86e-7];

H_ig_ref = [8330.789; 19459.101; 2.642; 9761.134; 19519.622; 29278.121; 29278.121; 39036.609; 39036.609; 48795.1026; ...
            58553.595; 68312.084; 78069.882; 86301.327; 105418.986; 131284.869; 158312.556; 180345.162; 209390.365; ...
            246519.548; 306655.264; 434614.185];

idx_C1  = 3;
idx_C7p = 11:22;

h_oil_base = [175, 204, 228, 327];
h_gas_base = [120, 100, 80, 60, 40, 20, 0];

vt_ids    = [0, 1, 2, 3, 4, 5, 6];
vt_names  = {'No VT', 'Peneloux', 'Magoulas-Tassios', 'Ungerer-Batut', 'Baled', 'Abudour', 'JY Tuned'};
n_vt = length(vt_ids);

results = struct();

for m = 1:n_vt
    vt_method = vt_ids(m);
    fprintf('\n=== Running Haase + %s (VT=%d) ===\n', vt_names{m}, vt_method);

    vt_params = struct();
    vt_params.Vc = Vc;
    vt_params.components = components;
    vt_params.Zc = Zc;
    if vt_method == 6
        vt_params.c_custom = c_JY;
    end

    [GOC_h, ~, ~, GOC_T] = detect_sgoc(1, comp_ref, press_ref, temp_ref, h_ref, dTdh, [80 175], Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
    fprintf('  GOC = %.2f m\n', GOC_h);

    h_oil = [GOC_h, GOC_h+2, h_oil_base];
    n_oil = length(h_oil);
    oil_P = zeros(n_oil,1); oil_Pbub = zeros(n_oil,1);
    oil_rho = zeros(n_oil,1); oil_comp = zeros(n_oil, n);

    for i = 1:n_oil
        [ch,Ph,Th,Pb,~] = main_hasse(h_oil(i), h_ref, comp_ref, press_ref, temp_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
        oil_P(i) = Ph; oil_Pbub(i) = Pb; oil_comp(i,:) = ch';
        try, oil_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params); catch, oil_rho(i) = NaN; end
    end

    [~, comp_gas_ref] = pressbub_multicomp_newton(oil_comp(1,:)', oil_P(1), GOC_T, Pc, Tc, acentric, BIP, tol, maxiter);

    h_gas = [GOC_h, h_gas_base(h_gas_base < GOC_h)];
    n_gas = length(h_gas);
    gas_P = zeros(n_gas,1); gas_Pdew = zeros(n_gas,1);
    gas_rho = zeros(n_gas,1); gas_comp = zeros(n_gas, n);

    for i = 1:n_gas
        [ch,Ph,Th,~,Pd] = main_hasse(h_gas(i), GOC_h, comp_gas_ref, oil_P(1), GOC_T, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
        gas_P(i) = Ph; gas_Pdew(i) = Pd; gas_comp(i,:) = ch';
        try, gas_rho(i) = calculate_density(ch, Ph, Th, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params); catch, gas_rho(i) = NaN; end
    end

    all_h    = [flipud(h_gas'); h_oil'];
    all_comp = [flipud(gas_comp); oil_comp];

    results(m).name = vt_names{m};
    results(m).GOC  = GOC_h;
    results(m).h    = all_h;
    results(m).P    = [flipud(gas_P);  oil_P];
    results(m).Psat = [flipud(gas_Pdew); oil_Pbub];
    results(m).rho  = [flipud(gas_rho); oil_rho];
    results(m).C1   = all_comp(:, idx_C1) * 100;
    results(m).C7p  = sum(all_comp(:, idx_C7p), 2) * 100;
end

%% Experimental (Pedersen 2015, Table 1)
exp_h    = [0; 175; 204; 228; 327];
exp_P    = [279; 284; 286; 287; 293];
exp_Psat = [270; 272; 267; 265; 242];
exp_C1   = [75.66; 50.04; 49.88; 48.89; 45.66];
exp_C7p  = [1.29+1.06+0.54+1.57; 3.50+3.75+2.28+15.88; 3.49+3.73+2.31+16.11; 3.62+3.85+2.36+16.70; 4.01+4.28+2.58+17.66];

%% ========================================================================
%  PUBLICATION FIGURE
%  ========================================================================

colors = {[0.00 0.45 0.74], ...   % blue     - No VT
          [0.85 0.33 0.10], ...   % red      - Peneloux
          [0.47 0.67 0.19], ...   % green    - Magoulas-Tassios
          [0.49 0.18 0.56], ...   % purple   - Ungerer-Batut
          [0.93 0.69 0.13], ...   % gold     - Baled
          [0.30 0.75 0.93], ...   % cyan     - Abudour
          [0.64 0.08 0.18]};      % crimson  - JY Tuned

lstyles = {'-', '--', '-.', ':', '-', '--', '-.'};
lw      = 1.6;
lw_sat  = 1.2;
fs      = 9;
fs_lab  = 10;
fs_pan  = 11;

fig = figure('Units', 'centimeters', 'Position', [2 2 24 20], ...
    'Color', 'w', 'PaperPositionMode', 'auto');

fmt_ax = @(ax) set(ax, 'YDir', 'reverse', 'FontSize', fs, 'FontName', 'Helvetica', ...
    'TickDir', 'in', 'LineWidth', 0.6, 'TickLength', [0.02 0.02], 'Box', 'on');

% ---- (a) Pressure & Psat ----
ax1 = subplot(2,2,1); hold on;
for m = 1:n_vt
    plot(results(m).P/1e5, results(m).h, lstyles{m}, 'Color', colors{m}, 'LineWidth', lw);
    v = ~isnan(results(m).Psat) & results(m).Psat > 0;
    plot(results(m).Psat(v)/1e5, results(m).h(v), lstyles{m}, 'Color', colors{m}, 'LineWidth', lw_sat);
end
plot(exp_P, exp_h, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
plot(exp_Psat, exp_h, 'k^', 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'LineWidth', 0.8);
fmt_ax(gca);
xlabel('Pressure (bar)', 'FontSize', fs_lab);
ylabel('Depth (m)', 'FontSize', fs_lab);
text(0.03, 0.03, '(a)', 'Units', 'normalized', 'FontSize', fs_pan, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');
xlim([215 305]);

% ---- (b) Density ----
ax2 = subplot(2,2,2); hold on;
for m = 1:n_vt
    v = ~isnan(results(m).rho) & results(m).rho > 0;
    plot(results(m).rho(v), results(m).h(v), lstyles{m}, 'Color', colors{m}, 'LineWidth', lw);
end
fmt_ax(gca);
xlabel('Density (kg/m^{3})', 'FontSize', fs_lab);
ylabel('Depth (m)', 'FontSize', fs_lab);
text(0.03, 0.03, '(b)', 'Units', 'normalized', 'FontSize', fs_pan, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');

% ---- (c) C1 mol% ----
ax3 = subplot(2,2,3); hold on;
for m = 1:n_vt
    plot(results(m).C1, results(m).h, lstyles{m}, 'Color', colors{m}, 'LineWidth', lw);
end
plot(exp_C1, exp_h, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
fmt_ax(gca);
xlabel('C_{1} (mol%)', 'FontSize', fs_lab);
ylabel('Depth (m)', 'FontSize', fs_lab);
text(0.03, 0.03, '(c)', 'Units', 'normalized', 'FontSize', fs_pan, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');

% ---- (d) C7+ mol% ----
ax4 = subplot(2,2,4); hold on;
for m = 1:n_vt
    plot(results(m).C7p, results(m).h, lstyles{m}, 'Color', colors{m}, 'LineWidth', lw);
end
plot(exp_C7p, exp_h, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
fmt_ax(gca);
xlabel('C_{7+} (mol%)', 'FontSize', fs_lab);
ylabel('Depth (m)', 'FontSize', fs_lab);
text(0.03, 0.03, '(d)', 'Units', 'normalized', 'FontSize', fs_pan, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');

% ---- Field GOC on all panels ----
GOC_field = 140;
for ax = [ax1, ax2, ax3, ax4]
    axes(ax);
    yline(GOC_field, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
end
axes(ax1);
xl1 = xlim;
text(xl1(1) + 0.02*(xl1(2)-xl1(1)), GOC_field - 4, 'Field GOC', ...
    'FontSize', 7, 'FontName', 'Helvetica', 'Color', [0.4 0.4 0.4]);

% ---- Shared legend ----
h_leg = gobjects(1, n_vt + 3);
axes(ax1);
for m = 1:n_vt
    h_leg(m) = plot(NaN, NaN, lstyles{m}, 'Color', colors{m}, 'LineWidth', lw);
end
h_leg(n_vt+1) = plot(NaN, NaN, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
h_leg(n_vt+2) = plot(NaN, NaN, 'k^', 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'LineWidth', 0.8);
h_leg(n_vt+3) = plot(NaN, NaN, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);

leg_str = cell(1, n_vt + 3);
for m = 1:n_vt
    leg_str{m} = sprintf('%s (GOC %.0f m)', results(m).name, results(m).GOC);
end
leg_str{n_vt+1} = 'Measured';
leg_str{n_vt+2} = 'Measured P_{sat}';
leg_str{n_vt+3} = 'Field GOC (140 m)';

lgd = legend(h_leg, leg_str, 'Orientation', 'horizontal', ...
    'FontSize', 6.5, 'NumColumns', 4, 'Box', 'off', 'FontName', 'Helvetica');
lgd.Position = [0.02 0.001 0.96 0.06];

annotation(fig, 'textbox', [0.10 0.96 0.80 0.04], ...
    'String', 'Reservoir 1 — Haase Model — VT Method Comparison (PR-EOS)', ...
    'FontSize', 11, 'FontWeight', 'bold', 'FontName', 'Helvetica', ...
    'HorizontalAlignment', 'center', 'EdgeColor', 'none');

