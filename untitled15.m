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
for i = 10:n
    BIP(1,i) = 0.0800; BIP(i,1) = 0.0800;
end
for i = 3:10
    BIP(2,i) = 0.1200; BIP(i,2) = 0.1200;
end
for i = 11:n
    BIP(2,i) = 0.100; BIP(i,2) = 0.100;
end

Vc = [89.8; 94.0; 99.0; 148.0; 203.0; 263.0; 255.0; 306.0; 304.0; 370.0; ...
      432; 492; 548; 600.22; 697.08; 810.90; 910.91; 980.42; 1061.51; ...
      1152.43; 1287.26; 1648.29];
Zc = Pc .* (Vc * 1e-6) ./ (R .* Tc);

h_ref = 175;
press_ref = 284e5;
temp_ref = 93 + 273.15;
dTdh = 0.025;
tau = 4.0;

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
   115.26,      3.795,   -1.48e-3,    1.86e-7
];

H_ig_ref = [8330.789; 19459.101; 2.642; 9761.134; 19519.622; 29278.121; 29278.121; 39036.609; 39036.609; 48795.1026; ...
            58553.595; 68312.084; 78069.882; 86301.327; 105418.986; 131284.869; 158312.556; 180345.162; 209390.365; ...
            246519.548; 306655.264; 434614.185];

vt_params.c_custom = c_JY;
vt_params.Vc = Vc;
vt_params.components = components;
vt_params.Zc = Zc;
vt_method = 0;

tol = 1e-10;
maxiter = 500;

% GOC_depth = 111.84 ; 

[GOC_depth, GOC_pressure, GOC_comp, GOC_temp] = detect_sgoc(3, comp_ref, press_ref, temp_ref, h_ref, dTdh, [80 175], Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
 
%% Oil zone
h_oil = [GOC_depth,120 , 175, 204, 228, 327];
n_oil = length(h_oil);
oil_P    = zeros(n_oil, 1);
oil_Pbub = zeros(n_oil, 1);
oil_T    = zeros(n_oil, 1);
oil_comp = zeros(n_oil, n);
oil_rho  = zeros(n_oil, 1);
oil_GOR  = zeros(n_oil, 1);
oil_Bo   = zeros(n_oil, 1);

for i = 1:n_oil
    fprintf('Oil point %d/%d: h = %.1f m\n', i, n_oil, h_oil(i));
    [comp_h, P_h, T_h, Pb, ~] = main_firoozabadi(h_oil(i), h_ref, comp_ref, press_ref, temp_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, tau, vt_method, vt_params);
    oil_comp(i,:) = comp_h';
    oil_P(i)    = P_h;
    oil_T(i)    = T_h;
    oil_Pbub(i) = Pb;
    try
        oil_rho(i) = calculate_density(comp_h, P_h, T_h, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params);
    catch
        oil_rho(i) = NaN;
    end
    try
        [oil_GOR(i), oil_Bo(i)] = calculate_GOR_STO(comp_h, P_h, T_h, Pc, Tc, acentric, BIP, M_gmol, 'vt_method', vt_method, 'vt_params', vt_params);
    catch
        oil_GOR(i) = NaN;
        oil_Bo(i)  = NaN;
    end
end

%% GOC -> vapor composition
[press_bub, comp_gas_ref] = pressbub_multicomp_newton(oil_comp(1,:)', oil_P(1), oil_T(1), Pc, Tc, acentric, BIP, tol, maxiter);

%% Gas zone
h_gas = [GOC_depth,110, 100, 80, 60, 40, 20, 0];
h_gas = h_gas(h_gas <= GOC_depth);
n_gas = length(h_gas);
gas_P    = zeros(n_gas, 1);
gas_Pdew = zeros(n_gas, 1);
gas_rho  = zeros(n_gas, 1);
gas_GOR  = NaN(n_gas, 1);
gas_Bo   = NaN(n_gas, 1);

for i = 1:n_gas
    fprintf('Gas point %d/%d: h = %.1f m\n', i, n_gas, h_gas(i));
    [comp_h, P_h, T_h, ~, Pd] = main_firoozabadi(h_gas(i), GOC_depth, comp_gas_ref, oil_P(1), oil_T(1), dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, tau, vt_method, vt_params);
    gas_P(i)    = P_h;
    gas_Pdew(i) = Pd;
    try
        gas_rho(i) = calculate_density(comp_h, P_h, T_h, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params);
    catch
        gas_rho(i) = NaN;
    end
end


all_h    = [flipud(h_gas'); h_oil'];
all_P    = [flipud(gas_P);  oil_P];
all_Psat = [flipud(gas_Pdew); oil_Pbub];
all_rho  = [flipud(gas_rho);  oil_rho];
all_GOR  = [flipud(gas_GOR);  oil_GOR];
all_Bo   = [flipud(gas_Bo);   oil_Bo];

%% Experimental
exp_h    = [0; 175; 204; 228; 327];
exp_P    = [279; 284; 286; 287; 293];
exp_Psat = [270; 272; 267; 265; 242];

%% Plot 1: Pressure & Psat
figure('Position', [50 100 500 450], 'Color', 'w'); hold on;
plot(all_P/1e5, all_h, 'b-', 'LineWidth', 2);
plot(all_Psat/1e5, all_h, 'b--', 'LineWidth', 2);
plot(exp_P, exp_h, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
plot(exp_Psat, exp_h, 'k^', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
yline(GOC_depth, 'r:', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse', 'FontSize', 11);
xlabel('Pressure [bar]'); ylabel('Depth [m]');
title(sprintf('Reservoir 1 — Firoozabadi (\\tau = %.1f), No VT', tau));
legend('P_{res}', 'P_{sat}', 'Exp P_{res}', 'Exp P_{sat}', 'GOC', 'Location', 'southeast');
grid on;

%% Plot 2: Density
figure('Position', [600 100 500 450], 'Color', 'w'); hold on;
plot(all_rho, all_h, 'b-o', 'LineWidth', 2, 'MarkerSize', 5, 'MarkerFaceColor', 'b');
yline(GOC_depth, 'r:', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse', 'FontSize', 11);
xlabel('Density [kg/m^3]'); ylabel('Depth [m]');
title(sprintf('Reservoir 1 — Firoozabadi (\\tau = %.1f), Density', tau));
legend('Density', 'GOC', 'Location', 'best');
grid on;

%% Plot 3: GOR (oil zone only)
% figure('Position', [50 600 500 450], 'Color', 'w'); hold on;
% plot(all_GOR, all_h, 'b-o', 'LineWidth', 2, 'MarkerSize', 5, 'MarkerFaceColor', 'b');
% yline(GOC_depth, 'r:', 'LineWidth', 1.5);
% set(gca, 'YDir', 'reverse', 'FontSize', 11);
% xlabel('GOR [Sm^3/Sm^3]'); ylabel('Depth [m]');
% title(sprintf('Reservoir 1 — Firoozabadi (\\tau = %.1f), GOR', tau));
% legend('GOR', 'GOC', 'Location', 'best');
% grid on;

