clear; clc;

R = 8.314472;

P = 0.101352e6;
T = 20 + 273.15;

comp = 1;
Tc = 647.096;
Pc = 22.064e6;
omega = 0.3449;
Zc = 0.229;
MW = 18.015;
Vc = Zc * R * Tc / Pc * 1e6;
components = {'H2O'};
BIP = 0;

rho_exp = 998.29;
V_exp = MW / rho_exp * 1000;


%% SRK EOS (Pina-Martinez Twu91 alpha with water-specific L,M,N)


[~, Z_srk] = fugacitycoef_multicomp(comp, P, T, Pc, Tc, omega, BIP, 'SRK');
if length(Z_srk) > 1, Z_srk = min(Z_srk); end
V_srk_cm3 = Z_srk * R * T / P * 1e6;
rho_srk = MW / V_srk_cm3 * 1000;

% SRK + Peneloux VT
mix_pen_srk = peneloux_volume_shift(Pc, Tc, omega, comp, Vc, 'SRK');
c_pen_srk = mix_pen_srk.c_i;
V_pen_srk_cm3 = V_srk_cm3 - c_pen_srk;
rho_pen_srk = MW / V_pen_srk_cm3 * 1000;

% SRK + Pina-Martinez VT (DIPPR table: c = 8.967 cm3/mol)
mix_pm_srk = pina_martinez_volume_shift(T, Pc, Tc, omega, components,comp);
c_pm_srk = mix_pm_srk.c_i;
V_pm_srk_cm3 = V_srk_cm3 - c_pm_srk;
rho_pm_srk = MW / V_pm_srk_cm3 * 1000;

% SRK + Chen-Li VT (table: c1=0.02425, c2=1.30564, c3=2.17549)
[~, ~, mix_cl] = chen_li_volume_shift(comp, P, T, Pc, Tc, omega, Zc, components, false);
c_cl = mix_cl.c_mix * 1e6;
V_cl_cm3 = mix_cl.V_VTEOS * 1e6;
rho_cl = MW / V_cl_cm3 * 1000;

% SRK + Baled VT (hydrocarbon correlation - not for water)
[~, ~, mix_baled_srk] = baled_volume_shift(T, Pc, Tc, omega, MW, components, comp, Vc, 'SRK');
c_baled_srk = mix_baled_srk.c_mix_cm3;
V_baled_srk_cm3 = V_srk_cm3 - c_baled_srk;
rho_baled_srk = MW / V_baled_srk_cm3 * 1000;

%% PR EOS (standard alpha)


[~, Z_pr] = fugacitycoef_multicomp(comp, P, T, Pc, Tc, omega, BIP, 'PR');
if length(Z_pr) > 1, Z_pr = min(Z_pr); end
V_pr_cm3 = Z_pr * R * T / P * 1e6;
rho_pr = MW / V_pr_cm3 * 1000;

% PR + Peneloux VT
mix_pen_pr = peneloux_volume_shift(Pc, Tc, omega, comp, Vc, 'PR');
c_pen_pr = mix_pen_pr.c_i;
V_pen_pr_cm3 = V_pr_cm3 - c_pen_pr;
rho_pen_pr = MW / V_pen_pr_cm3 * 1000;

% PR + Magoulas-Tassios VT
[~, ~, mix_mt] = magoulas_tassios_volume_shift(T, Pc, Tc, omega, comp, Vc, components);
c_mt = mix_mt.c_mix * 1e6;
V_mt_cm3 = V_pr_cm3 - c_mt;
rho_mt = MW / V_mt_cm3 * 1000;

% PR + Ungerer-Batut VT (hydrocarbon correlation)
[~, ~, mix_ub] = ungerer_batut_volume_shift(T, Pc, Tc, omega, MW, comp, Vc);
c_ub = mix_ub.c_mix_cm3;
V_ub_cm3 = V_pr_cm3 - c_ub;
rho_ub = MW / V_ub_cm3 * 1000;

% PR + Baled VT (hydrocarbon correlation)
[~, ~, mix_baled_pr] = baled_volume_shift(T, Pc, Tc, omega, MW, components, comp, Vc, 'PR');
c_baled_pr = mix_baled_pr.c_mix_cm3;
V_baled_pr_cm3 = V_pr_cm3 - c_baled_pr;
rho_baled_pr = MW / V_baled_pr_cm3 * 1000;

% PR + Abudour VT
[~, ~, mix_abu] = abudour_volume_shift(comp, P, T, Pc, Tc, omega, Vc, MW, components, BIP, 'PR');
c_abu = mix_abu.c_mix * 1e6;
V_abu_cm3 = mix_abu.V_VTEOS * 1e6;
rho_abu = MW / V_abu_cm3 * 1000;

%% OUTPUT RESULTS


fprintf('%-35s %10s %14s %12s %10s\n', 'Method', 'c [cm3/mol]', 'V [cm3/mol]', 'rho [kg/m3]', 'Error [%]');
fprintf('%s\n', repmat('=', 1, 90));
fprintf('%-35s %10s %14.4f %12.2f %10s\n', 'Experimental', '-', V_exp, rho_exp, '-');
fprintf('%s\n', repmat('=', 1, 90));

fprintf('\nSRK EOS (Twu91 alpha, water L,M,N):\n');
fprintf('%s\n', repmat('-', 1, 90));
fprintf('%-35s %10s %14.4f %12.2f %10.2f\n', 'SRK (no VT)', '-', V_srk_cm3, rho_srk, (rho_srk-rho_exp)/rho_exp*100);
fprintf('%-35s %10.4f %14.4f %12.2f %10.2f\n', 'SRK + Peneloux', c_pen_srk, V_pen_srk_cm3, rho_pen_srk, (rho_pen_srk-rho_exp)/rho_exp*100);
fprintf('%-35s %10.4f %14.4f %12.2f %10.2f\n', 'SRK + Pina-Martinez', c_pm_srk, V_pm_srk_cm3, rho_pm_srk, (rho_pm_srk-rho_exp)/rho_exp*100);
fprintf('%-35s %10.4f %14.4f %12.2f %10.2f\n', 'SRK + Chen-Li', c_cl, V_cl_cm3, rho_cl, (rho_cl-rho_exp)/rho_exp*100);
fprintf('%-35s %10.4f %14.4f %12.2f %10.2f\n', 'SRK + Baled*', c_baled_srk, V_baled_srk_cm3, rho_baled_srk, (rho_baled_srk-rho_exp)/rho_exp*100);

fprintf('\nPR EOS (standard alpha):\n');
fprintf('%s\n', repmat('-', 1, 90));
fprintf('%-35s %10s %14.4f %12.2f %10.2f\n', 'PR (no VT)', '-', V_pr_cm3, rho_pr, (rho_pr-rho_exp)/rho_exp*100);
fprintf('%-35s %10.4f %14.4f %12.2f %10.2f\n', 'PR + Peneloux', c_pen_pr, V_pen_pr_cm3, rho_pen_pr, (rho_pen_pr-rho_exp)/rho_exp*100);
fprintf('%-35s %10.4f %14.4f %12.2f %10.2f\n', 'PR + Magoulas-Tassios', c_mt, V_mt_cm3, rho_mt, (rho_mt-rho_exp)/rho_exp*100);
fprintf('%-35s %10.4f %14.4f %12.2f %10.2f\n', 'PR + Ungerer-Batut*', c_ub, V_ub_cm3, rho_ub, (rho_ub-rho_exp)/rho_exp*100);
fprintf('%-35s %10.4f %14.4f %12.2f %10.2f\n', 'PR + Baled*', c_baled_pr, V_baled_pr_cm3, rho_baled_pr, (rho_baled_pr-rho_exp)/rho_exp*100);
fprintf('%-35s %10.4f %14.4f %12.2f %10.2f\n', 'PR + Abudour', c_abu, V_abu_cm3, rho_abu, (rho_abu-rho_exp)/rho_exp*100);

