clear; clc;

R = 8.3144598;

P = 206 * 6894.757;
T = 288.75;

comp = [0.3094; 0.6906];
Tc = [305.322; 369.89];
Pc = [48.72e5; 42.48e5];
omega = [0.0995; 0.1521];
Zc = [0.2793; 0.2763];
MW = [30.069; 44.0956];
Vc = Zc .* R .* Tc ./ Pc * 1e6;
components = {'C2', 'C3'};
BIP = zeros(2);

V_exp = 83.573197;
rho_exp = 0.4757;
MW_mix = sum(comp .* MW);

%% PR EOS Methods
[~, Z_pr] = fugacitycoef_multicomp(comp, P, T, Pc, Tc, omega, BIP, 'PR');
if length(Z_pr) > 1, Z_pr = min(Z_pr); end
V_pr_cm3 = Z_pr * R * T / P * 1e6;

mix_pen_pr = peneloux_volume_shift(Pc, Tc, omega, comp, Vc, 'PR');
c_pen_pr = sum(comp .* mix_pen_pr.c_i);
V_pen_pr_cm3 = V_pr_cm3 - c_pen_pr;

[~, ~, mix_mt] = magoulas_tassios_volume_shift(T, Pc, Tc, omega, comp, Vc, components);
c_mt = mix_mt.c_mix * 1e6;
V_mt_cm3 = V_pr_cm3 - c_mt;

[~, ~, mix_ub] = ungerer_batut_volume_shift(T, Pc, Tc, omega, MW, comp, Vc);
c_ub = mix_ub.c_mix_cm3;
V_ub_cm3 = V_pr_cm3 - c_ub;

[~, ~, mix_baled_pr] = baled_volume_shift(T, Pc, Tc, omega, MW, components, comp, Vc, 'PR');
c_baled_pr = mix_baled_pr.c_mix_cm3;
V_baled_pr_cm3 = V_pr_cm3 - c_baled_pr;

[~, ~, mix_abu] = abudour_volume_shift(comp, P, T, Pc, Tc, omega, Vc, MW, components, BIP, 'PR');
c_abu = mix_abu.c_mix * 1e6;
V_abu_cm3 = mix_abu.V_VTEOS * 1e6;

%% SRK EOS Methods
[~, Z_srk] = fugacitycoef_multicomp(comp, P, T, Pc, Tc, omega, BIP, 'SRK');
if length(Z_srk) > 1, Z_srk = min(Z_srk); end
V_srk_cm3 = Z_srk * R * T / P * 1e6;

mix_pen_srk = peneloux_volume_shift(Pc, Tc, omega, comp, Vc, 'SRK');
c_pen_srk = sum(comp .* mix_pen_srk.c_i);
V_pen_srk_cm3 = V_srk_cm3 - c_pen_srk;

[~, ~, mix_cl] = chen_li_volume_shift(comp, P, T, Pc, Tc, omega, Zc, components, false);
c_cl = mix_cl.c_mix * 1e6;
V_cl_cm3 = mix_cl.V_VTEOS * 1e6;

[~, ~, mix_baled_srk] = baled_volume_shift(T, Pc, Tc, omega, MW, components, comp, Vc, 'SRK');
c_baled_srk = mix_baled_srk.c_mix_cm3;
V_baled_srk_cm3 = V_srk_cm3 - c_baled_srk;

mix_pm_srk = pina_martinez_volume_shift(T, Pc, Tc, omega, components,comp);
c_pm_srk = sum(comp .* mix_pm_srk.c_i);
V_pm_srk_cm3 = V_srk_cm3 - c_pm_srk;

%% Output Results
fprintf('======================================================================\n');
fprintf('VT MODEL VALIDATION: 206 psi C2/C3 Binary Mixture\n');
fprintf('======================================================================\n');
fprintf('Conditions: P = %.2f bar (206 psi), T = %.2f K\n', P/1e5, T);
fprintf('Composition: x_C2 = %.4f, x_C3 = %.4f\n', comp(1), comp(2));
fprintf('MW_mix = %.4f g/mol\n\n', MW_mix);

fprintf('%-25s %14s %12s %12s %12s\n', 'Method', 'V [cm3/mol]', 'rho [g/cm3]', 'c [cm3/mol]', 'Error [%]');
fprintf('%s\n', repmat('-', 1, 80));
fprintf('%-25s %14.4f %12.4f %12s %12s\n', 'Experimental', V_exp, rho_exp, '-', '-');
fprintf('%s\n', repmat('-', 1, 80));

fprintf('PR EOS METHODS:\n');
fprintf('%-25s %14.4f %12.4f %12s %12.2f\n', 'PR (untranslated)', V_pr_cm3, MW_mix/V_pr_cm3, '-', (V_pr_cm3-V_exp)/V_exp*100);
fprintf('%-25s %14.4f %12.4f %12.4f %12.2f\n', 'PR + Peneloux', V_pen_pr_cm3, MW_mix/V_pen_pr_cm3, c_pen_pr, (V_pen_pr_cm3-V_exp)/V_exp*100);
fprintf('%-25s %14.4f %12.4f %12.4f %12.2f\n', 'PR + Magoulas-Tassios', V_mt_cm3, MW_mix/V_mt_cm3, c_mt, (V_mt_cm3-V_exp)/V_exp*100);
fprintf('%-25s %14.4f %12.4f %12.4f %12.2f\n', 'PR + Ungerer-Batut', V_ub_cm3, MW_mix/V_ub_cm3, c_ub, (V_ub_cm3-V_exp)/V_exp*100);
fprintf('%-25s %14.4f %12.4f %12.4f %12.2f\n', 'PR + Baled', V_baled_pr_cm3, MW_mix/V_baled_pr_cm3, c_baled_pr, (V_baled_pr_cm3-V_exp)/V_exp*100);
fprintf('%-25s %14.4f %12.4f %12.4f %12.2f\n', 'PR + Abudour', V_abu_cm3, MW_mix/V_abu_cm3, c_abu, (V_abu_cm3-V_exp)/V_exp*100);
fprintf('%s\n', repmat('-', 1, 80));

fprintf('SRK EOS METHODS:\n');
fprintf('%-25s %14.4f %12.4f %12s %12.2f\n', 'SRK (untranslated)', V_srk_cm3, MW_mix/V_srk_cm3, '-', (V_srk_cm3-V_exp)/V_exp*100);
fprintf('%-25s %14.4f %12.4f %12.4f %12.2f\n', 'SRK + Peneloux', V_pen_srk_cm3, MW_mix/V_pen_srk_cm3, c_pen_srk, (V_pen_srk_cm3-V_exp)/V_exp*100);
fprintf('%-25s %14.4f %12.4f %12.4f %12.2f\n', 'SRK + Chen-Li', V_cl_cm3, MW_mix/V_cl_cm3, c_cl, (V_cl_cm3-V_exp)/V_exp*100);
fprintf('%-25s %14.4f %12.4f %12.4f %12.2f\n', 'SRK + Baled', V_baled_srk_cm3, MW_mix/V_baled_srk_cm3, c_baled_srk, (V_baled_srk_cm3-V_exp)/V_exp*100);
fprintf('%-25s %14.4f %12.4f %12.4f %12.2f\n', 'SRK + Pina-Martinez', V_pm_srk_cm3, MW_mix/V_pm_srk_cm3, c_pm_srk, (V_pm_srk_cm3-V_exp)/V_exp*100);
fprintf('%s\n', repmat('-', 1, 80));