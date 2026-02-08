function [GOC_depth, GOC_pressure, GOC_comp, GOC_temp] = detect_sgoc(model_type, comp_ref, P_ref, T_ref, h_ref, dTdh, depth_range, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params, tau)
%
% model_type : 0 = isothermal, 1 = haase, 2 = kempers, 3 = firoozabadi
% tau        : Firoozabadi parameter (default 4.0, only used for model 3)
%

if nargin < 17 || isempty(tau), tau = 4.0; end

scan_step = 2;
tol_depth = 0.1;
ssi_tol = 1e-10;
ssi_maxiter = 500;

comp_ref = comp_ref(:);
Pc = Pc(:); Tc = Tc(:); acentric = acentric(:); M_gmol = M_gmol(:);

if isnumeric(vt_method) && isscalar(vt_method) && vt_method >= 7
    eos_type = 'SRK';
else
    eos_type = 'PR';
end

model_names = {'Isothermal', 'Haase', 'Kempers', 'Firoozabadi'};

h_min = min(depth_range);
h_max = max(depth_range);

if h_ref <= h_min
    depths = h_ref:scan_step:h_max;
elseif h_ref >= h_max
    depths = h_ref:-scan_step:h_min;
else
    depths = [fliplr(h_ref:-scan_step:h_min), (h_ref+scan_step):scan_step:h_max];
end

GOC_depth = NaN;
GOC_pressure = NaN;
GOC_comp = comp_ref;
GOC_temp = T_ref;

stab_prev = NaN;
h_prev = NaN;

fprintf('=== GOC Detection: %s + SSI Stability ===\n', model_names{model_type+1});
fprintf('EOS = %s, VT = %d, ssi_tol = %.0e, ssi_maxiter = %d\n', eos_type, vt_method, ssi_tol, ssi_maxiter);
if model_type == 3, fprintf('tau = %.2f\n', tau); end
fprintf('%-8s %-10s %-10s %-10s %-12s %-8s\n', 'Depth', 'P [bar]', 'T [K]', 'Stability', 'TPD_min', 'Conv');
fprintf('%s\n', repmat('-', 1, 62));

for i = 1:length(depths)
    h = depths(i);

    try
        [comp_h, P_h, T_h] = call_grading(model_type, h, h_ref, comp_ref, P_ref, T_ref, dTdh, ...
            Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params, tau);
    catch ME
        fprintf('%-8.1f  ** grading failed: %s\n', h, ME.message);
        continue;
    end

    [stab, ~, ~, TPD_min, conv] = stability_analysis_ssi(comp_h, P_h, T_h, Pc, Tc, acentric, BIP, ssi_tol, ssi_maxiter, eos_type);

    fprintf('%-8.1f %-10.2f %-10.1f %-10d %-12.4e %-8d\n', h, P_h/1e5, T_h, stab, TPD_min, conv);

    if ~isnan(stab_prev) && stab_prev == 1 && stab == 2
        fprintf('\n*** Phase transition between h = %.1f (stable) and h = %.1f (unstable) ***\n', h_prev, h);

        [GOC_depth, GOC_pressure, GOC_comp, GOC_temp] = bisect_goc( ...
            model_type, h_prev, h, h_ref, comp_ref, P_ref, T_ref, dTdh, ...
            Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, ...
            vt_method, vt_params, tau, eos_type, ssi_tol, ssi_maxiter, tol_depth);

        fprintf('\n>>> GOC at h = %.2f m, P = %.2f bar, T = %.1f K\n', GOC_depth, GOC_pressure/1e5, GOC_temp);
        return;
    end

    stab_prev = stab;
    h_prev = h;
end

fprintf('\nNo GOC found in range [%.1f, %.1f] m.\n', h_min, h_max);

end


function [comp_h, P_h, T_h] = call_grading(model_type, h, h_ref, comp_ref, P_ref, T_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params, tau)

    switch model_type
        case 0
            [comp_h, P_h, ~, ~] = main(h, h_ref, comp_ref, P_ref, T_ref, ...
                Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params);
            T_h = T_ref;

        case 1
            [comp_h, P_h, T_h, ~, ~] = main_hasse(h, h_ref, comp_ref, P_ref, T_ref, dTdh, ...
                Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);

        case 2
            [comp_h, P_h, T_h, ~, ~] = main_kempers(h, h_ref, comp_ref, P_ref, T_ref, dTdh, ...
                Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);

        case 3
            [comp_h, P_h, T_h, ~, ~] = main_firoozabadi(h, h_ref, comp_ref, P_ref, T_ref, dTdh, ...
                Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, tau, vt_method, vt_params);

        otherwise
            error('model_type must be 0 (isothermal), 1 (haase), 2 (kempers), or 3 (firoozabadi)');
    end

end


function [h_goc, P_goc, comp_goc, T_goc] = bisect_goc(model_type, h_stable, h_unstable, h_ref, comp_ref, P_ref, T_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params, tau, eos_type, ssi_tol, ssi_maxiter, tol_depth)

h_lo = h_stable;
h_hi = h_unstable;

fprintf('\nBisection refinement:\n');

for iter = 1:50
    h_mid = (h_lo + h_hi) / 2;

    [comp_mid, P_mid, T_mid] = call_grading(model_type, h_mid, h_ref, comp_ref, P_ref, T_ref, dTdh, ...
        Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params, tau);

    [stab_mid, ~, ~, TPD_mid, ~] = stability_analysis_ssi(comp_mid, P_mid, T_mid, Pc, Tc, acentric, BIP, ssi_tol, ssi_maxiter, eos_type);

    fprintf('  Iter %2d: h = %8.3f m  [%.3f, %.3f]  stab = %d  TPD = %.4e\n', iter, h_mid, h_lo, h_hi, stab_mid, TPD_mid);

    if stab_mid == 1
        h_lo = h_mid;
    else
        h_hi = h_mid;
    end

    if abs(h_hi - h_lo) < tol_depth
        break;
    end
end

h_goc = (h_lo + h_hi) / 2;
[comp_goc, P_goc, T_goc] = call_grading(model_type, h_goc, h_ref, comp_ref, P_ref, T_ref, dTdh, ...
    Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params, tau);

end