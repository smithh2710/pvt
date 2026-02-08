function [GOC_depth, GOC_type, GOC_pressure, GOC_comp, info] = detect_GOC(comp_ref, press_ref, temp_ref, h_ref, depth_range, Pc, Tc, acentric, BIP, M_gmol, varargin)
% DETECT_GOC - Gas-Oil Contact Detection with Multiple Grading Models
% =========================================================================
% Comprehensive GOC detection supporting isothermal and non-isothermal
% compositional grading models. Based on Nikpoor (2014) and Hoier & Whitson
% (2000) methodology.
%
% GOC TYPES:
%   'saturated'     : P_res crosses P_sat (phase split occurs)
%   'undersaturated': P_bub = P_dew (critical transition, no phase split)
%   'none'          : No GOC found in depth range
%
% GRADING MODELS:
%   'isothermal'    : main.m - Isothermal compositional grading
%   'hasse'         : main_hasse.m - Haase/Kusochkova non-isothermal
%   'kempers'       : main_kempers.m - Kempers (1989) thermal diffusion
%   'firoozabadi'   : main_firoozabadi.m - Shukla-Firoozabadi (2000)
%
% =========================================================================
% VOLUME TRANSLATION METHODS (vt_method):
%
%   PR-EOS Based:
%     0  - PR, no volume translation
%     1  - PR + Peneloux (1982)
%     2  - PR + Magoulas-Tassios (1990)
%     3  - PR + Ungerer-Batut (1997)
%     4  - PR + Baled (2012)
%     5  - PR + Abudour (2012)
%     6  - PR + Custom (user-provided c values)
%
%   SRK-EOS Based:
%     7  - SRK, no volume translation
%     8  - SRK + Pina-Martinez (2018)
%     9  - SRK + Chen-Li (2020)
%     10 - SRK + Baled (2012)
%
% =========================================================================
% INPUTS (Required):
%   comp_ref    : Reference composition [mole fractions] (column vector)
%   press_ref   : Reference pressure [Pa]
%   temp_ref    : Reference temperature [K]
%   h_ref       : Reference depth [m]
%   depth_range : [h_min, h_max] - Depth range to search [m]
%   Pc          : Critical pressures [Pa]
%   Tc          : Critical temperatures [K]
%   acentric    : Acentric factors [-]
%   BIP         : Binary interaction parameter matrix
%   M_gmol      : Molecular weights [g/mol]
%
% INPUTS (Optional Name-Value Pairs):
%   'grading_model' : 'isothermal' (default), 'hasse', 'kempers', 'firoozabadi'
%   'dTdh'          : Temperature gradient [K/m] (required for non-isothermal)
%   'Cp_coeffs'     : Ideal gas Cp coefficients (required for hasse/kempers)
%   'H_ig_ref'      : Reference ideal gas enthalpies [J/mol]
%   'tau'           : Firoozabadi tau parameter (default: 4.0)
%   'vt_method'     : Volume translation method (0-10, default: 0)
%   'vt_params'     : VT parameters structure
%   'scan_step'     : Coarse scan step size [m] (default: 5)
%   'tol_depth'     : Depth convergence tolerance [m] (default: 0.5)
%   'verbose'       : Print detailed output (default: true)
%
% OUTPUTS:
%   GOC_depth   : Depth of GOC [m] (NaN if not found)
%   GOC_type    : 'saturated', 'undersaturated', or 'none'
%   GOC_pressure: Pressure at GOC [Pa]
%   GOC_comp    : Composition at GOC [mole fractions]
%   info        : Diagnostic structure with scan data
%
% =========================================================================
% USAGE EXAMPLES:
%
%   % Isothermal with PR-EOS (simplest)
%   [GOC_h, type, GOC_P, GOC_z, info] = detect_GOC(comp, P, T, h_ref, [0 200], ...
%       Pc, Tc, w, BIP, M);
%
%   % Isothermal with VT (PR + Abudour)
%   [GOC_h, type, GOC_P, GOC_z, info] = detect_GOC(comp, P, T, h_ref, [0 200], ...
%       Pc, Tc, w, BIP, M, 'vt_method', 5, 'vt_params', vt_params);
%
%   % SRK + Chen-Li
%   [GOC_h, type, GOC_P, GOC_z, info] = detect_GOC(comp, P, T, h_ref, [0 200], ...
%       Pc, Tc, w, BIP, M, 'vt_method', 9, 'vt_params', vt_params);
%
%   % Hasse non-isothermal
%   [GOC_h, type, GOC_P, GOC_z, info] = detect_GOC(comp, P, T, h_ref, [0 200], ...
%       Pc, Tc, w, BIP, M, 'grading_model', 'hasse', 'dTdh', 0.025, ...
%       'Cp_coeffs', Cp, 'H_ig_ref', H_ig);
%
%   % Kempers with VT
%   [GOC_h, type, GOC_P, GOC_z, info] = detect_GOC(comp, P, T, h_ref, [0 200], ...
%       Pc, Tc, w, BIP, M, 'grading_model', 'kempers', 'dTdh', 0.025, ...
%       'Cp_coeffs', Cp, 'H_ig_ref', H_ig, 'vt_method', 6, 'vt_params', vt_params);
%
%   % Firoozabadi with custom tau
%   [GOC_h, type, GOC_P, GOC_z, info] = detect_GOC(comp, P, T, h_ref, [0 200], ...
%       Pc, Tc, w, BIP, M, 'grading_model', 'firoozabadi', 'dTdh', 0.025, ...
%       'Cp_coeffs', Cp, 'H_ig_ref', H_ig, 'tau', 5.0);
%
% =========================================================================
% REFERENCES:
%   Nikpoor, M.H. (2014). MSc Thesis, University of Calgary, Section 4.4
%   Hoier, L. & Whitson, C.H. (2000). SPE 63085
%   Pedersen, K.S. et al. (2015). SPE-175085-MS
% =========================================================================

    %% Parse inputs
    p = inputParser;
    p.FunctionName = 'detect_GOC';
    
    addRequired(p, 'comp_ref', @isnumeric);
    addRequired(p, 'press_ref', @isnumeric);
    addRequired(p, 'temp_ref', @isnumeric);
    addRequired(p, 'h_ref', @isnumeric);
    addRequired(p, 'depth_range', @isnumeric);
    addRequired(p, 'Pc', @isnumeric);
    addRequired(p, 'Tc', @isnumeric);
    addRequired(p, 'acentric', @isnumeric);
    addRequired(p, 'BIP', @isnumeric);
    addRequired(p, 'M_gmol', @isnumeric);
    
    addParameter(p, 'grading_model', 'isothermal', @ischar);
    addParameter(p, 'dTdh', 0, @isnumeric);
    addParameter(p, 'Cp_coeffs', [], @isnumeric);
    addParameter(p, 'H_ig_ref', [], @isnumeric);
    addParameter(p, 'tau', 4.0, @isnumeric);
    addParameter(p, 'vt_method', 0);
    addParameter(p, 'vt_params', struct());
    addParameter(p, 'scan_step', 5, @isnumeric);
    addParameter(p, 'tol_depth', 0.5, @isnumeric);
    addParameter(p, 'verbose', true, @islogical);
    
    parse(p, comp_ref, press_ref, temp_ref, h_ref, depth_range, Pc, Tc, acentric, BIP, M_gmol, varargin{:});
    opts = p.Results;
    
    grading_model = lower(opts.grading_model);
    dTdh = opts.dTdh;
    Cp_coeffs = opts.Cp_coeffs;
    H_ig_ref = opts.H_ig_ref;
    tau = opts.tau;
    vt_method = opts.vt_method;
    vt_params = opts.vt_params;
    scan_step = opts.scan_step;
    tol_depth = opts.tol_depth;
    verbose = opts.verbose;
    
    % Determine EOS type from vt_method
    eos_type = get_eos_type(vt_method, length(comp_ref));
    
    %% Validate inputs for non-isothermal models
    if ismember(grading_model, {'hasse', 'kempers', 'firoozabadi'})
        if dTdh == 0
            warning('Non-isothermal model selected but dTdh = 0. Results equivalent to isothermal.');
        end
        if ismember(grading_model, {'hasse', 'kempers'}) && (isempty(Cp_coeffs) || isempty(H_ig_ref))
            error('Cp_coeffs and H_ig_ref are required for %s model.', grading_model);
        end
    end
    
    %% Initialize
    comp_ref = comp_ref(:);
    n = length(comp_ref);
    Pc = Pc(:);
    Tc = Tc(:);
    acentric = acentric(:);
    M_gmol = M_gmol(:);
    
    GOC_depth = NaN;
    GOC_type = 'none';
    GOC_pressure = NaN;
    GOC_comp = comp_ref;
    
    tol_sat = 1e-10;
    maxiter_sat = 500;
    max_bisection = 30;
    max_secant = 30;
    tol_DeltaK = 1e-6;
    
    h_min = min(depth_range);
    h_max = max(depth_range);
    
    if h_ref >= h_max
        dh = -scan_step;
        h_start = h_ref;
        h_end = h_min;
    else
        dh = scan_step;
        h_start = h_ref;
        h_end = h_max;
    end
    
    depths = h_start:dh:h_end;
    n_depths = length(depths);
    
    %% Create grading function handle
    grading_func = create_grading_function(grading_model, dTdh, Cp_coeffs, H_ig_ref, tau, vt_method, vt_params);
    
    %% Print header
    if verbose
        fprintf('=========================================================================\n');
        fprintf('GOC DETECTION - %s Model (%s EOS)\n', upper(grading_model), eos_type);
        fprintf('=========================================================================\n');
        fprintf('Reference: h = %.1f m, P = %.2f bar, T = %.1f K\n', h_ref, press_ref/1e5, temp_ref);
        if dTdh ~= 0
            fprintf('Temperature gradient: dT/dh = %.4f K/m\n', dTdh);
        end
        fprintf('VT method: %d\n', vt_method);
        fprintf('Scan: %.1f to %.1f m (step = %.1f m)\n', h_start, h_end, abs(dh));
        fprintf('=========================================================================\n\n');
    end
    
    %% Storage arrays
    info.depths = zeros(n_depths, 1);
    info.P_res = zeros(n_depths, 1);
    info.T = zeros(n_depths, 1);
    info.P_bub = zeros(n_depths, 1);
    info.P_dew = zeros(n_depths, 1);
    info.P_sat = zeros(n_depths, 1);
    info.Delta_P = zeros(n_depths, 1);
    info.Delta_K = zeros(n_depths, 1);
    info.sat_type = cell(n_depths, 1);
    info.comp = zeros(n_depths, n);
    info.grading_model = grading_model;
    info.eos_type = eos_type;
    info.vt_method = vt_method;
    
    %% Print table header
    if verbose
        fprintf('%-8s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n', ...
                'Depth', 'T[K]', 'P_res', 'P_bub', 'P_dew', 'P_sat', 'Sat Type', 'ΔP');
        fprintf('%s\n', repmat('-', 1, 88));
    end
    
    %% Scan depths
    GOC_found = false;
    GOC_bracket = [];
    prev_sign = NaN;
    prev_sat_type = '';
    prev_DeltaK = NaN;
    DeltaK_decreasing_count = 0;
    
    for i = 1:n_depths
        h = depths(i);
        info.depths(i) = h;
        
        try
            [comp_h, P_h, T_h, P_bub, P_dew] = grading_func(h, h_ref, comp_ref, press_ref, temp_ref, ...
                Pc, Tc, acentric, BIP, M_gmol);
            comp_h = comp_h(:);
        catch ME
            if verbose
                fprintf('Error at h = %.1f m: %s\n', h, ME.message);
            end
            continue;
        end
        
        info.T(i) = T_h;
        info.comp(i, :) = comp_h';
        
        P_bub_valid = P_bub;
        P_dew_valid = P_dew;
        if isnan(P_bub) || P_bub <= 0
            P_bub_valid = 0;
        end
        if isnan(P_dew) || P_dew <= 0
            P_dew_valid = 0;
        end
        
        if P_bub_valid >= P_dew_valid
            sat_type = 'Bubble';
            P_sat = P_bub_valid;
        else
            sat_type = 'Dew';
            P_sat = P_dew_valid;
        end
        
        if P_h > 0 && P_sat > 0
            Delta_P = (P_h - P_sat) / P_h;
        else
            Delta_P = NaN;
        end
        
        [Delta_K, ~] = calculate_DeltaK(comp_h, P_h, T_h, Pc, Tc, acentric, BIP, eos_type);
        
        info.P_res(i) = P_h;
        info.P_bub(i) = P_bub_valid;
        info.P_dew(i) = P_dew_valid;
        info.P_sat(i) = P_sat;
        info.Delta_P(i) = Delta_P;
        info.Delta_K(i) = Delta_K;
        info.sat_type{i} = sat_type;
        
        if verbose
            fprintf('%-8.1f %-10.1f %-10.2f %-10.2f %-10.2f %-10.2f %-10s %-10.4f\n', ...
                    h, T_h, P_h/1e5, P_bub_valid/1e5, P_dew_valid/1e5, P_sat/1e5, sat_type, Delta_P);
        end
        
        % Detection logic
        curr_sign = sign(Delta_P);
        
        % Check for SATURATED GOC (sign change in ΔP)
        if i > 1 && ~isnan(prev_sign) && ~isnan(curr_sign)
            if prev_sign * curr_sign < 0
                GOC_found = true;
                GOC_type = 'saturated';
                GOC_bracket = [depths(i-1), depths(i)];
                if verbose
                    fprintf('\n*** SATURATED GOC detected between %.1f and %.1f m ***\n', ...
                            GOC_bracket(1), GOC_bracket(2));
                end
                break;
            end
        end
        
        % Check for UNDERSATURATED GOC (Pb ≈ Pd, ΔK → 0)
        if ~isnan(Delta_K) && Delta_K < tol_DeltaK
            GOC_found = true;
            GOC_type = 'undersaturated';
            if i > 1
                GOC_bracket = [depths(i-1), depths(i)];
            else
                GOC_bracket = [depths(i), depths(i)];
            end
            if verbose
                fprintf('\n*** UNDERSATURATED GOC detected near %.1f m (ΔK = %.2e) ***\n', ...
                        h, Delta_K);
            end
            break;
        end
        
        % Check for saturation type transition (Bubble ↔ Dew)
        if i > 1 && ~isempty(prev_sat_type) && ~strcmp(sat_type, prev_sat_type)
            if abs(P_bub_valid - P_dew_valid) / max(P_bub_valid, P_dew_valid) < 0.05
                GOC_found = true;
                GOC_type = 'undersaturated';
                GOC_bracket = [depths(i-1), depths(i)];
                if verbose
                    fprintf('\n*** UNDERSATURATED GOC detected (Bubble-Dew transition) at %.1f m ***\n', h);
                end
                break;
            end
        end
        
        prev_sign = curr_sign;
        prev_sat_type = sat_type;
        prev_DeltaK = Delta_K;
    end
    
    %% Refine GOC location
    if GOC_found && ~isempty(GOC_bracket)
        if verbose
            fprintf('\n=========================================================================\n');
            fprintf('GOC Type: %s\n', upper(GOC_type));
            fprintf('Bracket: [%.2f, %.2f] m\n', GOC_bracket(1), GOC_bracket(2));
            fprintf('=========================================================================\n');
        end
        
        if strcmp(GOC_type, 'saturated')
            [GOC_depth, GOC_pressure, GOC_comp, GOC_temp] = refine_saturated_GOC(...
                GOC_bracket, h_ref, comp_ref, press_ref, temp_ref, ...
                Pc, Tc, acentric, BIP, M_gmol, grading_func, ...
                tol_depth, max_bisection, verbose);
        else
            [GOC_depth, GOC_pressure, GOC_comp, GOC_temp] = refine_undersaturated_GOC(...
                GOC_bracket, h_ref, comp_ref, press_ref, temp_ref, ...
                Pc, Tc, acentric, BIP, M_gmol, grading_func, eos_type, ...
                tol_depth, tol_DeltaK, max_secant, verbose);
        end
        
        info.GOC_temp = GOC_temp;
    else
        if verbose
            fprintf('\n*** No GOC found in depth range [%.1f, %.1f] m ***\n', h_min, h_max);
        end
    end
    
    %% Final summary
    if verbose && GOC_found
        fprintf('\n=========================================================================\n');
        fprintf('GOC RESULT\n');
        fprintf('=========================================================================\n');
        fprintf('  Depth    : %.2f m\n', GOC_depth);
        fprintf('  Type     : %s\n', GOC_type);
        fprintf('  Pressure : %.2f bar\n', GOC_pressure/1e5);
        if exist('GOC_temp', 'var')
            fprintf('  Temp     : %.1f K (%.1f °C)\n', GOC_temp, GOC_temp - 273.15);
        end
        fprintf('=========================================================================\n');
    end
    
end


%% ========================================================================
% HELPER: Get EOS type from vt_method
% =========================================================================
function eos_type = get_eos_type(vt_method, n)
    if isnumeric(vt_method) && length(vt_method) == n
        eos_type = 'PR';
    elseif isnumeric(vt_method) && vt_method >= 0 && vt_method <= 6
        eos_type = 'PR';
    elseif isnumeric(vt_method) && vt_method >= 7 && vt_method <= 10
        eos_type = 'SRK';
    else
        eos_type = 'PR';
    end
end


%% ========================================================================
% HELPER: Create grading function handle
% =========================================================================
function grading_func = create_grading_function(model, dTdh, Cp_coeffs, H_ig_ref, tau, vt_method, vt_params)

    switch lower(model)
        case 'isothermal'
            grading_func = @(h, h_ref, comp_ref, P_ref, T_ref, Pc, Tc, acentric, BIP, M_gmol) ...
                call_isothermal(h, h_ref, comp_ref, P_ref, T_ref, Pc, Tc, acentric, BIP, M_gmol, ...
                vt_method, vt_params);
            
        case 'hasse'
            grading_func = @(h, h_ref, comp_ref, P_ref, T_ref, Pc, Tc, acentric, BIP, M_gmol) ...
                call_hasse(h, h_ref, comp_ref, P_ref, T_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, ...
                Cp_coeffs, H_ig_ref, vt_method, vt_params);
            
        case 'kempers'
            grading_func = @(h, h_ref, comp_ref, P_ref, T_ref, Pc, Tc, acentric, BIP, M_gmol) ...
                call_kempers(h, h_ref, comp_ref, P_ref, T_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, ...
                Cp_coeffs, H_ig_ref, vt_method, vt_params);
            
        case 'firoozabadi'
            grading_func = @(h, h_ref, comp_ref, P_ref, T_ref, Pc, Tc, acentric, BIP, M_gmol) ...
                call_firoozabadi(h, h_ref, comp_ref, P_ref, T_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, ...
                Cp_coeffs, H_ig_ref, tau, vt_method, vt_params);
            
        otherwise
            error('Unknown grading model: %s', model);
    end
end


%% ========================================================================
% WRAPPER: Isothermal grading (main.m)
% =========================================================================
function [comp_h, P_h, T_h, P_bub, P_dew] = call_isothermal(h, h_ref, comp_ref, P_ref, T_ref, ...
    Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params)

    [comp_h, P_h, P_bub, P_dew] = main(h, h_ref, comp_ref, P_ref, T_ref, ...
        Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params);
    T_h = T_ref;
end


%% ========================================================================
% WRAPPER: Hasse non-isothermal (main_hasse.m)
% =========================================================================
function [comp_h, P_h, T_h, P_bub, P_dew] = call_hasse(h, h_ref, comp_ref, P_ref, T_ref, dTdh, ...
    Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params)

    [comp_h, P_h, T_h, P_bub, P_dew] = main_hasse(h, h_ref, comp_ref, P_ref, T_ref, dTdh, ...
        Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
end


%% ========================================================================
% WRAPPER: Kempers (main_kempers.m)
% =========================================================================
function [comp_h, P_h, T_h, P_bub, P_dew] = call_kempers(h, h_ref, comp_ref, P_ref, T_ref, dTdh, ...
    Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params)

    [comp_h, P_h, T_h, P_bub, P_dew] = main_kempers(h, h_ref, comp_ref, P_ref, T_ref, dTdh, ...
        Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params);
end


%% ========================================================================
% WRAPPER: Firoozabadi (main_firoozabadi.m)
% =========================================================================
function [comp_h, P_h, T_h, P_bub, P_dew] = call_firoozabadi(h, h_ref, comp_ref, P_ref, T_ref, dTdh, ...
    Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, tau, vt_method, vt_params)

    [comp_h, P_h, T_h, P_bub, P_dew] = main_firoozabadi(h, h_ref, comp_ref, P_ref, T_ref, dTdh, ...
        Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, tau, vt_method, vt_params);
end


%% ========================================================================
% HELPER: Calculate ΔK = Σ(ln K_i)² for undersaturated GOC detection
% =========================================================================
function [Delta_K, K_values] = calculate_DeltaK(comp, P, T, Pc, Tc, acentric, BIP, eos_type)

    n = length(comp);
    K_values = ones(n, 1);
    
    try
        K_wilson = wilsoneq(P, T, Pc, Tc, acentric);
        K_values = K_wilson(:);
        
        ln_K = log(K_values);
        Delta_K = sum(ln_K.^2);
    catch
        Delta_K = NaN;
    end
end


%% ========================================================================
% HELPER: Refine SATURATED GOC using Bisection
% =========================================================================
function [GOC_depth, GOC_pressure, GOC_comp, GOC_temp] = refine_saturated_GOC(...
    bracket, h_ref, comp_ref, press_ref, temp_ref, ...
    Pc, Tc, acentric, BIP, M_gmol, grading_func, ...
    tol_depth, max_iter, verbose)

    if verbose
        fprintf('\nRefining SATURATED GOC using BISECTION...\n');
    end
    
    h_low = min(bracket);
    h_high = max(bracket);
    
    [comp_low, P_low, T_low, Pb_low, Pd_low] = grading_func(h_low, h_ref, comp_ref, press_ref, temp_ref, ...
        Pc, Tc, acentric, BIP, M_gmol);
    Psat_low = max(Pb_low, Pd_low);
    f_low = P_low - Psat_low;
    
    [comp_high, P_high, T_high, Pb_high, Pd_high] = grading_func(h_high, h_ref, comp_ref, press_ref, temp_ref, ...
        Pc, Tc, acentric, BIP, M_gmol);
    Psat_high = max(Pb_high, Pd_high);
    f_high = P_high - Psat_high;
    
    if verbose
        fprintf('  h_low=%.2f m (f=%.2e), h_high=%.2f m (f=%.2e)\n', h_low, f_low, h_high, f_high);
    end
    
    for iter = 1:max_iter
        h_mid = (h_low + h_high) / 2;
        
        [comp_mid, P_mid, T_mid, Pb_mid, Pd_mid] = grading_func(h_mid, h_ref, comp_ref, press_ref, temp_ref, ...
            Pc, Tc, acentric, BIP, M_gmol);
        Psat_mid = max(Pb_mid, Pd_mid);
        f_mid = P_mid - Psat_mid;
        
        if verbose
            fprintf('  Iter %2d: h = %.3f m, P = %.2f bar, Psat = %.2f bar, f = %.2e\n', ...
                    iter, h_mid, P_mid/1e5, Psat_mid/1e5, f_mid);
        end
        
        if abs(h_high - h_low) < tol_depth || abs(f_mid) < 1e-3 * P_mid
            GOC_depth = h_mid;
            GOC_pressure = P_mid;
            GOC_comp = comp_mid;
            GOC_temp = T_mid;
            if verbose
                fprintf('  Converged! GOC at h = %.2f m\n', GOC_depth);
            end
            return;
        end
        
        if f_low * f_mid < 0
            h_high = h_mid;
            f_high = f_mid;
        else
            h_low = h_mid;
            f_low = f_mid;
        end
    end
    
    GOC_depth = h_mid;
    GOC_pressure = P_mid;
    GOC_comp = comp_mid;
    GOC_temp = T_mid;
    if verbose
        fprintf('  Max iterations. Best estimate: h = %.2f m\n', GOC_depth);
    end
end


%% ========================================================================
% HELPER: Refine UNDERSATURATED GOC using Secant Method
% =========================================================================
function [GOC_depth, GOC_pressure, GOC_comp, GOC_temp] = refine_undersaturated_GOC(...
    bracket, h_ref, comp_ref, press_ref, temp_ref, ...
    Pc, Tc, acentric, BIP, M_gmol, grading_func, eos_type, ...
    tol_depth, tol_DeltaK, max_iter, verbose)

    if verbose
        fprintf('\nRefining UNDERSATURATED GOC using SECANT method...\n');
    end
    
    h_0 = bracket(1);
    h_1 = bracket(2);
    
    [comp_0, P_0, T_0, Pb_0, Pd_0] = grading_func(h_0, h_ref, comp_ref, press_ref, temp_ref, ...
        Pc, Tc, acentric, BIP, M_gmol);
    [Delta_K_0, ~] = calculate_DeltaK(comp_0, P_0, T_0, Pc, Tc, acentric, BIP, eos_type);
    f_0 = Delta_K_0;
    
    [comp_1, P_1, T_1, Pb_1, Pd_1] = grading_func(h_1, h_ref, comp_ref, press_ref, temp_ref, ...
        Pc, Tc, acentric, BIP, M_gmol);
    [Delta_K_1, ~] = calculate_DeltaK(comp_1, P_1, T_1, Pc, Tc, acentric, BIP, eos_type);
    f_1 = Delta_K_1;
    
    GOC_depth = h_1;
    GOC_pressure = P_1;
    GOC_comp = comp_1;
    GOC_temp = T_1;
    
    for iter = 1:max_iter
        if abs(f_1 - f_0) < 1e-15
            break;
        end
        
        h_new = h_1 - f_1 * (h_1 - h_0) / (f_1 - f_0);
        
        h_new = max(min(bracket), min(max(bracket), h_new));
        
        [comp_new, P_new, T_new, Pb_new, Pd_new] = grading_func(h_new, h_ref, comp_ref, press_ref, temp_ref, ...
            Pc, Tc, acentric, BIP, M_gmol);
        [Delta_K_new, ~] = calculate_DeltaK(comp_new, P_new, T_new, Pc, Tc, acentric, BIP, eos_type);
        f_new = Delta_K_new;
        
        if verbose
            fprintf('  Iter %2d: h = %.3f m, ΔK = %.2e\n', iter, h_new, Delta_K_new);
        end
        
        if abs(h_new - h_1) < tol_depth || Delta_K_new < tol_DeltaK
            GOC_depth = h_new;
            GOC_pressure = P_new;
            GOC_comp = comp_new;
            GOC_temp = T_new;
            if verbose
                fprintf('  Converged! Undersaturated GOC at h = %.2f m\n', GOC_depth);
            end
            return;
        end
        
        h_0 = h_1;
        f_0 = f_1;
        h_1 = h_new;
        f_1 = f_new;
        
        GOC_depth = h_new;
        GOC_pressure = P_new;
        GOC_comp = comp_new;
        GOC_temp = T_new;
    end
    
    if verbose
        fprintf('  Max iterations. Best estimate: h = %.2f m\n', GOC_depth);
    end
end