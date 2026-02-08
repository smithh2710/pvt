function [comp_h, press_h, temp_h, pressbub_h, pressdew_h] = main_kempers(h_target, h_ref, comp_ref, P_ref, T_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_method, vt_params)
% MAIN_KEMPERS - Non-isothermal Compositional Grading (Kempers Model)
% =========================================================================
% Solves the non-isothermal grading equation using Kempers (1989) NHT:
%
%   Net Heat of Transport (Eq. 3.60 from Nikpoor 2014):
%     Q_i^net = (v_i / V_m) * H_m - H_i
%
%   Main equation:
%     ln(f_i^h / f_i^ref) = M_i*g*Δh/(R*T) - Q_i^net*ΔT/(R*T²)
%
%   Where:
%     v_i = partial molar volume [m³/mol]
%     V_m = mixture molar volume [m³/mol]
%     H_i = partial molar enthalpy [J/mol]
%     H_m = mixture molar enthalpy [J/mol]
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
%   Direct input:
%     Vector of c values [cm³/mol] with length = n_components (auto-detected)
%
% =========================================================================
% INPUTS:
%   h_target  : Target depth [m]
%   h_ref     : Reference depth [m]
%   comp_ref  : Reference composition [mole fractions]
%   P_ref     : Reference pressure [Pa]
%   T_ref     : Reference temperature [K]
%   dTdh      : Geothermal gradient [K/m] (positive = T increases with depth)
%   Pc        : Critical pressures [Pa]
%   Tc        : Critical temperatures [K]
%   acentric  : Acentric factors [-]
%   BIP       : Binary interaction parameter matrix
%   M_gmol    : Molecular weights [g/mol]
%   Cp_coeffs : Ideal gas Cp coefficients [n x 4], Cp = A + B*T + C*T² + D*T³
%   H_ig_ref  : Reference ideal gas enthalpies at 273.15 K [J/mol]
%   vt_method : Volume translation method (0-10) or direct c vector [cm³/mol]
%   vt_params : Structure with additional parameters:
%               .Vc         - Critical volumes [cm³/mol]
%               .components - Cell array of component names
%               .c_custom   - Custom c values for method 6 [cm³/mol]
%               .Zc         - Critical compressibility factors (for Chen-Li)
%
% OUTPUTS:
%   comp_h     : Composition at depth h [mole fractions]
%   press_h    : Pressure at depth h [Pa]
%   temp_h     : Temperature at depth h [K]
%   pressbub_h : Bubble point pressure at depth h [Pa]
%   pressdew_h : Dew point pressure at depth h [Pa]
%
% =========================================================================
% REFERENCES:
%   Kempers, L.J.T.M. (1989). J. Chem. Phys. 90(11), 6541-6548
%   Nikpoor, M.H. (2014). MSc Thesis, University of Calgary
% =========================================================================

    if nargin < 15
        vt_params = struct();
    end
    if nargin < 14
        vt_method = 0;
    end

    n = length(comp_ref);
    step_size = 1;
    
    comp_ref = comp_ref(:);
    comp_ref = comp_ref / sum(comp_ref);
    Pc = Pc(:);
    Tc = Tc(:);
    acentric = acentric(:);
    M_gmol = M_gmol(:);
    H_ig_ref = H_ig_ref(:);

    [vt_type, eos_type, vt_opts] = parse_vt_method(vt_method, vt_params, n, T_ref, Pc, Tc, acentric, M_gmol);

    tol = 1e-10;
    maxiter = 1500;

    delta_h_total = h_target - h_ref;
    n_steps = ceil(abs(delta_h_total) / step_size);
    
    if n_steps == 0
        n_steps = 1;
    end
    
    dh = delta_h_total / n_steps;

    h_current = h_ref;
    P_current = P_ref;
    T_current = T_ref;
    comp_current = comp_ref;

    for step = 1:n_steps
        h_next = h_current + dh;
        
        [comp_next, P_next, T_next] = kempers_single_step(...
            h_next, h_current, comp_current, P_current, T_current, dTdh, ...
            Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_type, eos_type, vt_opts);
        
        h_current = h_next;
        P_current = P_next;
        T_current = T_next;
        comp_current = comp_next;
    end

    comp_h = comp_current;
    press_h = P_current;
    temp_h = T_current;

    comp_h = max(comp_h, 1e-15);
    comp_h = comp_h / sum(comp_h);

    try
        [pressbub_h, ~] = pressbub_multicomp_newton(comp_h, press_h, temp_h, Pc, Tc, acentric, BIP, tol, maxiter, eos_type);
    catch
        pressbub_h = NaN;
    end
    if isnan(pressbub_h) || ~isreal(pressbub_h) || pressbub_h <= 0 || pressbub_h > press_h
        try
            [pressbub_h, ~] = pressbub_multicomp_ss(comp_h, press_h, temp_h, Pc, Tc, acentric, BIP, tol, maxiter);
        catch
            pressbub_h = NaN;
        end
    end
    if ~isnan(pressbub_h) && (pressbub_h <= 0 || pressbub_h > press_h)
        pressbub_h = NaN;
    end


    try
        [pressdew_h, ~] = pressdew_multicomp_newton(comp_h, press_h, temp_h, Pc, Tc, acentric, BIP, tol, maxiter, eos_type);
    catch
        pressdew_h = NaN;
    end
    if isnan(pressdew_h) || ~isreal(pressdew_h) || pressdew_h <= 0 || pressdew_h > press_h
        try
            [pressdew_h, ~] = pressdew_multicomp_ss(comp_h, press_h, temp_h, Pc, Tc, acentric, BIP, tol, maxiter);
        catch
            pressdew_h = NaN;
        end
    end
    if ~isnan(pressdew_h) && (pressdew_h <= 0 || pressdew_h > press_h)
        pressdew_h = NaN;
    end

end


function [comp_h, press_h, temp_h] = kempers_single_step(h, h_ref, comp_ref, P, T, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_type, eos_type, vt_opts)

    R = 8.3144598;
    g = 9.80665;
    n = length(comp_ref);

    comp_ref = comp_ref(:);
    comp_ref = comp_ref / sum(comp_ref);

    delta_h = h - h_ref;
    delta_T = dTdh * delta_h;
    temp_h = T + delta_T;

    c_ref = calc_vt(comp_ref, P, T, Pc, Tc, acentric, M_gmol, vt_type, eos_type, vt_opts);

    [fugcoef_ref_eos, ~] = fugacitycoef_multicomp(comp_ref, P, T, Pc, Tc, acentric, BIP, eos_type);
    
    if vt_type ~= 0 && vt_type ~= 7
        fugcoef_ref = fugcoef_ref_eos .* exp(-c_ref * P / (R * T));
    else
        fugcoef_ref = fugcoef_ref_eos;
    end
    
    f_ref = fugcoef_ref .* comp_ref * P;

    [H_partial, v_partial] = calc_partial_molar_properties(temp_h, P, comp_ref, ...
        Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_type, eos_type, vt_opts);

    H_mix_molar = sum(comp_ref .* H_partial);
    V_mix_molar = sum(comp_ref .* v_partial);

    term2 = zeros(n, 1);
    term3 = zeros(n, 1);

    for i = 1:n
        M_i_kg = M_gmol(i) / 1000;
        
        term2(i) = M_i_kg * g * delta_h / (R * T);
        
        Q_i_net = (v_partial(i) / V_mix_molar) * H_mix_molar - H_partial(i);
        term3(i) = Q_i_net * delta_T / (R * T^2);
    end

    params = struct();
    params.f_ref = f_ref;
    params.term2 = term2;
    params.term3 = term3;
    params.T_ref = T;
    params.Pc = Pc;
    params.Tc = Tc;
    params.acentric = acentric;
    params.BIP = BIP;
    params.M_gmol = M_gmol;
    params.n = n;
    params.vt_type = vt_type;
    params.eos_type = eos_type;
    params.vt_opts = vt_opts;
    params.R = R;

    initial_guess = [comp_ref; P];

    options = optimoptions('fsolve', ...
        'Display', 'none', ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12, ...
        'MaxIterations', 200);

    [solution, ~, ~] = fsolve(@(x) residual_kempers(x, params), initial_guess, options);

    comp_h = solution(1:n);
    press_h = solution(end);

    comp_h = max(comp_h, 1e-15);
    comp_h = comp_h / sum(comp_h);

end


function F = residual_kempers(x, params)

    f_ref = params.f_ref;
    term2 = params.term2;
    term3 = params.term3;
    T_ref = params.T_ref;
    Pc = params.Pc;
    Tc = params.Tc;
    acentric = params.acentric;
    BIP = params.BIP;
    M_gmol = params.M_gmol;
    n = params.n;
    vt_type = params.vt_type;
    eos_type = params.eos_type;
    vt_opts = params.vt_opts;
    R = params.R;

    z_h = x(1:n);
    P_h = x(end);

    z_h = max(z_h, 1e-15);
    z_h_norm = z_h / sum(z_h);
    P_h = max(P_h, 1e2);

    [fugcoef_h_eos, ~] = fugacitycoef_multicomp(z_h_norm, P_h, T_ref, Pc, Tc, acentric, BIP, eos_type);

    if vt_type ~= 0 && vt_type ~= 7
        c_h = calc_vt(z_h_norm, P_h, T_ref, Pc, Tc, acentric, M_gmol, vt_type, eos_type, vt_opts);
        fugcoef_h = fugcoef_h_eos .* exp(-c_h * P_h / (R * T_ref));
    else
        fugcoef_h = fugcoef_h_eos;
    end

    F = zeros(n+1, 1);

    for i = 1:n
        f_i_h = fugcoef_h(i) * z_h_norm(i) * P_h;
        F(i) = log(f_i_h / f_ref(i)) - term2(i) + term3(i);
    end

    F(n+1) = sum(z_h) - 1;

end


function [H_partial, v_partial] = calc_partial_molar_properties(T, P, comp, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, vt_type, eos_type, vt_opts)
% Partial molar enthalpy and volume
%
% Enthalpy: H_i = H_i^ig + H_i^res
%   H_i^res = -R*T² * d(ln φ_i)/dT
%
% Volume: v_i = RT/P + RT * d(ln φ_i)/dP
%   v_i^VT = v_i^EOS - c_i

    R = 8.3144598;
    T_ref_ig = 273.15;
    n = length(comp);
    comp = comp(:);
    
    % Residual enthalpy via numerical differentiation
    dT = max(0.1, 1e-4 * T);
    
    [phi_T_minus, ~] = fugacitycoef_multicomp(comp, P, T - dT, Pc, Tc, acentric, BIP, eos_type);
    [phi_T_plus, ~] = fugacitycoef_multicomp(comp, P, T + dT, Pc, Tc, acentric, BIP, eos_type);

    H_residual = zeros(n, 1);
    for i = 1:n
        if phi_T_plus(i) > 0 && phi_T_minus(i) > 0
            dln_phi_dT = (log(phi_T_plus(i)) - log(phi_T_minus(i))) / (2 * dT);
            H_residual(i) = -R * T^2 * dln_phi_dT;
        else
            H_residual(i) = 0;
        end
    end

    % Ideal gas enthalpy
    H_ig = zeros(n, 1);
    
    for i = 1:n
        A = Cp_coeffs(i, 1);
        B = Cp_coeffs(i, 2);
        C = Cp_coeffs(i, 3);
        D = Cp_coeffs(i, 4);
        if size(Cp_coeffs, 2) >= 5
            E = Cp_coeffs(i, 5);
        else
            E = 0;
        end
        
        delta_H_ig = A * (T - T_ref_ig) + ...
                     B/2 * (T^2 - T_ref_ig^2) + ...
                     C/3 * (T^3 - T_ref_ig^3) + ...
                     D/4 * (T^4 - T_ref_ig^4) + ...
                     E/5 * (T^5 - T_ref_ig^5);
        
        H_ig(i) = H_ig_ref(i) + delta_H_ig;
    end

    H_partial = H_ig + H_residual;

    % Partial molar volume via numerical differentiation of fugacity
    dP = max(1000, 1e-5 * P);
    [phi_P_minus, ~] = fugacitycoef_multicomp(comp, P - dP, T, Pc, Tc, acentric, BIP, eos_type);
    [phi_P_plus, ~] = fugacitycoef_multicomp(comp, P + dP, T, Pc, Tc, acentric, BIP, eos_type);

    v_partial_eos = zeros(n, 1);
    for i = 1:n
        if phi_P_plus(i) > 0 && phi_P_minus(i) > 0
            dln_phi_dP = (log(phi_P_plus(i)) - log(phi_P_minus(i))) / (2 * dP);
            v_partial_eos(i) = R * T / P + R * T * dln_phi_dP;
        else
            v_partial_eos(i) = R * T / P;
        end
    end

    % Apply volume translation
    if vt_type ~= 0 && vt_type ~= 7
        c = calc_vt(comp, P, T, Pc, Tc, acentric, M_gmol, vt_type, eos_type, vt_opts);
        v_partial = v_partial_eos - c;
    else
        v_partial = v_partial_eos;
    end

end


function c = calc_vt(comp, press, temp, Pc, Tc, acentric, M_gmol, vt_type, eos_type, vt_opts)
% Volume translation [m³/mol]
%
% Methods 0-6: PR-based
% Methods 7-10: SRK-based

    n = length(comp);
    comp = comp(:);
    
    switch vt_type
        
        case 0
            c = zeros(n, 1);
            
        case 1
            result = peneloux_volume_shift(Pc, Tc, acentric, comp, vt_opts.Vc, 'PR');
            c = result.c_i * 1e-6;
            
        case 2
            [~, c, ~] = magoulas_tassios_volume_shift(temp, Pc, Tc, acentric, comp, vt_opts.Vc, vt_opts.components);
            
        case 3
            [~, c, ~] = ungerer_batut_volume_shift(temp, Pc, Tc, acentric, M_gmol, comp, vt_opts.Vc);
            
        case 4
            [~, c, ~] = baled_volume_shift(temp, Pc, Tc, acentric, M_gmol, vt_opts.components, comp, vt_opts.Vc, 'PR');
            
        case 5
            BIP_local = zeros(n);
            [~, c, ~] = abudour_volume_shift(comp, press, temp, Pc, Tc, acentric, vt_opts.Vc, M_gmol, vt_opts.components, BIP_local, 'PR');
            
        case 6
            c = vt_opts.c_direct;
            
        case 7
            c = zeros(n, 1);
            
        case 8
            result = pina_martinez_volume_shift(temp, Pc, Tc, acentric, vt_opts.components, comp);
            c = result.c_i * 1e-6;
            
        case 9
            [~, c, ~] = chen_li_volume_shift(comp, press, temp, Pc, Tc, acentric, vt_opts.Zc, vt_opts.components, false);
            
        case 10
            [~, c, ~] = baled_volume_shift(temp, Pc, Tc, acentric, M_gmol, vt_opts.components, comp, vt_opts.Vc, 'SRK');
            
        case -1
            c = vt_opts.c_direct;
            
        otherwise
            c = zeros(n, 1);
    end
    
    c = c(:);
end


function [vt_type, eos_type, vt_opts] = parse_vt_method(vt_method, vt_params, n, temp, Pc, Tc, acentric, M_gmol)
% Parse VT method and determine EOS type
%
% Methods 0-6: PR-EOS
% Methods 7-10: SRK-EOS

    R = 8.3144598;
    vt_opts = struct();
    
    if isfield(vt_params, 'Vc') && ~isempty(vt_params.Vc)
        vt_opts.Vc = vt_params.Vc(:);
    else
        Zc_est = 0.2905 - 0.085 * acentric;
        vt_opts.Vc = Zc_est .* R .* Tc ./ Pc * 1e6;
    end
    
    if isfield(vt_params, 'Zc') && ~isempty(vt_params.Zc)
        vt_opts.Zc = vt_params.Zc(:);
    else
        vt_opts.Zc = 0.2905 - 0.085 * acentric(:);
    end
    
    if isfield(vt_params, 'components') && ~isempty(vt_params.components)
        vt_opts.components = vt_params.components;
    else
        vt_opts.components = cell(n, 1);
    end
    
    vt_opts.Pc = Pc;
    vt_opts.Tc = Tc;
    vt_opts.acentric = acentric;
    vt_opts.M_gmol = M_gmol;
    vt_opts.temp = temp;
    vt_opts.n = n;
    
    if isempty(vt_method)
        vt_type = 0;
        eos_type = 'PR';
        return;
    end
    
    if isnumeric(vt_method) && length(vt_method) == n
        vt_type = -1;
        eos_type = 'PR';
        vt_opts.c_direct = vt_method(:) * 1e-6;
        return;
    end
    
    if ischar(vt_method) || isstring(vt_method)
        vt_type = string_to_method_number(vt_method);
    else
        vt_type = vt_method;
    end
    
    if vt_type >= 0 && vt_type <= 6
        eos_type = 'PR';
    elseif vt_type >= 7 && vt_type <= 10
        eos_type = 'SRK';
    else
        vt_type = 0;
        eos_type = 'PR';
    end
    
    if vt_type == 6
        if isfield(vt_params, 'c_custom') && ~isempty(vt_params.c_custom)
            c_custom = vt_params.c_custom(:);
            if length(c_custom) ~= n
                error('c_custom must have %d elements (one per component)', n);
            end
            vt_opts.c_direct = c_custom * 1e-6;
        else
            error('For vt_method = 6, you must provide vt_params.c_custom [cm³/mol]');
        end
    end
    
end


function num = string_to_method_number(name)

    name = lower(char(name));
    
    switch name
        case {'pr', 'pr_novt', '0'}
            num = 0;
        case {'pr_peneloux', 'peneloux', 'pen', '1'}
            num = 1;
        case {'pr_mt', 'magoulas_tassios', 'magoulas', 'mt', '2'}
            num = 2;
        case {'pr_ub', 'ungerer_batut', 'ungerer', 'ub', '3'}
            num = 3;
        case {'pr_baled', 'baled_pr', '4'}
            num = 4;
        case {'pr_abudour', 'abudour', 'abu', '5'}
            num = 5;
        case {'pr_custom', 'custom', '6'}
            num = 6;
        case {'srk', 'srk_novt', '7'}
            num = 7;
        case {'srk_pm', 'pina_martinez', 'pina', 'pm', '8'}
            num = 8;
        case {'srk_cl', 'chen_li', 'chen', 'cl', '9'}
            num = 9;
        case {'srk_baled', 'baled_srk', '10'}
            num = 10;
        otherwise
            num = 0;
    end
end