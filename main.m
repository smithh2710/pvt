function [comp_h, press_h, pressbub_h, pressdew_h] = main(h, h_ref, comp_ref, press_ref, temp, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params)
% MAIN - Isothermal Compositional Grading (Schulte Model)
% =========================================================================
% Solves the isothermal compositional grading equations (Schulte 1980):
%
%   ln(f_i^h / f_i^ref) = M_i * g * (h - h_ref) / (R*T)
%   Σ z_i = 1
%
% This is the isothermal limit of non-isothermal grading (dT/dh = 0).
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
%   h         : Target depth [m]
%   h_ref     : Reference depth [m]
%   comp_ref  : Reference composition [mole fractions]
%   press_ref : Reference pressure [Pa]
%   temp      : Temperature [K] (constant for isothermal)
%   Pc        : Critical pressures [Pa]
%   Tc        : Critical temperatures [K]
%   acentric  : Acentric factors [-]
%   BIP       : Binary interaction parameter matrix
%   M_gmol    : Molecular weights [g/mol]
%   vt_method : Volume translation method (0-10) or direct c vector [cm³/mol]
%   vt_params : Structure with additional parameters:
%               .Vc         - Critical volumes [cm³/mol]
%               .components - Cell array of component names
%               .c_custom   - Custom c values for method 6 [cm³/mol]
%               .Zc         - Critical compressibility factors (for Chen-Li)
%
% OUTPUTS:
%   comp_h    : Composition at depth h [mole fractions]
%   press_h   : Pressure at depth h [Pa]
%   pressbub_h: Bubble point pressure at depth h [Pa]
%   pressdew_h: Dew point pressure at depth h [Pa]
%
% =========================================================================
% REFERENCES:
%   Schulte, A.M. (1980). SPE 9235
%   Whitson, C.H. & Brule, M.R. (2000). Phase Behavior, SPE Monograph
% =========================================================================

    R = 8.3144598;
    g = 9.80665;
    M = M_gmol / 1000;

    tol = 1e-10;
    maxiter = 1500;

    n = length(comp_ref);
    
    comp_ref = comp_ref(:);
    comp_ref = comp_ref / sum(comp_ref);
    M = M(:);
    M_gmol = M_gmol(:);
    Pc = Pc(:);
    Tc = Tc(:);
    acentric = acentric(:);

    if nargin < 12
        vt_params = struct();
    end
    if nargin < 11
        vt_method = 0;
    end

    [vt_type, eos_type, vt_opts] = parse_vt_method(vt_method, vt_params, n, temp, Pc, Tc, acentric, M_gmol);
    
    c_ref = calc_vt(comp_ref, press_ref, temp, Pc, Tc, acentric, M_gmol, vt_type, eos_type, vt_opts);

    [fugcoef_ref_eos, ~] = fugacitycoef_multicomp(comp_ref, press_ref, temp, Pc, Tc, acentric, BIP, eos_type);
    
    if vt_type ~= 0 && vt_type ~= 7
        fugcoef_ref = fugcoef_ref_eos .* exp(-c_ref * press_ref / (R * temp));
    else
        fugcoef_ref = fugcoef_ref_eos;
    end
    
    f_ref = fugcoef_ref .* comp_ref * press_ref;

    f_h = f_ref .* exp((M * g * (h - h_ref)) / (R * temp));

    initial_guess = [comp_ref; press_ref];

    fun = @(x) residual_fugacity(x(1:n), x(end), f_h, temp, Pc, Tc, acentric, BIP, M_gmol, vt_type, eos_type, vt_opts, R);
    
    options = optimoptions('fsolve', 'Display', 'none', 'FunctionTolerance', 1e-15, 'StepTolerance', 1e-15);
    solution = fsolve(fun, initial_guess, options);

    comp_h = solution(1:n);
    press_h = solution(end);

    comp_h = max(comp_h, 1e-15);
    comp_h = comp_h / sum(comp_h);

   try 
        [pressbub_h, ~] = pressbub_multicomp_newton(comp_h, press_h, temp, Pc, Tc, acentric, BIP, tol, maxiter, eos_type);
    catch 
       pressbub_h = NaN; 
    end
    if isnan(pressbub_h) || ~isreal(pressbub_h) || pressbub_h <= 0 || pressbub_h > press_h
        try
            [pressbub_h, ~] = pressbub_multicomp_ss(comp_h, press_h, temp, Pc, Tc, acentric, BIP, tol, maxiter);
        catch
            pressbub_h = NaN;
        end
    end
    if ~isnan(pressbub_h) && (pressbub_h <= 0 || pressbub_h > press_h)
        pressbub_h = NaN;
    end

    try 
        [pressdew_h, ~] = pressdew_multicomp_newton(comp_h, press_h, temp, Pc, Tc, acentric, BIP, tol, maxiter, eos_type);
    catch 
       pressdew_h = NaN; 
    end
    if isnan(pressdew_h) || ~isreal(pressdew_h) || pressdew_h <= 0 || pressdew_h > press_h
        try
            [pressdew_h, ~] = pressdew_multicomp_ss(comp_h, press_h, temp, Pc, Tc, acentric, BIP, tol, maxiter);
        catch
            pressdew_h = NaN;
        end
    end
    if ~isnan(pressdew_h) && (pressdew_h <= 0 || pressdew_h > press_h)
        pressdew_h = NaN;
    end

end


function F = residual_fugacity(comp_h, press_h, f_h, temp, Pc, Tc, acentric, BIP, M_gmol, vt_type, eos_type, vt_opts, R)

    n = length(f_h);
    
    comp_h = max(comp_h, 1e-15);
    comp_h = comp_h / sum(comp_h);
    press_h = max(press_h, 1e2);
    
    [fugcoef_h_eos, ~] = fugacitycoef_multicomp(comp_h, press_h, temp, Pc, Tc, acentric, BIP, eos_type);
    
    if vt_type ~= 0 && vt_type ~= 7
        c_shift = calc_vt(comp_h, press_h, temp, Pc, Tc, acentric, M_gmol, vt_type, eos_type, vt_opts);
        fugcoef_h = fugcoef_h_eos .* exp(-c_shift * press_h / (R * temp));
    else
        fugcoef_h = fugcoef_h_eos;
    end
    
    F = zeros(n+1, 1);
    for i = 1:n
        F(i) = comp_h(i) * press_h * fugcoef_h(i) - f_h(i);
    end
    F(n+1) = sum(comp_h) - 1;
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
            % c = -c ; 
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