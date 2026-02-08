function [rho, V_m, Z, c_mix] = calculate_density(comp, press, temp, Pc, Tc, acentric, BIP, M_gmol, vt_method, vt_params)
% CALCULATE_DENSITY - Calculate mixture density from EOS with volume translation
%
% INPUTS:
%   comp      : Composition [mole fractions]
%   press     : Pressure [Pa]
%   temp      : Temperature [K]
%   Pc        : Critical pressures [Pa]
%   Tc        : Critical temperatures [K]
%   acentric  : Acentric factors [-]
%   BIP       : Binary interaction parameter matrix
%   M_gmol    : Molecular weights [g/mol]
%   vt_method : Volume translation method (same as main_hasse):
%               0-6: PR-based, 7-10: SRK-based
%   vt_params : Structure with Vc, Zc, components (as needed)
%
% OUTPUTS:
%   rho       : Density [kg/m³]
%   V_m       : Molar volume [m³/mol]
%   Z         : Compressibility factor (from EOS, before VT)
%   c_mix     : Mixture volume translation [m³/mol]

    R = 8.3144598;
    
    if nargin < 10
        vt_params = struct();
    end
    if nargin < 9
        vt_method = 0;
    end
    
    comp = comp(:);
    comp = comp / sum(comp);
    n = length(comp);
    
    [vt_type, eos_type, vt_opts] = parse_vt_method(vt_method, vt_params, n, temp, Pc, Tc, acentric, M_gmol);
    
    [~, Z] = fugacitycoef_multicomp(comp, press, temp, Pc, Tc, acentric, BIP, eos_type);
    
    V_m_eos = Z * R * temp / press;
    
    if vt_type ~= 0 && vt_type ~= 7
        c_i = calc_vt(comp, press, temp, Pc, Tc, acentric, M_gmol, vt_type, eos_type, vt_opts);
        c_mix = sum(comp .* c_i);
    else
        c_mix = 0;
    end
    
    V_m = V_m_eos - c_mix;
    
    M_mix = sum(comp .* M_gmol) / 1000;
    
    rho = M_mix / V_m;

end


function c = calc_vt(comp, press, temp, Pc, Tc, acentric, M_gmol, vt_type, eos_type, vt_opts)

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
            vt_opts.c_direct = vt_params.c_custom(:) * 1e-6;
        else
            error('For vt_method = 6, provide vt_params.c_custom [cm³/mol]');
        end
    end
    
end


function num = string_to_method_number(name)
    name = lower(char(name));
    switch name
        case {'pr', '0'}, num = 0;
        case {'peneloux', '1'}, num = 1;
        case {'magoulas', 'mt', '2'}, num = 2;
        case {'ungerer', 'ub', '3'}, num = 3;
        case {'baled_pr', '4'}, num = 4;
        case {'abudour', '5'}, num = 5;
        case {'custom', '6'}, num = 6;
        case {'srk', '7'}, num = 7;
        case {'pina', 'pm', '8'}, num = 8;
        case {'chen_li', 'cl', '9'}, num = 9;
        case {'baled_srk', '10'}, num = 10;
        otherwise, num = 0;
    end
end