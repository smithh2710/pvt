function [GOR, Bo, Rs, info] = calculate_GOR_STO(comp, P_res, T_res, Pc, Tc, acentric, BIP, M_gmol, varargin)
% CALCULATE_GOR_STO - Calculate Gas-Oil Ratio via flash to stock tank conditions
% =========================================================================
% Performs a flash calculation from reservoir conditions to stock tank
% conditions (STO: 1 atm, 15°C) to determine:
%   - GOR: Gas-Oil Ratio [Sm³/Sm³]
%   - Bo:  Oil Formation Volume Factor [res m³/Sm³ STO]
%   - Rs:  Solution Gas-Oil Ratio [Sm³/Sm³]
%
% Method:
%   1. Flash reservoir fluid to stock tank conditions (1 atm, 288.15 K)
%   2. Calculate molar volumes of gas and liquid phases
%   3. Compute GOR = V_gas / V_oil at standard conditions
%
% =========================================================================
% INPUTS (Required):
%   comp     : Reservoir fluid composition [mole fractions]
%   P_res    : Reservoir pressure [Pa]
%   T_res    : Reservoir temperature [K]
%   Pc       : Critical pressures [Pa]
%   Tc       : Critical temperatures [K]
%   acentric : Acentric factors [-]
%   BIP      : Binary interaction parameters
%   M_gmol   : Molecular weights [g/mol]
%
% INPUTS (Optional Name-Value Pairs):
%   'vt_method'  : Volume translation method (0-6 or direct c vector, default: 0)
%   'vt_params'  : VT parameters structure
%   'eos_type'   : 'PR' (default) or 'SRK'
%   'T_STO'      : Stock tank temperature [K] (default: 288.15 = 15°C)
%   'P_STO'      : Stock tank pressure [Pa] (default: 101325 = 1 atm)
%   'verbose'    : Print output (default: false)
%
% OUTPUTS:
%   GOR  : Gas-Oil Ratio [Sm³/Sm³] (Inf for dry gas)
%   Bo   : Oil Formation Volume Factor [res m³/Sm³ STO]
%   Rs   : Solution Gas-Oil Ratio [Sm³/Sm³] (same as GOR for single-stage)
%   info : Structure with detailed calculation results
%
% =========================================================================
% USAGE EXAMPLES:
%
%   % Basic usage (PR-EOS, no VT)
%   [GOR, Bo, Rs] = calculate_GOR_STO(comp, P_res, T_res, Pc, Tc, w, BIP, M);
%
%   % With SRK-EOS
%   [GOR, Bo, Rs] = calculate_GOR_STO(comp, P_res, T_res, Pc, Tc, w, BIP, M, ...
%       'eos_type', 'SRK');
%
%   % With Abudour VT
%   vt_params.Vc = Vc;
%   vt_params.components = components;
%   [GOR, Bo, Rs] = calculate_GOR_STO(comp, P_res, T_res, Pc, Tc, w, BIP, M, ...
%       'vt_method', 5, 'vt_params', vt_params);
%
%   % With JY volume translation and verbose output
%   [GOR, Bo, Rs, info] = calculate_GOR_STO(comp, P_res, T_res, Pc, Tc, w, BIP, M, ...
%       'vt_method', c_JY, 'verbose', true);
%
%   % Custom STO conditions (60°F, 14.7 psia)
%   [GOR, Bo, Rs] = calculate_GOR_STO(comp, P_res, T_res, Pc, Tc, w, BIP, M, ...
%       'T_STO', 288.71, 'P_STO', 101325);
%
% =========================================================================
% REFERENCES:
%   Whitson, C.H. & Brule, M.R. (2000). Phase Behavior, SPE Monograph
%   Pedersen, K.S. & Christensen, P.L. (2007). Phase Behavior of Petroleum
% =========================================================================

    %% Parse inputs
    p = inputParser;
    p.FunctionName = 'calculate_GOR_STO';
    
    addRequired(p, 'comp', @isnumeric);
    addRequired(p, 'P_res', @isnumeric);
    addRequired(p, 'T_res', @isnumeric);
    addRequired(p, 'Pc', @isnumeric);
    addRequired(p, 'Tc', @isnumeric);
    addRequired(p, 'acentric', @isnumeric);
    addRequired(p, 'BIP', @isnumeric);
    addRequired(p, 'M_gmol', @isnumeric);
    
    addParameter(p, 'vt_method', 0);
    addParameter(p, 'vt_params', struct());
    addParameter(p, 'eos_type', 'PR', @ischar);
    addParameter(p, 'T_STO', 288.15, @isnumeric);
    addParameter(p, 'P_STO', 101325, @isnumeric);
    addParameter(p, 'verbose', false, @islogical);
    
    parse(p, comp, P_res, T_res, Pc, Tc, acentric, BIP, M_gmol, varargin{:});
    opts = p.Results;
    
    vt_method = opts.vt_method;
    vt_params = opts.vt_params;
    eos_type = upper(opts.eos_type);
    T_STO = opts.T_STO;
    P_STO = opts.P_STO;
    verbose = opts.verbose;
    
    %% Validate EOS type
    if ~ismember(eos_type, {'PR', 'SRK'})
        error('eos_type must be ''PR'' or ''SRK''. Got: %s', eos_type);
    end
    
    %% Initialize
    R = 8.3144598;
    
    comp = comp(:);
    n = length(comp);
    Pc = Pc(:);
    Tc = Tc(:);
    acentric = acentric(:);
    M_gmol = M_gmol(:);
    
    comp = comp / sum(comp);
    
    GOR = NaN;
    Bo = NaN;
    Rs = NaN;
    
    info = struct();
    info.eos_type = eos_type;
    info.P_STO = P_STO;
    info.T_STO = T_STO;
    info.P_res = P_res;
    info.T_res = T_res;
    
    if verbose
        fprintf('=========================================================================\n');
        fprintf('GOR CALCULATION - Flash to Stock Tank (%s EOS)\n', eos_type);
        fprintf('=========================================================================\n');
        fprintf('Reservoir: P = %.2f bar, T = %.1f K (%.1f °C)\n', P_res/1e5, T_res, T_res-273.15);
        fprintf('Stock Tank: P = %.4f bar, T = %.1f K (%.1f °C)\n', P_STO/1e5, T_STO, T_STO-273.15);
        fprintf('=========================================================================\n\n');
    end
    
    %% Step 1: Flash to Stock Tank Conditions
    try
        [K_sto, comp_vap_sto, comp_liq_sto, V_frac] = vaporliquideq(P_STO, T_STO, ...
            comp, Pc, Tc, acentric, BIP, 1e-10, 200, eos_type);
        
        flash_converged = true;
        
        if V_frac < 1e-10
            V_frac = 0;
            comp_liq_sto = comp;
            comp_vap_sto = comp;
        elseif V_frac > (1 - 1e-10)
            V_frac = 1;
            comp_vap_sto = comp;
            comp_liq_sto = comp;
        end
        
    catch ME
        flash_converged = false;
        if verbose
            fprintf('Flash failed: %s\n', ME.message);
            fprintf('Using Wilson K-values as fallback...\n');
        end
        
        K_wilson = wilsoneq(P_STO, T_STO, Pc, Tc, acentric);
        V_frac = solve_RR(comp, K_wilson);
        
        if V_frac < 0
            V_frac = 0;
        elseif V_frac > 1
            V_frac = 1;
        end
        
        comp_vap_sto = K_wilson .* comp ./ (1 + V_frac * (K_wilson - 1));
        comp_liq_sto = comp ./ (1 + V_frac * (K_wilson - 1));
    end
    
    comp_vap_sto = comp_vap_sto(:);
    comp_liq_sto = comp_liq_sto(:);
    comp_vap_sto = comp_vap_sto / sum(comp_vap_sto);
    comp_liq_sto = comp_liq_sto / sum(comp_liq_sto);
    
    info.V_frac = V_frac;
    info.comp_vap_sto = comp_vap_sto;
    info.comp_liq_sto = comp_liq_sto;
    info.flash_converged = flash_converged;
    
    if verbose
        fprintf('Flash Results:\n');
        fprintf('  Vapor fraction: %.4f\n', V_frac);
        fprintf('  Liquid fraction: %.4f\n', 1 - V_frac);
        if V_frac > 0 && V_frac < 1
            fprintf('  Vapor C1: %.2f mol%%\n', comp_vap_sto(3)*100);
            fprintf('  Liquid C1: %.2f mol%%\n', comp_liq_sto(3)*100);
        end
        fprintf('\n');
    end
    
    %% Step 2: Calculate Molar Volumes at STO
    
    V_gas_molar = R * T_STO / P_STO;
    
    try
        [~, Z_liq] = fugacitycoef_multicomp_liquid(comp_liq_sto, P_STO, T_STO, ...
            Pc, Tc, acentric, BIP, eos_type);
        V_liq_molar_EOS = Z_liq * R * T_STO / P_STO;
        
        c_VT = calc_volume_translation(comp_liq_sto, P_STO, T_STO, Pc, Tc, acentric, ...
            vt_method, vt_params, n, eos_type, M_gmol);
        c_mix = sum(comp_liq_sto .* c_VT);
        V_liq_molar = V_liq_molar_EOS - c_mix;
        
        if V_liq_molar <= 0
            V_liq_molar = V_liq_molar_EOS;
        end
        
    catch ME
        if verbose
            fprintf('EOS liquid volume failed: %s\n', ME.message);
            fprintf('Using density estimate...\n');
        end
        
        MW_liq = sum(comp_liq_sto .* M_gmol);
        rho_liq = 700;
        V_liq_molar = MW_liq / (rho_liq * 1000);
    end
    
    info.V_gas_molar = V_gas_molar;
    info.V_liq_molar = V_liq_molar;
    
    if verbose
        fprintf('Molar Volumes at STO:\n');
        fprintf('  Gas: %.4f m³/mol (ideal gas)\n', V_gas_molar);
        fprintf('  Liquid: %.6f m³/mol (EOS + VT)\n', V_liq_molar);
        fprintf('\n');
    end
    
    %% Step 3: Calculate GOR
    
    n_gas = V_frac;
    n_liq = 1 - V_frac;
    
    V_gas_std = n_gas * V_gas_molar;
    V_liq_std = n_liq * V_liq_molar;
    
    if V_liq_std > 1e-15
        GOR = V_gas_std / V_liq_std;
    else
        GOR = Inf;
    end
    
    info.n_gas = n_gas;
    info.n_liq = n_liq;
    info.V_gas_std = V_gas_std;
    info.V_liq_std = V_liq_std;
    
    %% Step 4: Calculate Oil Formation Volume Factor (Bo)
    
    try
        [~, Z_res] = fugacitycoef_multicomp_liquid(comp, P_res, T_res, ...
            Pc, Tc, acentric, BIP, eos_type);
        V_res_molar_EOS = Z_res * R * T_res / P_res;
        
        c_VT_res = calc_volume_translation(comp, P_res, T_res, Pc, Tc, acentric, ...
            vt_method, vt_params, n, eos_type, M_gmol);
        c_mix_res = sum(comp .* c_VT_res);
        V_res_molar = V_res_molar_EOS - c_mix_res;
        
        if V_res_molar <= 0
            V_res_molar = V_res_molar_EOS;
        end
        
    catch ME
        if verbose
            fprintf('EOS reservoir volume failed: %s\n', ME.message);
        end
        
        MW_res = sum(comp .* M_gmol);
        rho_res = 600;
        V_res_molar = MW_res / (rho_res * 1000);
    end
    
    if V_liq_std > 1e-15
        Bo = V_res_molar / V_liq_molar;
    else
        Bo = NaN;
    end
    
    info.V_res_molar = V_res_molar;
    
    Rs = GOR;
    
    %% Final output
    if verbose
        fprintf('=========================================================================\n');
        fprintf('RESULTS\n');
        fprintf('=========================================================================\n');
        fprintf('  GOR : %.2f Sm³/Sm³\n', GOR);
        fprintf('  Bo  : %.4f res m³/Sm³ STO\n', Bo);
        fprintf('  Rs  : %.2f Sm³/Sm³\n', Rs);
        fprintf('=========================================================================\n');
    end
    
end


%% =========================================================================
% Helper: Solve Rachford-Rice equation
% =========================================================================
function V = solve_RR(z, K)

    z = z(:);
    K = K(:);
    
    K_max = max(K);
    K_min = min(K);
    
    if K_max < 1
        V = 0;
        return;
    end
    if K_min > 1
        V = 1;
        return;
    end
    
    V_min = 1 / (1 - K_max);
    V_max = 1 / (1 - K_min);
    
    V_min = max(0, V_min + 1e-10);
    V_max = min(1, V_max - 1e-10);
    
    if V_min >= V_max
        V = 0.5;
        return;
    end
    
    tol = 1e-12;
    for iter = 1:100
        V = (V_min + V_max) / 2;
        f = sum(z .* (K - 1) ./ (1 + V * (K - 1)));
        
        if abs(f) < tol
            break;
        end
        
        if f > 0
            V_min = V;
        else
            V_max = V;
        end
    end
end


%% =========================================================================
% Helper: Calculate volume translation
% =========================================================================
function c = calc_volume_translation(comp, press, temp, Pc, Tc, acentric, vt_method, vt_params, n, eos_type, M_gmol)

    c = zeros(n, 1);
    
    if isnumeric(vt_method) && length(vt_method) > 1
        c = vt_method(:) * 1e-6;
        if length(c) ~= n
            error('Volume translation vector length (%d) must match number of components (%d)', length(c), n);
        end
        return;
    end
    
    method = vt_method;
    if method == 0
        return;
    end
    
    switch method
        case 1  % Peneloux
            result = peneloux_volume_shift(Pc, Tc, acentric, comp, vt_params.Vc, eos_type);
            c = result.c_i * 1e-6;
            
        case 2  % Magoulas-Tassios (returns m³/mol)
            [~, c, ~] = magoulas_tassios_volume_shift(temp, Pc, Tc, acentric, comp, vt_params.Vc, vt_params.components);
            
        case 3  % Ungerer-Batut (returns m³/mol)
            [~, c, ~] = ungerer_batut_volume_shift(temp, Pc, Tc, acentric, M_gmol, comp, vt_params.Vc);
            
        case 4  % Baled (returns m³/mol)
            [~, c, ~] = baled_volume_shift(temp, Pc, Tc, acentric, M_gmol, vt_params.components, comp, vt_params.Vc, eos_type);
            
        case 5  % Abudour (returns m³/mol)
            if ~isfield(vt_params, 'Vc') || ~isfield(vt_params, 'components')
                error('Abudour method requires vt_params.Vc and vt_params.components');
            end
            Vc = vt_params.Vc;
            components = vt_params.components;
            BIP_local = zeros(n);
            [~, c] = abudour_volume_shift(comp, press, temp, Pc, Tc, acentric, Vc, M_gmol, components, BIP_local, eos_type);
            
        case 6  % Custom (input in cm³/mol)
            if ~isfield(vt_params, 'c_custom')
                error('Method 6 requires vt_params.c_custom');
            end
            c = vt_params.c_custom(:) * 1e-6;
            
        otherwise
            warning('Unknown VT method %d. Using no volume translation.', method);
    end
end