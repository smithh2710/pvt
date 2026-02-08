function [H_mix_specific, H_partial_specific, H_ig_MR, H_res_MR, H_total_MR] = calculate_absolute_enthalpy(T, P, comp, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref, eos_type)
% CALCULATE_ABSOLUTE_ENTHALPY - Partial molar enthalpies for Haase model
%
% Implements Pedersen et al. (2015) equations:
%   H_i_partial = H_i^ig + H_i^res_partial                    (Eq. 7)
%   H_mix = sum(z_i * H_i_partial)                            (Eq. 8)
%   H_i^res_partial = -R*T^2 * (d ln(phi_i)/dT)_{P,x}         (from Eq. 9)
%   H_i^ig(T) = H_i^ig(T_ref) + integral(Cp_i dT, T_ref, T)   (Eq. 10)
%
% INPUTS:
%   T          : Temperature [K]
%   P          : Pressure [Pa]
%   comp       : Mole fractions [-] (column vector)
%   Pc         : Critical pressures [Pa]
%   Tc         : Critical temperatures [K]
%   acentric   : Acentric factors [-]
%   BIP        : Binary interaction parameter matrix
%   M_gmol     : Molecular weights [g/mol]
%   Cp_coeffs  : Ideal gas Cp coefficients [n x 4], Cp = C1 + C2*T + C3*T^2 + C4*T^3 [J/(mol·K)]
%   H_ig_ref   : Ideal gas enthalpy at T_ref=273.15 K [J/mol]
%   eos_type   : 'PR' (default) or 'SRK'
%
% OUTPUTS:
%   H_mix_specific     : Mixture specific enthalpy H_mix/M_mix [J/g]
%   H_partial_specific : Partial molar specific enthalpy H_i/M_i [J/g] for each component
%   H_ig_MR            : Ideal gas specific enthalpy H_i^ig/(M_i*R) [K] for each component
%   H_res_MR           : Residual specific enthalpy H_i^res/(M_i*R) [K] for each component
%   H_total_MR         : Total specific enthalpy H_i/(M_i*R) [K] for each component

    R = 8.3144598;    % Universal gas constant [J/(mol·K)]
    T_ref = 273.15;     % Reference temperature [K]

    % Handle optional input
    if nargin < 11 || isempty(eos_type)
        eos_type = 'PR';
    end

    n = length(comp);
    comp = comp(:);
    comp = comp / sum(comp);

    Pc = Pc(:);
    Tc = Tc(:);
    acentric = acentric(:);
    M_gmol = M_gmol(:);
    H_ig_ref = H_ig_ref(:);

    %% Ideal gas enthalpy: H_i^ig(T) = H_i^ig(T_ref) + integral(Cp dT)
    H_ig = zeros(n, 1);
    for i = 1:n
        C1 = Cp_coeffs(i, 1);
        C2 = Cp_coeffs(i, 2);
        C3 = Cp_coeffs(i, 3);
        C4 = Cp_coeffs(i, 4);
        
        % Analytical integration of Cp from T_ref to T
        delta_H_ig = C1 * (T - T_ref) + ...
                     C2 / 2 * (T^2 - T_ref^2) + ...
                     C3 / 3 * (T^3 - T_ref^3) + ...
                     C4 / 4 * (T^4 - T_ref^4);
        
        H_ig(i) = H_ig_ref(i) + delta_H_ig;
    end

    %% Partial molar residual enthalpy: H_i^res = -R*T^2 * (d ln phi_i / dT)
    H_res = calculate_partial_molar_residual_enthalpy(T, P, comp, Pc, Tc, acentric, BIP, R, eos_type);

    %% Total partial molar enthalpy [J/mol]
    H_partial = H_ig + H_res;
    
    %% Mixture enthalpy [J/mol] (Eq. 8)
    H_mix = sum(comp .* H_partial);
    
    %% Mixture molecular weight [g/mol]
    M_mix = sum(comp .* M_gmol);
    
    %% Output 1: H_mix/M_mix [J/g]
    H_mix_specific = H_mix / M_mix;
    
    %% Output 2: H_i/M_i [J/g]
    H_partial_specific = H_partial ./ M_gmol;
    
    %% Output 3: H_i^ig/(M_i*R) [K]
    H_ig_MR = H_ig ./ (M_gmol * R);
    
    %% Output 4: H_i^res/(M_i*R) [K]
    H_res_MR = H_res ./ (M_gmol * R);
    
    %% Output 5: H_i/(M_i*R) [K]
    H_total_MR = H_partial ./ (M_gmol * R);

end


function H_res = calculate_partial_molar_residual_enthalpy(T, P, comp, Pc, Tc, acentric, BIP, R, eos_type)
% Partial molar residual enthalpy via numerical differentiation
%   H_i^res = -R * T^2 * (d ln(phi_i) / dT)_{P,x}

    n = length(comp);
    H_res = zeros(n, 1);

    % Adaptive step size: larger of 0.1 K or 0.01% of T
    dT = max(0.1, 1e-4 * T);

    % Central difference
    [phi_minus, ~] = fugacitycoef_multicomp(comp, P, T - dT, Pc, Tc, acentric, BIP, eos_type);
    [phi_plus, ~]  = fugacitycoef_multicomp(comp, P, T + dT, Pc, Tc, acentric, BIP, eos_type);

    for i = 1:n
        if phi_plus(i) > 1e-15 && phi_minus(i) > 1e-15
            dln_phi_dT = (log(phi_plus(i)) - log(phi_minus(i))) / (2 * dT);
            H_res(i) = -R * T^2 * dln_phi_dT;
        else
            % Fallback to forward difference
            [phi_0, ~] = fugacitycoef_multicomp(comp, P, T, Pc, Tc, acentric, BIP, eos_type);
            [phi_fwd, ~] = fugacitycoef_multicomp(comp, P, T + dT, Pc, Tc, acentric, BIP, eos_type);
            if phi_0(i) > 1e-15 && phi_fwd(i) > 1e-15
                dln_phi_dT = (log(phi_fwd(i)) - log(phi_0(i))) / dT;
                H_res(i) = -R * T^2 * dln_phi_dT;
            else
                warning('Component %d has invalid fugacity coefficient', i);
                H_res(i) = 0;
            end
        end
    end

    % Check for numerical issues
    bad_idx = ~isfinite(H_res);
    if any(bad_idx)
        warning('Non-finite residual enthalpies for %d components', sum(bad_idx));
        H_res(bad_idx) = 0;
    end
end