function [H_abs_components, H_abs_mixture, H_ig_specific, H_res_specific, H_abs_specific] = calculate_absolute_enthalpy(temp, press, comp, Pc, Tc, acentric, BIP, Mw, Cp_coeffs, H_ig_ref_option)
% Calculate absolute enthalpy for each component and mixture
% Based on Pedersen & Hjermstad (2006) SPE 101275
%
% Inputs:
% temp - Temperature in Kelvin
% press - Pressure in Pa
% comp - Mole fraction vector for all components (will be normalized)
% Pc - Critical pressure vector (Pa)
% Tc - Critical temperature vector (K)
% acentric - Acentric factor vector
% BIP - Binary interaction parameters matrix
% Mw - Molecular weights vector (g/mol)
% Cp_coeffs - Heat capacity coefficients matrix (n_comp x 4) [C1, C2, C3, C4]
%             where Cp = C1 + C2*T + C3*T^2 + C4*T^3 (J/mol/K)
% H_ig_ref_option - Either:
%                   1) Vector of H_ig(273.15K) in J/mol for each component
%                   2) Empty [] to use Pedersen correlation (Eq. 8)
%                   3) Vector of H_ig/(M*R) values in K/g units
%
% Outputs:
% H_abs_components - Absolute enthalpy for each component (J/mol)
% H_abs_mixture - Mixture absolute enthalpy (J/mol)
% H_ig_specific - Ideal gas specific enthalpy H_ig/M (J/MOL)
% H_res_specific - Residual specific enthalpy H_res/M (J/MOL)
% H_abs_specific - Absolute specific enthalpy H_abs/M (J/MOL)  ( J/g )

% Constants
R = 8.314462618; % J/mol/K
T_ref = 273.15; % Reference temperature (K)

% Get dimensions and normalize composition
n_comp = length(comp);
comp = comp(:) / sum(comp);
Pc = Pc(:);
Tc = Tc(:);
acentric = acentric(:);
Mw = Mw(:);

% Step 1: Determine reference ideal gas enthalpy at 273.15 K
if isempty(H_ig_ref_option)
    % Use Pedersen correlation (Equation 8 from paper)
    H_ig_ref = R * (-1342 + 8.367 * Mw); % J/mol

elseif length(H_ig_ref_option) == n_comp
    % if max(abs(H_ig_ref_option)) < 1000
    %     % Input is likely H_ig/(M*R) in K/g units
    %     % H_ig_ref = H_ig_ref_option .* Mw * R; % Convert to J/mol
    % else
    %     % Input is already in J/mol
        H_ig_ref = H_ig_ref_option;
    % end
else
    error('H_ig_ref_option must be empty or have same length as composition vector');
end

% Step 2: Calculate ideal gas enthalpy at temperature T
H_ig_components = zeros(n_comp, 1);
for i = 1:n_comp
    % Integration of Cp from T_ref to T
    delta_H = Cp_coeffs(i, 1) * (temp - T_ref) + ...
              Cp_coeffs(i, 2)/2 * (temp^2 - T_ref^2) + ...
              Cp_coeffs(i, 3)/3 * (temp^3 - T_ref^3) + ...
              Cp_coeffs(i, 4)/4 * (temp^4 - T_ref^4);
    
    H_ig_components(i) = H_ig_ref(i) + delta_H;
end

% Step 3: Calculate partial molar residual enthalpies
H_res_components = calculate_residual_enthalpy_components(temp, press, comp, Pc, Tc, acentric, BIP);

% Step 4: Calculate absolute enthalpies
H_abs_components = H_ig_components + H_res_components;

% Step 5: Calculate mixture properties
H_abs_mixture = sum(comp .* H_abs_components);

% Step 6: Calculate specific enthalpies (per unit mass)
H_ig_specific = H_ig_components ./ Mw;  % J/g
H_res_specific = H_res_components ./ Mw; % J/g
H_abs_specific = H_abs_components ./ Mw; % J/g

end

%% Residual enthalpy calculation function
function H_res_components = calculate_residual_enthalpy_components(temp, press, comp, Pc, Tc, acentric, BIP)
% Calculate partial molar residual enthalpies using numerical differentiation
% Based on H̃_i^res = -R*T^2 * (∂ln φ_i/∂T)_P,x

R = 8.314462618; % J/mol/K
n_comp = length(comp);
H_res_components = zeros(n_comp, 1);

% Adaptive step size for numerical differentiation
dT = min(0.5, 0.001 * temp); % Use 0.5 K or 0.1% of temperature

% Calculate fugacity coefficients at perturbed temperatures
try
    [phi_minus, ~] = fugacitycoef_multicomp(comp, press, temp - dT, Pc, Tc, acentric, BIP);
    [phi_plus, ~] = fugacitycoef_multicomp(comp, press, temp + dT, Pc, Tc, acentric, BIP);
    
    % Use central difference for better accuracy
    for i = 1:n_comp
        if phi_plus(i) > 0 && phi_minus(i) > 0
            dln_phi_dT = (log(phi_plus(i)) - log(phi_minus(i))) / (2 * dT);
            H_res_components(i) = -R * temp^2 * dln_phi_dT;
        else
            % Fallback to forward difference if needed
            [phi_0, ~] = fugacitycoef_multicomp(comp, press, temp, Pc, Tc, acentric, BIP);
            [phi_forward, ~] = fugacitycoef_multicomp(comp, press, temp + dT, Pc, Tc, acentric, BIP);
            if phi_0(i) > 0 && phi_forward(i) > 0
                dln_phi_dT = (log(phi_forward(i)) - log(phi_0(i))) / dT;
                H_res_components(i) = -R * temp^2 * dln_phi_dT;
            else
                H_res_components(i) = 0; % Default if calculation fails
            end
        end
    end
catch
    % If fugacity calculation fails, try with smaller step
    dT = dT / 10;
    try
        [phi_minus, ~] = fugacitycoef_multicomp(comp, press, temp - dT, Pc, Tc, acentric, BIP);
        [phi_plus, ~] = fugacitycoef_multicomp(comp, press, temp + dT, Pc, Tc, acentric, BIP);
        
        for i = 1:n_comp
            if phi_plus(i) > 0 && phi_minus(i) > 0
                dln_phi_dT = (log(phi_plus(i)) - log(phi_minus(i))) / (2 * dT);
                H_res_components(i) = -R * temp^2 * dln_phi_dT;
            else
                H_res_components(i) = 0;
            end
        end
    catch
        % Final fallback - return zeros
        H_res_components = zeros(n_comp, 1);
    end
end

% Ensure finite values
H_res_components(~isfinite(H_res_components)) = 0;

end
