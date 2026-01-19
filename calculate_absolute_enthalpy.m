function [H_abs_components, H_abs_mixture, H_ig_specific, H_res_specific, H_abs_specific] = ...
    calculate_absolute_enthalpy(temp, press, comp, Pc, Tc, acentric, BIP, Mw, Cp_coeffs, H_ig_ref_input)

R = 8.314462618;
T_ref = 273.15;

n_comp = length(comp);
comp = comp(:) / sum(comp);
Pc = Pc(:);
Tc = Tc(:);
acentric = acentric(:);
Mw = Mw(:);

H_ig_ref = zeros(n_comp, 1);

if isempty(H_ig_ref_input)
    H_ig_over_MR = -1342./Mw + 8.367;
    H_ig_ref = H_ig_over_MR .* Mw * R;
elseif isnumeric(H_ig_ref_input) && length(H_ig_ref_input) == n_comp
    H_ig_over_MR = H_ig_ref_input(:);
    H_ig_ref = H_ig_over_MR .* Mw * R;
else
    error('H_ig_ref_input must be empty or a vector of length n_comp');
end

H_ig_components = zeros(n_comp, 1);
for i = 1:n_comp
    C1 = Cp_coeffs(i, 1);
    C2 = Cp_coeffs(i, 2);
    C3 = Cp_coeffs(i, 3);
    C4 = Cp_coeffs(i, 4);
    
    integral_Cp = C1 * (temp - T_ref) + ...
                  C2/2 * (temp^2 - T_ref^2) + ...
                  C3/3 * (temp^3 - T_ref^3) + ...
                  C4/4 * (temp^4 - T_ref^4);
    
    H_ig_components(i) = H_ig_ref(i) + integral_Cp;
end

H_res_components = calculate_partial_molar_residual_enthalpy(temp, press, comp, Pc, Tc, acentric, BIP, R);

H_abs_components = H_ig_components + H_res_components;

H_abs_mixture = sum(comp .* H_abs_components);

H_ig_specific = H_ig_components ./ Mw;
H_res_specific = H_res_components ./ Mw;
H_abs_specific = H_abs_components ./ Mw;

end


function H_res = calculate_partial_molar_residual_enthalpy(temp, press, comp, Pc, Tc, acentric, BIP, R)

n_comp = length(comp);
H_res = zeros(n_comp, 1);

dT = min(0.5, 0.001 * temp);

try
    [phi_minus, ~] = fugacitycoef_multicomp(comp, press, temp - dT, Pc, Tc, acentric, BIP);
    [phi_plus, ~] = fugacitycoef_multicomp(comp, press, temp + dT, Pc, Tc, acentric, BIP);
    
    for i = 1:n_comp
        if phi_plus(i) > 0 && phi_minus(i) > 0
            dln_phi_dT = (log(phi_plus(i)) - log(phi_minus(i))) / (2 * dT);
            H_res(i) = -R * temp^2 * dln_phi_dT;
        else
            H_res(i) = 0;
        end
    end
catch
    H_res = zeros(n_comp, 1);
end

for i = 1:n_comp
    if ~isfinite(H_res(i))
        H_res(i) = 0;
    end
end

end