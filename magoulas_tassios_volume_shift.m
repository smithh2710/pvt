function [s, c, mix] = magoulas_tassios_volume_shift(T, Pc, Tc, acentric, comp, Vc, components)
%MAGOULAS_TASSIOS_VOLUME_SHIFT Exponential temperature-dependent VT for PR EOS
%
% Reference:
%   Magoulas, K., & Tassios, D. (1990)
%   "Thermophysical properties of n-alkanes from C1 to C20 and their 
%    prediction for higher ones"
%   Fluid Phase Equilibria, 56, 119-140
%
% VALID RANGE: n-alkanes C1 to C20
%
% Volume translation (Eq. 17):
%   c = -[c_m0 + (δ_c - c_m0) * exp(β * |1 - T/Tc|)]
%
% where:
%   c_m0 = (R*Tc/Pc) * (k_m0 + k_m1*ω + k_m2*ω² + k_m3*ω³ + k_m4*ω⁴)  (18)
%   β = l_m0 + l_m1*ω²                                                  (19)
%   δ_c = (R*Tc/Pc) * (Zc_PR - Zc_exp)                                 (20)
%   Zc_exp = 0.289 - 0.0701*ω - 0.0207*ω²                              (21)
%
% Parameters:
%   k_m0 to k_m4: -0.014471, 0.067498, -0.084852, 0.067298, -0.017366
%   l_m0 = -10.2447, l_m1 = -28.6312
%   Zc_PR = 0.3074
%
% Physical behavior:
%   - At Tr → 1 (critical): c → -δ_c (negative, adds volume)
%   - At Tr → 0 (low T):    c → -c_m0 (positive, reduces volume)
%   - Near-critical: gives appropriate NEGATIVE c for light HC
%
% Modified alpha function (Eq. 22):
%   m = d_m0 + d_m1*ω + d_m2*ω² + d_m3*ω³ + d_m4*ω⁴
%   d_m0 to d_m4: 0.384401, 1.52276, -0.213808, 0.034616, -0.001976
%
% Inputs:
%   T        : Temperature [K]
%   Pc       : Critical pressures [Pa]
%   Tc       : Critical temperatures [K]
%   acentric : Acentric factors [-]
%   comp     : Mole fractions [-]
%   Vc       : (Optional) Critical volumes [cm³/mol]
%   components : (Optional) Cell array of component names
%
% Outputs:
%   s   : Dimensionless volume shift [-]
%   c   : Volume shift [m³/mol]
%   mix : Mixture properties structure

R = 8.3144598;
Omega_b = 0.07780;
Zc_PR = 0.3074;

k_m = [-0.014471, 0.067498, -0.084852, 0.067298, -0.017366];
l_m0 = -10.2447;
l_m1 = -28.6312;

d_m = [0.384401, 1.52276, -0.213808, 0.034616, -0.001976];

Pc = Pc(:);
Tc = Tc(:);
acentric = acentric(:);
comp = comp(:);
comp = comp / sum(comp);

ncomp = length(Pc);

if nargin < 6 || isempty(Vc)
    Zc_est = 0.29056 - 0.08775 * acentric;
    Vc = Zc_est .* R .* Tc ./ Pc * 1e6;
else
    Vc = Vc(:);
end

if nargin < 7 || isempty(components)
    components = cell(ncomp, 1);
elseif ischar(components)
    components = {components};
end

Tr = T ./ Tc;

b_i = Omega_b * R * Tc ./ Pc;

Zc_exp = 0.289 - 0.0701 * acentric - 0.0207 * acentric.^2;

delta_c = R * Tc ./ Pc .* (Zc_PR - Zc_exp);

c_m0_factor = k_m(1) + k_m(2)*acentric + k_m(3)*acentric.^2 + ...
              k_m(4)*acentric.^3 + k_m(5)*acentric.^4;
c_m0 = R * Tc ./ Pc .* c_m0_factor;

beta = l_m0 + l_m1 * acentric.^2;

exp_term = exp(beta .* abs(1 - Tr));
c = -1 * (c_m0 + (delta_c - c_m0) .* exp_term);

s = c ./ b_i;

m_modified = d_m(1) + d_m(2)*acentric + d_m(3)*acentric.^2 + ...
             d_m(4)*acentric.^3 + d_m(5)*acentric.^4;

Vc_23 = Vc.^(2/3);
seta = (comp .* Vc_23) / sum(comp .* Vc_23);

c_mix = sum(seta .* c);
s_mix = sum(seta .* s);

mix = struct();
mix.Zc_exp = Zc_exp;
mix.Zc_PR = Zc_PR;
mix.delta_c = delta_c;
mix.delta_c_cm3 = delta_c * 1e6;
mix.c_m0 = c_m0;
mix.c_m0_cm3 = c_m0 * 1e6;
mix.beta = beta;
mix.exp_term = exp_term;
mix.s = s;
mix.c = c;
mix.c_cm3 = c * 1e6;
mix.m_modified = m_modified;
mix.Tr = Tr;
mix.seta = seta;
mix.c_mix = c_mix;
mix.c_mix_cm3 = c_mix * 1e6;
mix.s_mix = s_mix;

end