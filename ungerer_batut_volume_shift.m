function [s, c, mix] = ungerer_batut_volume_shift(T, Pc, Tc, acentric, MW, comp, Vc)
%UNGERER_BATUT_VOLUME_SHIFT Temperature-dependent volume translation for PR EOS
%
% Reference:
%   Ungerer, P., & Batut, C. (1997)
%   "Prédiction des propriétés volumétriques des hydrocarbures par une 
%    translation de volume améliorée"
%   Oil & Gas Science and Technology - Rev. IFP, 52(6), 609-623
%
% VALID RANGE: C6-C40 hydrocarbons (MW ≈ 86-563 g/mol)
%              Temperature: 273-573 K
%              Pressure: 5-120 MPa (optimized at 50 MPa)
%
% Volume translation (Eq. 11):
%   c(T) = (p1 + p2·MW)·T + (p3 + p4·MW)
%
% Parameters:
%   p1 = 0.023 cm³/(mol·K)
%   p2 = -0.00056 cm³/(g·K)
%   p3 = -34.5 cm³/mol
%   p4 = 0.4666 cm³/g
%
% Simplified form:
%   c(T) = (0.023 - 0.00056·MW)·T + (-34.5 + 0.4666·MW)
%
% Alternative (Eq. 9) for n-alkanes using carbon number:
%   c(T) = (0.021 - 0.008n)·T + (-64 + 9.6n)
%
% Mixing rule (Eq. 13):
%   c_mix = Σ xi·ci (linear mole-fraction weighted)
%
% For components OUTSIDE valid range (MW < 86), Peneloux fallback is used.
%
% Inputs:
%   T        : Temperature [K]
%   Pc       : Critical pressures [Pa]
%   Tc       : Critical temperatures [K]
%   acentric : Acentric factors [-]
%   MW       : Molecular weights [g/mol]
%   comp     : Mole fractions [-]
%   Vc       : (Optional) Critical volumes [cm³/mol]
%
% Outputs:
%   s   : Dimensionless volume shift [-]
%   c   : Volume shift [m³/mol]
%   mix : Mixture properties structure

R = 8.3144598;
Omega_b = 0.07780;

p1 = 0.023;
p2 = -0.00056;
p3 = -34.5;
p4 = 0.4666;

MW_min = 86;

Pc = Pc(:);
Tc = Tc(:);
acentric = acentric(:);
MW = MW(:);
comp = comp(:);
comp = comp / sum(comp);

ncomp = length(Pc);

if nargin < 7 || isempty(Vc)
    Zc_est = 0.29056 - 0.08775 * acentric;
    Vc = Zc_est .* R .* Tc ./ Pc * 1e6;
else
    Vc = Vc(:);
end

b_i = Omega_b * R * Tc ./ Pc;

Zra = 0.29056 - 0.08775 * acentric;

A_coef = zeros(ncomp, 1);
B_coef = zeros(ncomp, 1);
c = zeros(ncomp, 1);
c_cm3 = zeros(ncomp, 1);
valid = true(ncomp, 1);
source = cell(ncomp, 1);

for i = 1:ncomp
    if MW(i) < MW_min
        valid(i) = false;
        c_pen = 0.40768 * (0.29441 - Zra(i)) * R * Tc(i) / Pc(i);
        c(i) = c_pen;
        c_cm3(i) = c_pen * 1e6;
        A_coef(i) = NaN;
        B_coef(i) = NaN;
        source{i} = 'Peneloux';
    else
        A_coef(i) = p1 + p2 * MW(i);
        B_coef(i) = p3 + p4 * MW(i);
        c_cm3(i) = A_coef(i) * T + B_coef(i);
        c(i) = c_cm3(i) / 1e6;
        source{i} = 'Ungerer-Batut';
    end
end

s = c ./ b_i;

c_mix_cm3 = sum(comp .* c_cm3);
c_mix = c_mix_cm3 / 1e6;
s_mix = sum(comp .* s);

mix = struct();
mix.p1 = p1;
mix.p2 = p2;
mix.p3 = p3;
mix.p4 = p4;
mix.A_coef = A_coef;
mix.B_coef = B_coef;
mix.valid = valid;
mix.source = source;
mix.s = s;
mix.c = c;
mix.c_cm3 = c_cm3;
mix.c_mix = c_mix;
mix.c_mix_cm3 = c_mix_cm3;
mix.s_mix = s_mix;

end