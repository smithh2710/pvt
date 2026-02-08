function [s, c, mix] = abudour_volume_shift(comp, press, temp, Pc, Tc, acentric, Vc, MW, components, BIP, eos_type)

R = 8.3144598;

if nargin < 10 || isempty(BIP)
    BIP = zeros(length(Pc));
end
if nargin < 11 || isempty(eos_type)
    eos_type = 'PR';
end

if strcmp(eos_type, 'SRK')
    Zc_EOS = 1/3;
else
    Zc_EOS = 0.3074;
end

nc = length(Pc);
comp = comp(:);
Pc = Pc(:);
Tc = Tc(:);
acentric = acentric(:);
Vc = Vc(:);
MW = MW(:);

c1 = zeros(nc, 1);
Zc = zeros(nc, 1);

for i = 1:nc
    if i <= length(components) && ~isempty(components{i})
        [c1_db, zc_db, ~, ~, ~, ~, found] = abudour_parameters(strtrim(components{i}));
        if found
            c1(i) = c1_db;
            Zc(i) = zc_db;
        else
            [c1(i), Zc(i)] = abudour_generalized_c1(acentric(i));
        end
    else
        [c1(i), Zc(i)] = abudour_generalized_c1(acentric(i));
    end
end

a_coef = zeros(nc, 1);
b_coef = zeros(nc, 1);
for i = 1:nc
    [a_coef(i), b_coef(i)] = calc_ab_EOS(Tc(i), Pc(i), acentric(i), temp, eos_type);
end

a_mix = 0;
b_mix = 0;
for i = 1:nc
    for j = 1:nc
        a_ij = sqrt(a_coef(i) * a_coef(j)) * (1 - BIP(i,j));
        a_mix = a_mix + comp(i) * comp(j) * a_ij;
        b_ij = (b_coef(i) + b_coef(j)) / 2;
        b_mix = b_mix + comp(i) * comp(j) * b_ij;
    end
end

Vc_23 = Vc.^(2/3);
sum_Vc23 = sum(comp .* Vc_23);
if sum_Vc23 > 0
    seta = (comp .* Vc_23) / sum_Vc23;
else
    seta = comp;
end

Tc_mix = sum(seta .* Tc);
Vc_mix = sum(seta .* Vc);
Vc_mix_m3 = Vc_mix * 1e-6;
omega_mix = sum(comp .* acentric);
Zc_mix_est = 0.2905 - 0.085 * omega_mix;
Pc_mix = Zc_mix_est * R * Tc_mix / Vc_mix_m3;
Zc_mix = sum(comp .* Zc);
c1_mix = sum(comp .* c1);
MW_mix = sum(comp .* MW);

[V_EOS, ~] = calc_EOS_volume(press, temp, a_mix, b_mix, R, eos_type);
d = calc_distance_function(V_EOS, temp, a_mix, b_mix, R, Tc_mix, eos_type);
c_base = (R * Tc_mix / Pc_mix) * (c1_mix - (0.004 + c1_mix) * exp(-2 * d));

vc_i_EOS = Zc_EOS * R * Tc ./ Pc;
Vc_mix_EOS = sum(comp .* vc_i_EOS);
Vc_mix_true = Vc_mix_m3;
delta_c = Vc_mix_EOS - Vc_mix_true;

f_d = 0.35 / (0.35 + d);
c_eff =  - ( c_base - delta_c * f_d ) ;

V_VTEOS = V_EOS - c_eff;
rho_mol = 1 / V_VTEOS;
rho_kg = rho_mol * MW_mix / 1000;

c = zeros(nc, 1);
s = zeros(nc, 1);
for i = 1:nc
    [V_i, ~] = calc_EOS_volume(press, temp, a_coef(i), b_coef(i), R, eos_type);
    d_i = calc_distance_function(V_i, temp, a_coef(i), b_coef(i), R, Tc(i), eos_type);
    c_i = (R * Tc(i) / Pc(i)) * (c1(i) - (0.004 + c1(i)) * exp(-2 * d_i));
    Vc_EOS_i = Zc_EOS * R * Tc(i) / Pc(i);
    Vc_true_i = Vc(i) * 1e-6;
    delta_c_i = Vc_EOS_i - Vc_true_i;
    f_d_i = 0.35 / (0.35 + d_i);
    c(i) = c_i - delta_c_i * f_d_i;
    s(i) = c(i) / b_coef(i);
end

mix.c_mix = c_eff;
mix.c_mix_cm3 = c_eff * 1e6;
mix.c_base = c_base;
mix.delta_c = delta_c;
mix.f_d = f_d;
mix.d = d;
mix.c1_mix = c1_mix;
mix.c1 = c1;
mix.Zc = Zc;
mix.Tc_mix = Tc_mix;
mix.Pc_mix = Pc_mix;
mix.Vc_mix = Vc_mix;
mix.Vc_mix_EOS = Vc_mix_EOS * 1e6;
mix.Zc_mix = Zc_mix;
mix.omega_mix = omega_mix;
mix.MW_mix = MW_mix;
mix.V_EOS = V_EOS;
mix.V_VTEOS = V_VTEOS;
mix.rho_mol = rho_mol;
mix.rho_kg = rho_kg;
mix.seta = seta;
mix.eos_type = eos_type;

end


function [a, b] = calc_ab_EOS(Tc, Pc, omega, T, eos_type)

R = 8.3144598;

if strcmp(eos_type, 'SRK')
    Omega_a = 0.42748;
    Omega_b = 0.08664;
    m = 0.48 + 1.574*omega - 0.176*omega^2;
else
    Omega_a = 0.45724;
    Omega_b = 0.07780;
    if omega > 0.49
        m = 0.379642 + 1.48503*omega - 0.164423*omega^2 + 0.016666*omega^3;
    else
        m = 0.37464 + 1.54226*omega - 0.26992*omega^2;
    end
end

b = Omega_b * R * Tc / Pc;
Tr = T / Tc;
alpha = (1 + m * (1 - sqrt(Tr)))^2;
a = Omega_a * (R * Tc)^2 / Pc * alpha;

end


function [V, rho] = calc_EOS_volume(P, T, a, b, R, eos_type)

A = a * P / (R * T)^2;
B = b * P / (R * T);

if strcmp(eos_type, 'SRK')
    p2 = -1;
    p1 = A - B - B^2;
    p0 = -A*B;
else
    p2 = -(1 - B);
    p1 = A - 3*B^2 - 2*B;
    p0 = -(A*B - B^2 - B^3);
end

roots_Z = roots([1, p2, p1, p0]);
real_roots = roots_Z(abs(imag(roots_Z)) < 1e-10);
real_roots = real(real_roots);

valid_roots = real_roots(real_roots > B);
if isempty(valid_roots)
    Z = max(real_roots);
else
    Z = min(valid_roots);
end

V = Z * R * T / P;
rho = 1 / V;

end


function d = calc_distance_function(V, T, a, b, R, Tc, eos_type)

if strcmp(eos_type, 'SRK')
    denom1 = V - b;
    denom2 = V * (V + b);
    dPdV = -R*T / denom1^2 + a * (2*V + b) / denom2^2;
else
    denom1 = V - b;
    denom2 = V*(V + b) + b*(V - b);
    dPdV = -R*T / denom1^2 + a * (2*V + 2*b) / denom2^2;
end

dPdrho = -V^2 * dPdV;
d = abs(dPdrho) / (R * Tc);

end