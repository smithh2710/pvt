function [s, c, mix] = chen_li_volume_shift(comp, press, temp, Pc, Tc, acentric, Zc, components, verbose)
% CHEN_LI_VOLUME_SHIFT - Volume translation for SRK-EOS using Chen & Li (2020)
% Reference: Chen & Li (2020), Fluid Phase Equilibria 521, 112724
% Alpha function: Pina-Martinez Twu91 with component-specific L, M, N
%
% Uses get_tc_srk_parameters() for L, M, N lookup

R = 8.314472;

nc = length(comp);
comp = comp(:);
Pc = Pc(:);
Tc = Tc(:);
acentric = acentric(:);
Zc = Zc(:);

if nargin < 8 || isempty(components)
    components = cell(nc, 1);
end

if nargin < 9 || isempty(verbose)
    verbose = false;
end

Zc_SRK = 1/3;
Omega_a = 0.42748;
Omega_b = 0.08664;

%% Chen-Li c1, c2, c3 table (these are VT-specific, not in get_tc_srk_parameters)
table_c123 = {
    'Methane',      0.00195,  0.79540,  2.13497;
    'C1',           0.00195,  0.79540,  2.13497;
    'Ethane',       0.00270,  0.85431,  2.59463;
    'C2',           0.00270,  0.85431,  2.59463;
    'Propane',      0.00492,  0.89221,  2.75570;
    'C3',           0.00492,  0.89221,  2.75570;
    'Butane',       0.00594,  0.93036,  2.60453;
    'nC4',          0.00594,  0.93036,  2.60453;
    'n-Butane',     0.00594,  0.93036,  2.60453;
    'Isobutane',    0.00554,  0.89812,  2.61896;
    'iC4',          0.00554,  0.89812,  2.61896;
    'Pentane',      0.00831,  0.99529,  3.00384;
    'nC5',          0.00831,  0.99529,  3.00384;
    'n-Pentane',    0.00831,  0.99529,  3.00384;
    '2-Methylbutane', 0.00688, 0.89542, 3.07548;
    'iC5',          0.00688,  0.89542,  3.07548;
    'Isopentane',   0.00688,  0.89542,  3.07548;
    'Hexane',       0.00923,  0.86394,  2.07261;
    'nC6',          0.00923,  0.86394,  2.07261;
    'n-Hexane',     0.00923,  0.86394,  2.07261;
    'Heptane',      0.01032,  1.18249,  2.08864;
    'nC7',          0.01032,  1.18249,  2.08864;
    'n-Heptane',    0.01032,  1.18249,  2.08864;
    'Octane',       0.01226,  1.19624,  2.35762;
    'nC8',          0.01226,  1.19624,  2.35762;
    'n-Octane',     0.01226,  1.19624,  2.35762;
    'Nonane',       0.01444,  1.12707,  2.66252;
    'nC9',          0.01444,  1.12707,  2.66252;
    'n-Nonane',     0.01444,  1.12707,  2.66252;
    'Decane',       0.01574,  1.25347,  3.03108;
    'nC10',         0.01574,  1.25347,  3.03108;
    'n-Decane',     0.01574,  1.25347,  3.03108;
    'Cyclohexane',  0.00459,  1.09262,  2.21252;
    'Benzene',      0.00680,  1.00489,  2.55901;
    'Toluene',      0.01024,  1.05282,  2.87184;
    'CO2',          0.00608,  0.92912,  2.65917;
    'Carbon Dioxide', 0.00608, 0.92912, 2.65917;
    'N2',           0.00252,  0.75199,  2.19566;
    'Nitrogen',     0.00252,  0.75199,  2.19566;
    'H2S',          0.00144,  0.97009,  2.45887;
    'Hydrogen Sulfide', 0.00144, 0.97009, 2.45887;
    'Water',        0.02425,  1.30564,  2.17549;
    'H2O',          0.02425,  1.30564,  2.17549;
    'Ammonia',      0.02004,  1.14567,  2.55131;
    'NH3',          0.02004,  1.14567,  2.55131;
};

%% Get L, M, N for each component from centralized function
L = zeros(nc, 1);
M = zeros(nc, 1);
N = zeros(nc, 1);
source_LMN = cell(nc, 1);

for i = 1:nc
    found = false;
    
    if i <= length(components) && ~isempty(components{i})
        [params, found] = get_tc_srk_parameters(components{i});
        if found
            L(i) = params.L;
            M(i) = params.M;
            N(i) = params.N;
            source_LMN{i} = 'Table';
        end
    end
    
    if ~found
        omega_i = acentric(i);
        L(i) = 0.1359 + 0.7535 * omega_i + 0.0611 * omega_i^2;
        M(i) = 0.8787 - 0.2063 * omega_i + 0.1709 * omega_i^2;
        N(i) = 2.0;
        source_LMN{i} = 'Generalized';
    end
end

%% Get c1, c2, c3 for each component
c1 = zeros(nc, 1);
c2 = zeros(nc, 1);
c3 = zeros(nc, 1);
source_c = cell(nc, 1);

for i = 1:nc
    found = false;
    
    if i <= length(components) && ~isempty(components{i})
        comp_name = strtrim(components{i});
        
        for j = 1:size(table_c123, 1)
            if strcmpi(comp_name, table_c123{j, 1})
                c1(i) = table_c123{j, 2};
                c2(i) = table_c123{j, 3};
                c3(i) = table_c123{j, 4};
                source_c{i} = 'Table';
                found = true;
                break;
            end
        end
    end
    
    if ~found
        omega_i = acentric(i);
        Zc_i = Zc(i);
        
        c1(i) = -3.90812e-4 + 0.03274 * omega_i;
        
        if Zc_i < 0.2215
            Zc_i = 0.2215;
        elseif Zc_i > 0.2864
            Zc_i = 0.2864;
        end
        
        c2(i) = 3.06048 - 7.64314 * Zc_i;
        c3(i) = 8.34576 - 21.07619 * Zc_i;
        
        source_c{i} = 'Generalized';
    end
end

%% Calculate EOS parameters with component-specific Twu91 alpha
Tr = temp ./ Tc;
alpha = Tr.^(N.*(M-1)) .* exp(L .* (1 - Tr.^(N.*M)));

a_i = Omega_a * R^2 * Tc.^2 ./ Pc .* alpha;
b_i = Omega_b * R * Tc ./ Pc;

a_mix = 0;
for i = 1:nc
    for j = 1:nc
        a_mix = a_mix + comp(i) * comp(j) * sqrt(a_i(i) * a_i(j));
    end
end
b_mix = sum(comp .* b_i);

%% Solve SRK cubic for Z
A = a_mix * press / (R^2 * temp^2);
B = b_mix * press / (R * temp);

coeffs = [1, -1, (A - B - B^2), -A*B];
Z_roots = roots(coeffs);
Z_roots = Z_roots(imag(Z_roots) == 0 & real(Z_roots) > 0);

if isempty(Z_roots)
    Z = 0.3;
else
    Z = min(real(Z_roots));
end

V_EOS = Z * R * temp / press;

%% For pure component or mixture
if nc == 1
    Tc_use = Tc(1);
    Pc_use = Pc(1);
    Zc_exp = Zc(1);
    seta = 1;
    Vc_mix = Zc(1) * R * Tc(1) / Pc(1);
else
    Vc_est = Zc .* R .* Tc ./ Pc;
    Vc_23 = Vc_est.^(2/3);
    sum_xVc23 = sum(comp .* Vc_23);
    if sum_xVc23 > 1e-15
        seta = (comp .* Vc_23) / sum_xVc23;
    else
        seta = comp;
    end
    Tc_use = sum(seta .* Tc);
    Vc_mix = sum(seta .* Vc_est);
    omega_mix = sum(comp .* acentric);
    Pc_use = (0.2905 - 0.085 * omega_mix) * R * Tc_use / Vc_mix;
    Zc_exp = 0.2905 - 0.085 * omega_mix;
end

%% Distance function (Eq. 5)
term1 = (V_EOS^2 * temp) / (Tc_use * (V_EOS - b_mix)^2);
term2 = (a_mix * (2*V_EOS + b_mix)) / (R * Tc_use * (V_EOS + b_mix)^2);
d_mix = term1 - term2;

%% Volume translation parameters
c1_mix = sum(comp .* c1);
c2_mix = sum(comp .* c2);
c3_mix = sum(comp .* c3);

%% delta_c (Eq. 7)
delta_c = (R * Tc_use / Pc_use) * (Zc_SRK - Zc_exp);

%% Total volume translation (Eq. 9)
c_far = c1_mix * R * Tc_use / Pc_use;
c_near = delta_c / (c2_mix + c3_mix * d_mix);
c_total = c_far + c_near;

V_VTEOS = V_EOS - c_total;

%% Output
c = c1 .* R .* Tc ./ Pc;
s = c ./ b_i;

mix.c_mix = c_total;
mix.c_mix_cm3 = c_total * 1e6;
mix.delta_c = delta_c;
mix.d = d_mix;
mix.c1 = c1;
mix.c2 = c2;
mix.c3 = c3;
mix.c1m = c1_mix;
mix.c2m = c2_mix;
mix.c3m = c3_mix;
mix.L = L;
mix.M = M;
mix.N = N;
mix.alpha = alpha;
mix.Tc_mix = Tc_use;
mix.Pc_mix = Pc_use;
mix.Vc_mix = Vc_mix;
mix.V_EOS = V_EOS;
mix.V_VTEOS = V_VTEOS;
mix.Z = Z;
mix.source_LMN = source_LMN;
mix.source_c = source_c;
mix.seta = seta;
mix.a = a_mix;
mix.b = b_mix;

if verbose
    fprintf('Chen-Li Volume Shift Results:\n');
    fprintf('  V_EOS = %.6e m3/mol\n', V_EOS);
    fprintf('  c_total = %.6e m3/mol (%.4f cm3/mol)\n', c_total, c_total*1e6);
    fprintf('  V_VTEOS = %.6e m3/mol\n', V_VTEOS);
    for i = 1:nc
        fprintf('  Comp %d: L=%.4f, M=%.4f, N=%.4f (%s)\n', i, L(i), M(i), N(i), source_LMN{i});
        fprintf('          c1=%.5f, c2=%.5f, c3=%.5f (%s)\n', c1(i), c2(i), c3(i), source_c{i});
    end
end

end