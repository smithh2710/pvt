function result = pina_martinez_volume_shift(T, Pc, Tc, omega, components, comp)
% PINA_MARTINEZ_VOLUME_SHIFT - Volume translation for tc-SRK EOS
% Reference: Pina-Martinez et al. (2018) J. Chem. Eng. Data 63, 3980-3988
%
% Uses get_tc_srk_parameters() for c lookup (DIPPR-fitted values)
%
% INPUTS:
%   T          - Temperature [K] (kept for interface consistency)
%   Pc         - Critical pressure [Pa]
%   Tc         - Critical temperature [K]
%   omega      - Acentric factor
%   components - Cell array of component NAMES
%   comp       - Mole fractions (optional, defaults to equal)
%
% OUTPUT:
%   result.c_i        - Component c values [cm3/mol]
%   result.c_mix_cm3  - Mixture c [cm3/mol]
%   result.c_mix      - Mixture c [m3/mol]
%   result.source     - Cell array ('Table' or 'Correlation')

R = 8.314472;

nc = length(components);
Pc = Pc(:);
Tc = Tc(:);
omega = omega(:);

if nargin < 6 || isempty(comp)
    comp = ones(nc, 1) / nc;
else
    comp = comp(:);
end

c_i = zeros(nc, 1);
source = cell(nc, 1);

for i = 1:nc
    found = false;
    
    if i <= length(components) && ~isempty(components{i})
        [params, found] = get_tc_srk_parameters(components{i});
        if found
            c_i(i) = params.c;
            source{i} = 'Table';
        end
    end
    
    if ~found
        Zra = 0.29056 - 0.08775 * omega(i);
        c_i(i) = 0.40768 * (0.29441 - Zra) * R * Tc(i) / Pc(i) * 1e6;
        source{i} = 'Correlation';
    end
end

c_mix_cm3 = sum(comp .* c_i);
c_mix = c_mix_cm3 * 1e-6;

result.c_i = c_i;
result.c_mix_cm3 = c_mix_cm3;
result.c_mix = c_mix;
result.source = source;

end