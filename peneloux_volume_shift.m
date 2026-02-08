function result = peneloux_volume_shift(Pc, Tc, omega, comp, Vc, eos_type)
% PENELOUX_VOLUME_SHIFT - Generalized Peneloux volume translation
%
% References:
%   SRK: Peneloux et al. (1982), Fluid Phase Equilibria, 8, 7-23
%   PR:  Jhaveri & Youngren (1988), SPE Reservoir Engineering, 3(3), 1033-1040
%
% INPUTS:
%   Pc       - Critical pressure vector [Pa]
%   Tc       - Critical temperature vector [K]
%   omega    - Acentric factor vector
%   comp     - Mole fraction vector
%   Vc       - Critical volume vector [cm³/mol] (used for mixing rule)
%   eos_type - 'PR' or 'SRK' (optional, default = 'PR')
%
% OUTPUT:
%   result.c_i       - Component volume shifts [cm³/mol]
%   result.c_mix_cm3 - Mixture volume shift [cm³/mol]
%   result.c_mix     - Mixture volume shift [m³/mol]
%   result.eos_type  - EOS type used

R = 8.314;  % J/(mol·K)
nc = length(Pc);

if nargin < 6 || isempty(eos_type)
    eos_type = 'PR';
end

eos_type = upper(eos_type);

% EOS-specific coefficients
% Formula: c = coef1 × (coef2 - Zra) × R × Tc / Pc
switch eos_type
    case 'SRK'
        % Peneloux et al. (1982)
        % c = 0.40768 × (0.29441 - Zra) × R × Tc / Pc
        coef1 = 0.40768;
        coef2 = 0.29441;
    case 'PR'
        % Jhaveri & Youngren (1988)
        % c = 0.50033 × (0.25969 - Zra) × R × Tc / Pc
        coef1 = 0.50033;
        coef2 = 0.25969;
    otherwise
        error('Unknown EOS type: %s. Use ''PR'' or ''SRK''', eos_type);
end

% Rackett compressibility factor
Zra = 0.29056 - 0.08775 * omega;

% Component volume shifts [m³/mol]
c_i_m3 = coef1 * (coef2 - Zra) .* R .* Tc ./ Pc;

% Convert to cm³/mol
c_i = c_i_m3 * 1e6;

% Mixing rule: Vc^(2/3) weighted
Vc_23 = Vc.^(2/3);
c_mix_cm3 = sum(comp .* c_i .* Vc_23) / sum(comp .* Vc_23);

result.c_i = c_i;
result.c_mix_cm3 = c_mix_cm3;
result.c_mix = c_mix_cm3 * 1e-6;  % m³/mol
result.eos_type = eos_type;

end