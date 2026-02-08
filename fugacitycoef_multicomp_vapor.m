%% CALCULATE THE FUGACITY COEFFICIENT AND Z-FACTOR OF MULTI-COMPONENT SYSTEMS
% -------------------------------------------------------------------------
% INPUTS:
%   comp       - composition (mole fractions)
%   press      - pressure [Pa]
%   temp       - temperature [K]
%   pressc     - critical pressure [Pa]
%   tempc      - critical temperature [K]
%   acentric   - acentric factor
%   BIP        - binary interaction parameter matrix
%   eos_type   - 'PR' (default) or 'SRK'
%   components - (optional) cell array of component names for tc-SRK lookup
%
% OUTPUTS:
%   fugcoef    - fugacity coefficient
%   zfactor    - compressibility factor (maximum root for vapor)
% -------------------------------------------------------------------------

function [fugcoef, zfactor] = fugacitycoef_multicomp_vapor(comp, press, temp, pressc, tempc, acentric, BIP, eos_type, components)

if nargin < 8 || isempty(eos_type)
    eos_type = 'PR';
end

if nargin < 9
    components = {};
end

[A, B] = calcab_multicomp(press, temp, pressc, tempc, acentric, eos_type, components);

[Amix, Bmix, Amix2] = calcabmix(comp, A, B, BIP);

zfactor = calczfactor(Amix, Bmix, eos_type);

if (size(zfactor,1) > 1)
    zfactor = max(zfactor);
end

fugcoef = calcfugcoef_multicomp(zfactor, A, B, Amix, Bmix, Amix2, eos_type);

end