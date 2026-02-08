function [s, c, mix] = baled_volume_shift(T, Pc, Tc, acentric, M, components, comp, Vc, eos_type, verbose)
%BALED_VOLUME_SHIFT Temperature-dependent volume translation for PR/SRK EOS
%
% Reference: Baled et al. (2012), Fluid Phase Equilibria 317, 65-76
%            "Prediction of hydrocarbon densities at extreme conditions
%             using volume-translated SRK and PR equations of state fit
%             to high temperature, high pressure PVT data"
%
% Applicable range: 278-533 K, 7-276 MPa (HTHP conditions)
%
% Volume translation equation (Eq. 10):
%   c = A + B * Tr   [cm³/mol]
%
% For KNOWN components: A and B from Table 1
% For PSEUDO/UNKNOWN: A and B from generalized correlation (Eq. 14, Table 2)
%
% Generalized correlation:
%   A = k0 + k1*exp(-1/(k2*M*ω)) + k3*exp(-1/(k4*M*ω)) + k5*exp(-1/(k6*M*ω))
%   B = k0 + k1*exp(-1/(k2*M*ω)) + k3*exp(-1/(k4*M*ω)) + k5*exp(-1/(k6*M*ω))
%
% Performance (from paper):
%   HTHP VT-SRK: Overall MAPD = 1.47%
%   HTHP VT-PR:  Overall MAPD = 2.01%
%
% Inputs:
%   T          : Temperature [K]
%   Pc         : Critical pressures [Pa]
%   Tc         : Critical temperatures [K]
%   acentric   : Acentric factors [-]
%   M          : Molecular weights [g/mol]
%   components : Cell array of component names
%   comp       : (Optional) Mole fractions [-]
%   Vc         : (Optional) Critical volumes [cm³/mol]
%   eos_type   : (Optional) 'PR' (default) or 'SRK'
%   verbose    : (Optional) Print parameter table (default: false)
%
% Outputs:
%   s   : Dimensionless volume shift [-]
%   c   : Volume shift [m³/mol]
%   mix : Mixture properties structure

R = 8.3144598;

Pc = Pc(:);
Tc = Tc(:);
acentric = acentric(:);
M = M(:);

ncomp = length(Pc);

if nargin < 6 || isempty(components)
    components = cell(ncomp, 1);
elseif ischar(components)
    components = {components};
end

if nargin < 7
    comp = [];
end

if nargin < 8
    Vc = [];
end

if nargin < 9 || isempty(eos_type)
    eos_type = 'PR';
end

if nargin < 10 || isempty(verbose)
    verbose = false;
end

switch upper(eos_type)
    case 'PR'
        Omega_b = 0.0778;
    case 'SRK'
        Omega_b = 0.08664;
    otherwise
        error('Unknown eos_type: %s. Use ''PR'' or ''SRK''.', eos_type);
end

[A, B, source] = get_baled_parameters(components, M, acentric, ncomp, eos_type, verbose);

Tr = T ./ Tc;

c = zeros(ncomp, 1);
s = zeros(ncomp, 1);
c_cm3 = zeros(ncomp, 1);

for i = 1:ncomp
    c_cm3(i) = A(i) + B(i) * Tr(i);
    c(i) = c_cm3(i) * 1e-6;
    b_i = Omega_b * R * Tc(i) / Pc(i);
    s(i) = c(i) / b_i;
end

mix = struct();
mix.A = A;
mix.B = B;
mix.c_cm3 = c_cm3;
mix.source = source;
mix.eos_type = eos_type;

if ~isempty(comp)
    comp = comp(:);
    comp = comp / sum(comp);
    
    A_mix = sum(comp .* A);
    B_mix = sum(comp .* B);
    
    if ~isempty(Vc)
        Vc = Vc(:);
        Vc_23 = Vc.^(2/3);
        seta = (comp .* Vc_23) / sum(comp .* Vc_23);
        Tc_mix = sum(seta .* Tc);
        Vc_mix = sum(seta .* Vc);
        mix.seta = seta;
        mix.Vc_mix = Vc_mix;
    else
        seta = comp;
        Tc_mix = sum(comp .* Tc);
    end
    
    Tr_mix = T / Tc_mix;
    c_mix_cm3 = A_mix + B_mix * Tr_mix;
    c_mix = c_mix_cm3 * 1e-6;
    
    mix.c_mix = c_mix;
    mix.c_mix_cm3 = c_mix_cm3;
    mix.Tc_mix = Tc_mix;
    mix.Tr_mix = Tr_mix;
    mix.A_mix = A_mix;
    mix.B_mix = B_mix;
end

end


function [A, B, source] = get_baled_parameters(components, M, acentric, ncomp, eos_type, verbose)

switch upper(eos_type)
    case 'PR'
        table1 = {
            'C1',          -3.047,   -0.610
            'C2',          -3.328,   -3.189
            'C3',          -3.328,   -3.189
            'nC5',          7.181,  -13.89
            'Cyclohexane',  3.864,  -15.02
            'nC7',         11.24,   -14.57
            'nC8',         20.70,   -23.73
            'iC8',          7.824,  -19.51
            'Cyclooctane',  9.066,  -20.72
            'nC10',        33.71,   -30.91
            'nC13',        62.23,   -45.39
            'nC16',        88.55,   -55.34
            'nC18',       109.0,    -72.80
            'nC20',       116.5,    -60.70
            'nC30',       250.3,   -150.6
            'nC40',       750.5,   -246.9
            'Benzene',      2.074,   -8.227
            'Toluene',     12.17,   -15.37
        };

        kA = [-4.1034, 31.723, 0.0531, 188.68, 0.0057, 20196, 0.0003];
        kB = [-0.3489, -28.547, 0.0687, -817.73, 0.0007, -65.067, 0.0076];
        
    case 'SRK'
       table1 = {
        'C1',          0.233,   -0.420
        'C2',          2.977,   -1.225
        'C3',          2.977,   -1.225
        'nC5',        17.95,   -12.39
        'Cyclohexane',13.52,   -11.65
        'nC7',        26.21,   -11.82
        'nC8',        36.80,   -20.15
        'iC8',        23.92,   -17.46
        'Cyclooctane',23.48,   -19.22
        'nC10',       54.85,   -26.90
        'nC13',       90.21,   -38.01
        'nC16',      127.5,    -52.69
        'nC18',      155.1,    -73.00
        'nC20',      169.4,    -62.91
        'nC30',      325.8,   -146.7
        'nC40',      881.2,   -201.1
        'Benzene',    11.51,    -6.490
        'Toluene',    20.57,   -12.66
    };
        kA = [0.2300, 46.843, 0.0571, 23161, 0.0003, 267.40, 0.0053];
        kB = [-0.3471, -29.748, 0.0644, -347.04, 0.0010, -88.547, 0.0048];
end

A = zeros(ncomp, 1);
B = zeros(ncomp, 1);
source = cell(ncomp, 1);

if verbose
    fprintf('\n--- Baled Volume Shift (%s): Parameter Lookup ---\n', eos_type);
    fprintf('%-12s %12s %12s %12s\n', 'Component', 'A [cm³/mol]', 'B [cm³/mol]', 'Source');
    fprintf('%s\n', repmat('-', 1, 52));
end

for i = 1:ncomp
    found = false;
    comp_name = '';
    
    if i <= length(components) && ~isempty(components{i})
        comp_name = strtrim(components{i});
        
        for j = 1:size(table1, 1)
            if strcmpi(comp_name, table1{j, 1})
                A(i) = table1{j, 2};
                B(i) = table1{j, 3};
                source{i} = 'Table1';
                found = true;
                break;
            end
        end
    end
    
    if ~found
        Mw = M(i) * acentric(i);
        
        if Mw < 0.01
            Mw = 0.01;
        end
        
        A(i) = kA(1) + kA(2)*exp(-1/(kA(3)*Mw)) + ...
               kA(4)*exp(-1/(kA(5)*Mw)) + kA(6)*exp(-1/(kA(7)*Mw));
        
        B(i) = kB(1) + kB(2)*exp(-1/(kB(3)*Mw)) + ...
               kB(4)*exp(-1/(kB(5)*Mw)) + kB(6)*exp(-1/(kB(7)*Mw));
        
        source{i} = 'Generalized';
        
        if isempty(comp_name)
            comp_name = sprintf('Pseudo%d', i);
        end
    end
    
    if verbose
        fprintf('%-12s %12.4f %12.4f %12s\n', comp_name, A(i), B(i), source{i});
    end
end

if verbose
    fprintf('%s\n', repmat('-', 1, 52));
end

end