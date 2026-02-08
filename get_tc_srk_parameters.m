function [params, found] = get_tc_srk_parameters(comp_name)
% GET_TC_SRK_PARAMETERS - Lookup table for tc-SRK (Pina-Martinez 2018)
%
% Returns L, M, N (Twu91 alpha) and c (volume translation) for known components
%
% USAGE:
%   [params, found] = get_tc_srk_parameters('H2O')
%   [params, found] = get_tc_srk_parameters('C1')
%
% OUTPUT:
%   params.L    - Twu91 L parameter
%   params.M    - Twu91 M parameter
%   params.N    - Twu91 N parameter
%   params.c    - Volume translation [cm3/mol]
%   params.name - Matched component name
%   found       - true if component found in table
%
% REFERENCES:
%   Pina-Martinez et al. (2018) J. Chem. Eng. Data 63, 3980-3988
%   Table S2: tc-RK (SRK) parameters

% tc-SRK parameters from Pina-Martinez (2018) Supporting Information Table S2
% Format: {Name, L, M, N, c [cm3/mol]}
table = {
    % Inorganics
    'Nitrogen',         0.1901,  0.8900,  2.0107,   1.3475;
    'N2',               0.1901,  0.8900,  2.0107,   1.3475;
    'Carbon Dioxide',   0.2806,  0.8684,  2.2782,   4.1585;
    'CO2',              0.2806,  0.8684,  2.2782,   4.1585;
    'Hydrogen Sulfide', 0.1748,  0.8686,  2.2761,   3.0181;
    'H2S',              0.1748,  0.8686,  2.2761,   3.0181;
    'Water',            0.4171,  0.8758,  2.1818,   8.9670;
    'H2O',              0.4171,  0.8758,  2.1818,   8.9670;
    
    % n-Alkanes C1-C10
    'Methane',          0.2170,  0.9082,  1.8172,   2.0509;
    'C1',               0.2170,  0.9082,  1.8172,   2.0509;
    'Ethane',           0.2968,  0.8812,  1.7252,   4.6079;
    'C2',               0.2968,  0.8812,  1.7252,   4.6079;
    'Propane',          0.5427,  0.8811,  1.1904,   7.7000;
    'C3',               0.5427,  0.8811,  1.1904,   7.7000;
    'n-Butane',         0.3515,  0.8609,  1.8323,  11.0274;
    'Butane',           0.3515,  0.8609,  1.8323,  11.0274;
    'nC4',              0.3515,  0.8609,  1.8323,  11.0274;
    'n-Pentane',        0.2950,  0.8513,  2.2388,  16.2371;
    'Pentane',          0.2950,  0.8513,  2.2388,  16.2371;
    'nC5',              0.2950,  0.8513,  2.2388,  16.2371;
    'n-Hexane',         0.2982,  0.8491,  2.4179,  22.0561;
    'Hexane',           0.2982,  0.8491,  2.4179,  22.0561;
    'nC6',              0.2982,  0.8491,  2.4179,  22.0561;
    'n-Heptane',        0.3529,  0.8368,  2.2644,  27.8526;
    'Heptane',          0.3529,  0.8368,  2.2644,  27.8526;
    'nC7',              0.3529,  0.8368,  2.2644,  27.8526;
    'n-Octane',         0.3597,  0.8321,  2.3840,  34.8744;
    'Octane',           0.3597,  0.8321,  2.3840,  34.8744;
    'nC8',              0.3597,  0.8321,  2.3840,  34.8744;
    'n-Nonane',         0.4046,  0.8255,  2.3023,  41.4328;
    'Nonane',           0.4046,  0.8255,  2.3023,  41.4328;
    'nC9',              0.4046,  0.8255,  2.3023,  41.4328;
    'n-Decane',         0.3786,  0.8255,  2.5740,  48.7183;
    'Decane',           0.3786,  0.8255,  2.5740,  48.7183;
    'nC10',             0.3786,  0.8255,  2.5740,  48.7183;
    
    % Branched alkanes
    'Isobutane',        0.6853,  0.8842,  1.0305,  10.8558;
    'iC4',              0.6853,  0.8842,  1.0305,  10.8558;
    '2-Methylpropane',  0.6853,  0.8842,  1.0305,  10.8558;
    'Isopentane',       0.2507,  0.8548,  2.3951,  13.8762;
    'iC5',              0.2507,  0.8548,  2.3951,  13.8762;
    '2-Methylbutane',   0.2507,  0.8548,  2.3951,  13.8762;
};

params = struct('L', [], 'M', [], 'N', [], 'c', [], 'name', '');
found = false;

if isempty(comp_name)
    return;
end

comp_name = strtrim(comp_name);

for i = 1:size(table, 1)
    if strcmpi(comp_name, table{i, 1})
        params.L = table{i, 2};
        params.M = table{i, 3};
        params.N = table{i, 4};
        params.c = table{i, 5};
        params.name = table{i, 1};
        found = true;
        return;
    end
end

end