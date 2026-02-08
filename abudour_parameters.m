function [c1, zc, omega, Tc, Pc, dipole, found] = abudour_parameters(component_name)

% ABUDOUR_PARAMETERS - Optimized Abudour (2012) parameters for PR-EOS volume translation
% Reference: Abudour et al. (2012), Fluid Phase Equilibria 335, 74-87, Table 1
% Database: {Tc[K], Pc[bar], omega, zc, dipole[D], c1_opt}
% 65 fluids from Table 1

db = struct();

% Inorganic gases
db.CO2 = [304.13, 73.773, 0.22394, 0.2746, 0.000, 0.00652];
db.N2 = [126.19, 33.958, 0.037, 0.2894, 0.000, 0.01386];
db.H2O = [647.14, 220.640, 0.3443, 0.2294, 1.855, -0.01416];
db.O2 = [154.58, 50.430, 0.0222, 0.2879, 0.000, 0.01314];
db.Ar = [150.69, 48.630, -0.00219, 0.2895, 0.000, 0.01433];
db.Xe = [289.73, 58.420, 0.00363, 0.2895, 0.000, 0.01208];
db.F2 = [144.41, 51.724, 0.0449, 0.2880, 0.000, 0.01215];
db.CO = [132.86, 34.935, 0.05, 0.2915, 0.100, 0.01330];
db.SO2 = [430.64, 78.840, 0.2557, 0.2687, 1.600, 0.00333];
db.H2S = [373.10, 90.000, 0.1, 0.2847, 0.970, 0.01091];
db.SF6 = [318.73, 37.546, 0.21, 0.2782, 0.000, 0.00989];
db.NH3 = [405.40, 113.330, 0.25601, 0.2440, 1.470, -0.00777];
db.COS = [378.77, 63.700, 0.0978, 0.2731, 0.715, 0.00877];

% n-Alkanes C1-C20
db.C1 = [190.56, 45.992, 0.011, 0.2863, 0.000, 0.01313];
db.C2 = [305.33, 48.718, 0.0993, 0.2776, 0.000, 0.00993];
db.C3 = [369.83, 42.477, 0.1524, 0.2769, 0.084, 0.00778];
db.nC4 = [425.13, 37.960, 0.201, 0.2738, 0.050, 0.00642];
db.nC5 = [469.70, 33.700, 0.251, 0.2684, 0.370, 0.00434];
db.nC6 = [507.82, 30.340, 0.299, 0.2656, 0.050, 0.00306];
db.nC7 = [540.13, 27.360, 0.349, 0.2632, 0.000, 0.00095];
db.nC8 = [569.32, 24.970, 0.393, 0.2565, 0.070, 0.000195];
db.nC9 = [594.55, 22.810, 0.443, 0.2550, 0.070, -0.00196];
db.nC10 = [617.70, 21.030, 0.488, 0.2501, 0.070, -0.00230];
db.nC12 = [658.10, 18.170, 0.574, 0.2492, 0.000, -0.00428];
db.nC14 = [693.00, 15.700, 0.643, 0.2400, 0.000, -0.00772];
db.nC20 = [768.00, 11.600, 0.9069, 0.2430, 0.000, -0.00913];

% Branched alkanes
db.iC4 = [407.81, 36.290, 0.184, 0.2759, 0.132, 0.00706];
db.iC5 = [460.35, 33.957, 0.2296, 0.2712, 0.100, 0.00689];
db.x2MP = [497.70, 30.400, 0.28, 0.2706, 0.000, 0.00414];  % 2-methylpentane

% Alkenes
db.C2H4 = [282.35, 50.418, 0.0866, 0.2813, 0.000, 0.00995];
db.C3H6 = [364.90, 46.000, 0.142, 0.2800, 0.366, 0.00834];
db.propyne = [402.38, 56.260, 0.204, 0.2751, 0.781, 0.00426];

% Cyclic compounds
db.cyclohexane = [553.64, 40.750, 0.20926, 0.2729, 0.300, 0.00711];
db.cyclopropane = [398.30, 55.797, 0.1305, 0.2743, 0.000, 0.00631];

% Aromatics
db.benzene = [562.05, 48.940, 0.20921, 0.2686, 0.000, 0.00526];
db.toluene = [591.75, 41.263, 0.266, 0.2647, 0.360, 0.00288];

% Alcohols
db.methanol = [512.64, 80.970, 0.5625, 0.2240, 1.700, -0.01264];
db.ethanol = [513.92, 61.480, 0.649, 0.2430, 1.700, -0.00424];
db.propanol = [536.78, 51.750, 0.629, 0.2540, 1.700, -0.00049];
db.isopropanol = [508.30, 47.620, 0.665, 0.2480, 1.700, -0.00220];
db.butanol = [563.05, 44.230, 0.59, 0.2600, 1.800, 0.00073];
db.isobutanol = [547.78, 43.000, 0.59, 0.2580, 1.700, 0.00067];
db.pentanol = [588.15, 39.090, 0.579, 0.2620, 1.700, 0.00224];

% Ketones
db.acetone = [508.20, 47.010, 0.3065, 0.2330, 2.900, -0.00880];

% Refrigerants
db.R11 = [471.11, 44.076, 0.18875, 0.2790, 0.450, 0.00820];
db.R12 = [385.12, 41.361, 0.17948, 0.2764, 0.510, 0.00830];
db.R13 = [302.00, 38.790, 0.1723, 0.2768, 0.510, 0.00914];
db.R14 = [227.51, 37.750, 0.1785, 0.2807, 0.000, 0.01152];
db.R21 = [451.48, 51.812, 0.2061, 0.2701, 1.370, 0.00571];
db.R22 = [369.30, 49.900, 0.22082, 0.2683, 1.458, 0.00425];
db.R32 = [351.26, 57.820, 0.2769, 0.2429, 1.978, -0.00930];
db.R41 = [317.28, 58.970, 0.2012, 0.2400, 1.851, -0.00708];
db.R113 = [487.21, 33.922, 0.25253, 0.2740, 0.803, 0.00614];
db.R114 = [418.83, 32.570, 0.2523, 0.2757, 0.658, 0.00752];
db.R115 = [353.10, 31.200, 0.252, 0.2730, 0.520, 0.00797];
db.R116 = [293.03, 30.480, 0.257, 0.2815, 0.000, 0.00911];
db.R123 = [456.83, 36.618, 0.28192, 0.2681, 1.356, 0.00405];
db.R124 = [395.43, 36.243, 0.2881, 0.2687, 1.469, 0.00507];
db.R125 = [339.17, 36.177, 0.3052, 0.2684, 1.563, 0.00470];
db.R134a = [374.21, 40.593, 0.32684, 0.2601, 2.058, 0.00075];
db.R141b = [477.50, 42.120, 0.22, 0.2706, 2.014, 0.00508];
db.R142b = [410.26, 40.550, 0.232, 0.2679, 2.140, 0.00347];
db.R143a = [345.86, 37.610, 0.2615, 0.2550, 2.340, -0.00221];
db.R152a = [386.41, 45.168, 0.27521, 0.2520, 2.262, -0.00359];
db.R218 = [345.02, 26.400, 0.317, 0.2755, 0.140, 0.00825];
db.R236ea = [412.44, 33.564, 0.3794, 0.2641, 1.129, 0.00343];

% Alias map: alternative names -> canonical name
aliases = struct();
aliases.carbon_dioxide = 'CO2';
aliases.nitrogen = 'N2';
aliases.water = 'H2O';
aliases.oxygen = 'O2';
aliases.argon = 'Ar';
aliases.xenon = 'Xe';
aliases.fluorine = 'F2';
aliases.carbon_monoxide = 'CO';
aliases.sulfur_dioxide = 'SO2';
aliases.hydrogen_sulfide = 'H2S';
aliases.sulfur_hexafluoride = 'SF6';
aliases.ammonia = 'NH3';
aliases.carbonyl_sulfide = 'COS';

aliases.methane = 'C1';
aliases.ethane = 'C2';
aliases.propane = 'C3';
aliases.butane = 'nC4';
aliases.n_butane = 'nC4';
aliases.pentane = 'nC5';
aliases.n_pentane = 'nC5';
aliases.hexane = 'nC6';
aliases.n_hexane = 'nC6';
aliases.heptane = 'nC7';
aliases.n_heptane = 'nC7';
aliases.octane = 'nC8';
aliases.n_octane = 'nC8';
aliases.nonane = 'nC9';
aliases.n_nonane = 'nC9';
aliases.decane = 'nC10';
aliases.n_decane = 'nC10';
aliases.dodecane = 'nC12';
aliases.n_dodecane = 'nC12';
aliases.tetradecane = 'nC14';
aliases.n_tetradecane = 'nC14';
aliases.eicosane = 'nC20';
aliases.n_eicosane = 'nC20';

aliases.isobutane = 'iC4';
aliases.i_butane = 'iC4';
aliases.isopentane = 'iC5';
aliases.i_pentane = 'iC5';
aliases.x2_methylpentane = 'x2MP';

aliases.ethylene = 'C2H4';
aliases.propylene = 'C3H6';

aliases.x1_propanol = 'propanol';
aliases.x2_propanol = 'isopropanol';
aliases.x1_butanol = 'butanol';
aliases.x1_pentanol = 'pentanol';

aliases.CF4 = 'R14';
aliases.tetrafluoromethane = 'R14';
aliases.C2F6 = 'R116';
aliases.hexafluoroethane = 'R116';

% Handle 'list' command
if strcmpi(component_name, 'list')
    fprintf('\n=== Abudour (2012) Table 1 Database ===\n\n');
    fprintf('%-20s %8s %8s %8s %8s %8s %10s\n', ...
            'Component', 'Tc[K]', 'Pc[bar]', 'omega', 'zc', 'dipole', 'c1');
    fprintf('%s\n', repmat('-', 1, 80));
    
    names = fieldnames(db);
    for i = 1:length(names)
        d = db.(names{i});
        fprintf('%-20s %8.2f %8.3f %8.4f %8.4f %8.3f %10.5f\n', ...
                names{i}, d(1), d(2), d(3), d(4), d(5), d(6));
    end
    
    fprintf('\nGeneralized correlation: c1 = 0.4266*Zc - 0.1101 (Eq. 11)\n\n');
    
    c1 = []; zc = []; omega = []; Tc = []; Pc = []; dipole = []; found = false;
    return;
end

% Normalize input name
name = lower(strtrim(component_name));
name = strrep(name, '-', '_');
name = strrep(name, ' ', '_');
name = strrep(name, '1-', 'x1_');
name = strrep(name, '2-', 'x2_');

% Check aliases first
if isfield(aliases, name)
    name = aliases.(name);
end

% Try direct lookup (case-insensitive for db keys)
found = false;
dbnames = fieldnames(db);
for i = 1:length(dbnames)
    if strcmpi(name, dbnames{i})
        d = db.(dbnames{i});
        Tc = d(1);
        Pc = d(2);
        omega = d(3);
        zc = d(4);
        dipole = d(5);
        c1 = d(6);
        found = true;
        return;
    end
end

% Not found
c1 = []; zc = []; omega = []; Tc = []; Pc = []; dipole = [];

end