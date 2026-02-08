function [s, c_cm3] = jhaveri_youngren_volume_shift(Pc, Tc, acentric, MW)
% Jhaveri & Youngren (1988) volume translation for PR-EOS
%
% INPUTS:
%   Pc       : Critical pressures [bar]
%   Tc       : Critical temperatures [K]
%   acentric : Acentric factors [-]
%   MW       : Molecular weights [g/mol]
%
% OUTPUTS:
%   s     : Dimensionless shift factors (c/b) [-]
%   c_cm3 : Volume translation [cm³/mol]

    R = 8.3144598;  % J/(mol·K) = Pa·m³/(mol·K)

    nc = length(Pc);
    Pc = Pc(:);
    Tc = Tc(:);
    acentric = acentric(:);
    MW = MW(:);

    % PR covolume b [cm³/mol]
    % b = 0.07780 * R * Tc / Pc  (with Pc in Pa, gives m³/mol, then * 1e6)
    b = 0.07780 * R * Tc ./ (Pc * 1e5) * 1e6;

    % =====================================================================
    % Hardcoded s = c/b for defined components (JY 1988, Table 1)
    % Rows: [MW_target, acentric_target, s_value]
    % acentric used to disambiguate isomers (iC4/nC4, iC5/nC5)
    % =====================================================================
    %         MW      omega     s = c/b
    sTable = [
        28.0,   0.040,  -0.1600    % N2
        44.0,   0.225,  -0.0800    % CO2
        16.0,   0.008,  -0.1540    % C1
        30.1,   0.098,  -0.1002    % C2
        44.1,   0.152,  -0.0850    % C3
        58.1,   0.176,  -0.0794    % iC4
        58.1,   0.193,  -0.0641    % nC4
        72.2,   0.227,  -0.0435    % iC5
        72.2,   0.251,  -0.0418    % nC5
        86.2,   0.296,  -0.0148    % C6
    ];

    % C7+ correlation: c = (1 - d/MW^e) * b
    % n-alkane defaults (update for different PNA)
    d_heavy = 2.258;
    e_heavy = 0.1823;

    % =====================================================================
    % Assign s values
    % =====================================================================
    s = zeros(nc, 1);

    for i = 1:nc
        matched = false;
        for j = 1:size(sTable, 1)
            if abs(MW(i) - sTable(j,1)) < 1.0 && abs(acentric(i) - sTable(j,2)) < 0.05
                s(i) = sTable(j,3);
                matched = true;
                break;
            end
        end

        if ~matched
            s(i) = 1 - d_heavy / MW(i)^e_heavy;
        end
    end

    c_cm3 = s .* b;

end 