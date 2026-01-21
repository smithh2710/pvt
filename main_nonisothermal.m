function [comp_h, press_h, temp_h, pressbub_h, pressdew_h] = main_nonisothermal(h, h_ref, comp_ref, press_ref, temp_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref)

%   RT ln(φᵢʰ zᵢʰ Pʰ) - RT ln(φᵢʰ° zᵢʰ° Pʰ°) = Mᵢg(h-h°) - Mᵢ(H^abs/M - H̃ᵢ^abs/Mᵢ)(ΔT/T)
%
% Rearranged:
%   ln(fᵢʰ) = ln(fᵢʰ°) + [Mᵢg(h-h°)]/(R·T_h) - [Mᵢ(H/M - H̃ᵢ/Mᵢ)ΔT]/(R·T_h·T_ref)
%
% Where:
%   - fᵢ = φᵢ zᵢ P is the fugacity of component i
%   - H^abs/M is the mixture absolute specific enthalpy [J/g]
%   - H̃ᵢ^abs/Mᵢ is the partial molar specific enthalpy of component i [J/g]
%   - ΔT = T(h) - T(h°) is the temperature difference
%   - dT/dh is the vertical temperature gradient (typically 0.025 K/m)
%
% =========================================================================
% INPUTS:
%   h           : Target depth [m]
%   h_ref       : Reference depth [m]
%   comp_ref    : Reference composition [mole fractions]
%   press_ref   : Reference pressure [Pa]
%   temp_ref    : Reference temperature [K]
%   dTdh        : Temperature gradient [K/m] (positive = T increases with depth)
%                 Typical value: 0.025 K/m (from paper)
%   Pc          : Critical pressures [Pa]
%   Tc          : Critical temperatures [K]
%   acentric    : Acentric factors [-]
%   BIP         : Binary interaction parameter matrix
%   M_gmol      : Molecular weights [g/mol]
%   Cp_coeffs   : Ideal gas Cp coefficients [n x 4], Cp in [J/(mol·K)]
%                 Cp = C1 + C2*T + C3*T^2 + C4*T^3
%   H_ig_ref    : Ideal gas enthalpy at T_ref=273.15 K [J/mol]
%                 Convert from Table 6: H_ig_ref = (H^ig/(M*R)) * M * R
%
% OUTPUTS:
%   comp_h     : Composition at depth h [mole fractions]
%   press_h    : Pressure at depth h [Pa]
%   temp_h     : Temperature at depth h [K]
%   pressbub_h : Bubble point pressure at depth h [Pa]
%   pressdew_h : Dew point pressure at depth h [Pa]
%
% =========================================================================
% USAGE EXAMPLE (Reservoir 1):
%
%   % Table 6 values: H^ig/(M*R) in [K/g]
%   H_ig_per_mass_R = [-20; 20; 0; 7.5; 15; 17; 17; 25; 25; 33; ...];
%   
%   % Convert to J/mol
%   R = 8.314462618;
%   H_ig_ref = H_ig_per_mass_R .* M_gmol * R;
%   
%   % Temperature gradient from paper
%   dTdh = 0.025;  % K/m
%   
%   [comp_h, P_h, T_h, Pb, Pd] = main_nonisothermal(h, h_ref, comp, P_ref, ...
%       T_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref);
%
% =========================================================================
% SIGN CONVENTION:
%   - Depth h increases downward (deeper = larger h)
%   - dTdh > 0 means temperature increases with depth (normal geothermal)
%   - ΔT = dTdh * (h - h_ref)
%     * If h > h_ref (deeper): ΔT > 0 (warmer)
%     * If h < h_ref (shallower): ΔT < 0 (cooler)
%
% =========================================================================


    R = 8.3144598;        % Universal gas constant [J/(mol·K)]
    g = 9.80665;          % Gravitational acceleration [m/s²]
    

    n = length(comp_ref);
    

    comp_ref = comp_ref(:);
    comp_ref = comp_ref / sum(comp_ref) ; 
    
    M_gmol = M_gmol(:);
    Pc = Pc(:);
    Tc = Tc(:);
    acentric = acentric(:);
    H_ig_ref = H_ig_ref(:);
    
    M_kgmol = M_gmol / 1000;
    

    delta_h = h - h_ref;
    temp_h = temp_ref + dTdh * delta_h;
    delta_T = temp_h - temp_ref;

    tol = 1e-10;
    maxiter = 1500;
    

    [fugcoef_ref, ~] = fugacitycoef_multicomp(comp_ref, press_ref, temp_ref, Pc, Tc, acentric, BIP);
    
    
    f_ref = fugcoef_ref .* comp_ref * press_ref;
    
    [~, H_mix_ref, ~, ~, H_abs_specific_ref] = calculate_absolute_enthalpy(temp_ref, press_ref, comp_ref, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref);
    
   
    M_mix_ref = sum(comp_ref .* M_gmol);  % Mixture molecular weight [g/mol]
    H_specific_mix_ref = H_mix_ref / M_mix_ref;  % H^abs/M [J/g]
    
    %% Calculate Target Fugacity Using Haase Equation (Eq. 6)
    %
    % ln(fᵢʰ) = ln(fᵢʰ°) + [Mᵢ·g·(h-h°)]/(R·T_h) - [Mᵢ·(H/M - H̃ᵢ/Mᵢ)·ΔT]/(R·T_h·T_ref)
    %
    % Gravitational term: drives heavy components to bottom
    % Thermal term: components with H̃ᵢ/Mᵢ > H/M prefer warmer region (bottom)
    
    f_h_target = zeros(n, 1);
    
    for i = 1:n
        % Gravitational term (same as isothermal)
        grav_term = (M_kgmol(i) * g * delta_h) / (R * temp_h);
        
        % Thermal term from Haase model
        % (H^abs/M - H̃ᵢ^abs/Mᵢ) where both are in [J/g]
        enthalpy_diff = H_specific_mix_ref - H_abs_specific_ref(i);  % [J/g]
        
        % Full thermal term: Mᵢ[g/mol] * enthalpy_diff[J/g] * ΔT[K] / (R * T_h * T_ref)
        % Units: [g/mol] * [J/g] * [K] / ([J/(mol·K)] * [K] * [K]) = dimensionless
        thermal_term = (M_gmol(i) * enthalpy_diff * delta_T) / (R * temp_h * temp_ref);
        
        % Target fugacity at depth h
        ln_f_h = log(f_ref(i)) + grav_term - thermal_term;
        f_h_target(i) = exp(ln_f_h);
    end
    
    %% Solve for Composition and Pressure at Depth h
    
    % Initial guess
    initial_guess = [comp_ref; press_ref];
    
    % Define residual function
    residual_fun = @(x) residual_haase(x, f_h_target, temp_h, Pc, Tc, acentric, BIP, n);
    
    % Solve using fsolve
    options = optimoptions('fsolve', ...
        'Display', 'none', ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12, ...
        'OptimalityTolerance', 1e-12, ...
        'MaxIterations', maxiter, ...
        'MaxFunctionEvaluations', maxiter * (n+1));
    
    [solution, ~, exitflag] = fsolve(residual_fun, initial_guess, options);
    
    if exitflag <= 0
        warning('main_nonisothermal: fsolve did not converge (exitflag = %d)', exitflag);
    end
    
    % Extract solution
    comp_h = solution(1:n);
    press_h = solution(end);
    
    % Ensure valid composition
    comp_h = max(comp_h, 0);
    comp_h = comp_h / sum(comp_h);
    
    %% Calculate Saturation Pressures at Depth h
    
    % Bubble point pressure
    try
        pressbub_ini = 260e5;
        if isnan(pressbub_ini) || pressbub_ini <= 0
            pressbub_ini = press_h;
        end
        [pressbub_h, ~] = pressbub_multicomp_newton(comp_h, pressbub_ini, temp_h, Pc, Tc, acentric, BIP, tol, maxiter);
        if ~isreal(pressbub_h) || pressbub_h <= 0 || ~isfinite(pressbub_h)
            pressbub_h = NaN;
        end
    catch
        pressbub_h = NaN;
    end
    
    % Dew point pressure
    try
        pressdew_ini = 260e5;
        if isnan(pressdew_ini) || pressdew_ini <= 0
            pressdew_ini = press_h;
        end
        [pressdew_h, ~] = pressdew_multicomp_newton(comp_h, pressdew_ini, temp_h, Pc, Tc, acentric, BIP, tol, maxiter);
        if ~isreal(pressdew_h) || pressdew_h <= 0 || ~isfinite(pressdew_h)
            pressdew_h = NaN;
        end
    catch
        pressdew_h = NaN;
    end

end



function F = residual_haase(x, f_target, T, Pc, Tc, acentric, BIP, nc)
% RESIDUAL_HAASE - Residual equations for non-isothermal compositional grading
%
% Equations:
%   F(1:n) = φᵢ(z,P,T) * zᵢ * P - fᵢ_target  (fugacity matching)
%   F(n+1) = Σzᵢ - 1                          (mole fraction constraint)

    % Extract variables
    z = x(1:nc);
    P = x(end);
    
    % Ensure positive values
    z = max(z, 1e-15);
    P = max(P, 1e5);
    
    % Normalize for fugacity calculation
    z_norm = z / sum(z);
    
    % Calculate fugacity coefficients
    [phi, ~] = fugacitycoef_multicomp(z_norm, P, T, Pc, Tc, acentric, BIP);
    
    % Current fugacity
    f_current = phi .* z_norm * P;
    
    % Residual equations
    F = zeros(nc+1, 1);
    F(1:nc) = f_current - f_target;
    F(nc+1) = sum(z) - 1;

end
