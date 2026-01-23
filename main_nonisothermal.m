function [comp_h, press_h, temp_h, pressbub_h, pressdew_h] = main_nonisothermal(h, h_ref, comp_ref, press_ref, temp_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref)

% INPUTS:
%   h           : Target depth [m]
%   h_ref       : Reference depth [m]
%   comp_ref    : Reference composition [mole fractions]
%   press_ref   : Reference pressure [Pa]
%   temp_ref    : Reference temperature [K]
%   dTdh        : Temperature gradient [K/m] (positive = T increases with depth)
%   Pc          : Critical pressures [Pa]
%   Tc          : Critical temperatures [K]
%   acentric    : Acentric factors [-]
%   BIP         : Binary interaction parameter matrix
%   M_gmol      : Molecular weights [g/mol]
%   Cp_coeffs   : Ideal gas Cp coefficients [n x 4], Cp in [J/(mol·K)]
%   H_ig_ref    : Ideal gas enthalpy at T_ref=273.15 K [J/mol]
% OUTPUTS:
%   comp_h     : Composition at depth h [mole fractions]
%   press_h    : Pressure at depth h [Pa]
%   temp_h     : Temperature at depth h [K]
%   pressbub_h : Bubble point pressure at depth h [Pa]
%   pressdew_h : Dew point pressure at depth h [Pa]

    R = 8.3144598;        % Universal gas constant [J/(mol·K)]
    g = 9.80665;          % Gravitational acceleration [m/s²]

    n = length(comp_ref);

    % Ensure column vectors and normalize
    comp_ref = comp_ref(:);
    comp_ref = comp_ref / sum(comp_ref);
    M_gmol = M_gmol(:);
    Pc = Pc(:);
    Tc = Tc(:);
    acentric = acentric(:);
    H_ig_ref = H_ig_ref(:);

    M_kgmol = M_gmol / 1000;

    % Temperature at target depth
    delta_h = h - h_ref;
    temp_h = temp_ref + dTdh * delta_h;
    delta_T = temp_h - temp_ref;

    tol = 1e-10;
    maxiter = 1500;

    % Calculate fugacity at reference conditions
    [fugcoef_ref, ~] = fugacitycoef_multicomp(comp_ref, press_ref, temp_ref, Pc, Tc, acentric, BIP);
    f_ref = fugcoef_ref .* comp_ref * press_ref;

    % Initial guess
    initial_guess = [comp_ref; press_ref];

    % Create parameter structure to pass to residual function
    params.f_ref = f_ref;
    params.temp_ref = temp_ref;
    params.temp_h = temp_h;
    params.delta_T = delta_T;
    params.delta_h = delta_h;
    params.Pc = Pc;
    params.Tc = Tc;
    params.acentric = acentric;
    params.BIP = BIP;
    params.M_gmol = M_gmol;
    params.M_kgmol = M_kgmol;
    params.Cp_coeffs = Cp_coeffs;
    params.H_ig_ref = H_ig_ref;
    params.R = R;
    params.g = g;
    params.n = n;

    % Solve using fsolve
    residual_fun = @(x) residual_haase_iterative(x, params);

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
        [pressdew_h, ~] = pressdew_multicomp_newton(comp_h, pressdew_ini, temp_h, Pc, Tc, acentric, BIP, tol, maxiter);
        if ~isreal(pressdew_h) || pressdew_h <= 0 || ~isfinite(pressdew_h)
            pressdew_h = NaN;
        end
    catch
        pressdew_h = NaN;
    end

end


function F = residual_haase_iterative(x, params)
% RESIDUAL_HAASE_ITERATIVE - Residual function with iterative enthalpy calculation
%
% The key improvement: enthalpy terms are recalculated at each iteration
% based on the current trial composition and pressure.
%
% Haase equilibrium:
%   ln(f_i^h) = ln(f_i^ref) + grav_term - thermal_term
%
% where thermal_term depends on H^abs(z^h, P^h, T^h) and H_i^abs(z^h, P^h, T^h)

    % Extract parameters
    f_ref = params.f_ref;
    temp_ref = params.temp_ref;
    temp_h = params.temp_h;
    delta_T = params.delta_T;
    delta_h = params.delta_h;
    Pc = params.Pc;
    Tc = params.Tc;
    acentric = params.acentric;
    BIP = params.BIP;
    M_gmol = params.M_gmol;
    M_kgmol = params.M_kgmol;
    Cp_coeffs = params.Cp_coeffs;
    H_ig_ref = params.H_ig_ref;
    R = params.R;
    g = params.g;
    n = params.n;

    % Extract current trial values
    z_h = x(1:n);
    P_h = x(end);

    % Ensure positive and normalized composition
    z_h = max(z_h, 1e-15);
    z_h_norm = z_h / sum(z_h);

    % Ensure positive pressure
    P_h = max(P_h, 1e5);

    %% Calculate enthalpies at CURRENT iteration conditions (z_h, P_h, T_h)
    % This is the key fix - enthalpies depend on the current state
    [H_abs_i, H_mix, ~, ~, H_abs_specific] = calculate_absolute_enthalpy(temp_h, P_h, z_h_norm, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref);

    % Mixture specific enthalpy [J/g]
    M_mix = sum(z_h_norm .* M_gmol);  % [g/mol]
    H_specific_mix = H_mix / M_mix;   % [J/g]


    %% Calculate fugacity coefficients at current conditions
    [fugcoef_h, ~] = fugacitycoef_multicomp(z_h_norm, P_h, temp_h, Pc, Tc, acentric, BIP);

    % Current fugacity
    f_h_current = fugcoef_h .* z_h_norm * P_h;

    %% Calculate target fugacity using Haase equation with CURRENT enthalpies
    f_h_target = zeros(n, 1);

    for i = 1:n
        % Gravitational term (same as isothermal)
        grav_term = (M_kgmol(i) * g * delta_h) / (R * temp_h);

        % Thermal term using CURRENT enthalpies (this is the fix)
        % H_abs_specific(i) is H_i^abs/M_i in [J/g]
        enthalpy_diff = H_specific_mix - H_abs_specific(i);  % [J/g]

        % thermal_term = M_i[g/mol] * (H/M - H_i/M_i)[J/g] * dT[K] / (R * T_h * T_ref)
        thermal_term = (M_gmol(i) * enthalpy_diff * delta_T) / (R * temp_h * temp_ref);

        % Target fugacity
        ln_f_target = log(f_ref(i)) + grav_term - thermal_term;
        f_h_target(i) = exp(ln_f_target);
    end

    %% Residual equations
    F = zeros(n+1, 1);

    % Fugacity matching: f_i^current - f_i^target = 0
    F(1:n) = f_h_current - f_h_target;

    % Mole fraction constraint: sum(z_i) - 1 = 0
    F(n+1) = sum(z_h) - 1;

end
