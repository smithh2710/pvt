function [stability, x_trial1, x_trial2, TPD_min] = stability_analysis_qnss(comp, press, temp, Pc, Tc, acentric, BIP, tol, maxiter, eos_type)
% ENHANCED STABILITY ANALYSIS USING TANGENT PLANE DISTANCE (TPD)
% Based on Michelsen (1982) and Nghiem & Li (1984) QNSS method
%
% Inputs:
%   comp     : Overall composition (mole fractions)
%   press    : Pressure [Pa]
%   temp     : Temperature [K]
%   Pc       : Critical pressure vector [Pa]
%   Tc       : Critical temperature vector [K]
%   acentric : Acentric factor vector
%   BIP      : Binary interaction parameter matrix
%   tol      : Convergence tolerance (default: 1e-8)
%   maxiter  : Maximum iterations (default: 100)
%   eos_type : 'PR' (default) or 'SRK'
%
% Outputs:
%   stability : 1 = stable (single phase), 2 = unstable (will split)
%   x_trial1  : First trial phase composition (if unstable)
%   x_trial2  : Second trial phase composition (if unstable)
%   TPD_min   : Minimum tangent plane distance value

% Set default values if not provided
if nargin < 8 || isempty(tol)
    tol = 1e-8;
end
if nargin < 9 || isempty(maxiter)
    maxiter = 100;
end
if nargin < 10 || isempty(eos_type)
    eos_type = 'PR';
end

% Normalize composition
comp = comp(:) / sum(comp);
ncomp = length(comp);

% Protect against zero compositions
comp = max(comp, 1e-15);
comp = comp / sum(comp);

% Calculate reference fugacity coefficients for overall composition
[fugcoef_ref, ~] = fugacitycoef_multicomp(comp, press, temp, Pc, Tc, acentric, BIP, eos_type);

% Check for valid fugacity coefficients
if any(~isfinite(fugcoef_ref)) || any(fugcoef_ref <= 0)
    stability = 1;
    x_trial1 = comp;
    x_trial2 = comp;
    TPD_min = 0;
    return;
end

lnfugcoef_ref = log(fugcoef_ref);

% Initialize storage for multiple trial calculations
n_trials = 5;  % 4 Wilson-based + 1 pure component trial
TPD = inf(n_trials, 1);
u_min = zeros(n_trials, ncomp);
x_min = zeros(n_trials, ncomp);
converged_flags = false(n_trials, 1);

%% TRIAL 1-4: Wilson equation based initializations
for trial = 1:4
    % Different initialization strategies
    switch trial
        case 1  % Vapor-like initialization
            lnK = wilson_K_factors(press, temp, Pc, Tc, acentric, 1.0);
        case 2  % Intermediate vapor
            lnK = wilson_K_factors(press, temp, Pc, Tc, acentric, 0.33);
        case 3  % Liquid-like initialization
            lnK = wilson_K_factors(press, temp, Pc, Tc, acentric, -1.0);
        case 4  % Intermediate liquid
            lnK = wilson_K_factors(press, temp, Pc, Tc, acentric, -0.33);
    end
    
    % Perform QNSS iteration for both vapor and liquid-like phases
    [u_vapor, conv_v] = qnss_iteration(comp, lnK, press, temp, Pc, Tc, ...
                                       acentric, BIP, lnfugcoef_ref, tol, maxiter, 'vapor', eos_type);
    [u_liquid, conv_l] = qnss_iteration(comp, lnK, press, temp, Pc, Tc, ...
                                        acentric, BIP, lnfugcoef_ref, tol, maxiter, 'liquid', eos_type);
    
    % Choose the more promising result (higher sum = more likely unstable)
    if sum(u_vapor) > sum(u_liquid) && conv_v
        u_min(trial, :) = u_vapor;
        converged_flags(trial) = conv_v;
    elseif conv_l
        u_min(trial, :) = u_liquid;
        converged_flags(trial) = conv_l;
    elseif conv_v
        u_min(trial, :) = u_vapor;
        converged_flags(trial) = conv_v;
    end
    
    % Calculate TPD if converged
    if converged_flags(trial) && sum(u_min(trial, :)) > 1e-15
        x_min(trial, :) = u_min(trial, :) / sum(u_min(trial, :));
        x_trial = x_min(trial, :)';
        x_trial = max(x_trial, 1e-15);
        x_trial = x_trial / sum(x_trial);
        
        try
            [fugcoef_trial, ~] = fugacitycoef_multicomp(x_trial, press, temp, Pc, Tc, acentric, BIP, eos_type);
            if all(isfinite(fugcoef_trial)) && all(fugcoef_trial > 0)
                TPD(trial) = calculate_tpd(x_min(trial, :)', comp, fugcoef_trial, fugcoef_ref);
            end
        catch
            % Skip this trial
        end
    end
end

%% TRIAL 5: Pure component initialization (most volatile component)
[~, most_volatile] = max(Pc ./ Tc .* sqrt(Pc));

lnK = zeros(ncomp, 1);
for i = 1:ncomp
    if i == most_volatile
        lnK(i) = log((1 - 1e-10) / max(comp(i), 1e-15));
    else
        lnK(i) = log(1e-10 / max((ncomp - 1), 1) / max(comp(i), 1e-15));
    end
end

% Perform stability test
[u_vapor, conv_v] = qnss_iteration(comp, lnK, press, temp, Pc, Tc, ...
                                   acentric, BIP, lnfugcoef_ref, tol, maxiter, 'vapor', eos_type);
[u_liquid, conv_l] = qnss_iteration(comp, lnK, press, temp, Pc, Tc, ...
                                    acentric, BIP, lnfugcoef_ref, tol, maxiter, 'liquid', eos_type);

if sum(u_vapor) > sum(u_liquid) && conv_v
    u_min(5, :) = u_vapor;
    converged_flags(5) = conv_v;
elseif conv_l
    u_min(5, :) = u_liquid;
    converged_flags(5) = conv_l;
elseif conv_v
    u_min(5, :) = u_vapor;
    converged_flags(5) = conv_v;
end

if converged_flags(5) && sum(u_min(5, :)) > 1e-15
    x_min(5, :) = u_min(5, :) / sum(u_min(5, :));
    x_trial = x_min(5, :)';
    x_trial = max(x_trial, 1e-15);
    x_trial = x_trial / sum(x_trial);
    
    try
        [fugcoef_trial, ~] = fugacitycoef_multicomp(x_trial, press, temp, Pc, Tc, acentric, BIP, eos_type);
        if all(isfinite(fugcoef_trial)) && all(fugcoef_trial > 0)
            TPD(5) = calculate_tpd(x_min(5, :)', comp, fugcoef_trial, fugcoef_ref);
        end
    catch
        % Skip
    end
end

%% Determine stability and output results
valid_idx = converged_flags & isfinite(TPD);

if any(valid_idx)
    [TPD_min, idx_min] = min(TPD);
    
    % Find a second distinct minimum if it exists
    TPD_temp = TPD;
    TPD_temp(idx_min) = inf;
    [TPD_second, idx_second] = min(TPD_temp);
    
    % Stability criterion: TPD < -tolerance indicates instability
    if TPD_min < -1e-6
        stability = 2;  % Unstable (will split into phases)
        x_trial1 = x_min(idx_min, :)';
        x_trial1 = max(x_trial1, 1e-15);
        x_trial1 = x_trial1 / sum(x_trial1);
        
        % Find second trial phase that is sufficiently different
        if isfinite(TPD_second) && norm(x_min(idx_second, :)' - x_trial1) > 0.05
            x_trial2 = x_min(idx_second, :)';
            x_trial2 = max(x_trial2, 1e-15);
            x_trial2 = x_trial2 / sum(x_trial2);
        else
            % Use complementary phase estimate
            x_trial2 = estimate_complementary_phase(x_trial1, comp);
        end
    else
        stability = 1;  % Stable (single phase)
        x_trial1 = comp;
        x_trial2 = comp;
    end
else
    % No valid results - assume stable
    stability = 1;
    x_trial1 = comp;
    x_trial2 = comp;
    TPD_min = 0;
end

end

%% QNSS ITERATION FUNCTION
function [u, converged] = qnss_iteration(z, lnK_init, press, temp, Pc, Tc, acentric, BIP, lnfugcoef_z, tol, maxiter, phase_type, eos_type)
% Quasi-Newton Successive Substitution iteration for stability analysis

ncomp = length(z);
lnK = lnK_init(:);

% Initial guess based on phase type
if strcmp(phase_type, 'vapor')
    u = exp(lnK) .* z;
else
    u = z ./ exp(lnK);
end

% Protect against overflow/underflow
u = max(u, 1e-15);
u = min(u, 1e10);

% Normalize
x = u / sum(u);
x = max(x, 1e-15);
x = x / sum(x);

% Calculate initial fugacity coefficients
try
    [fugcoef_x, ~] = fugacitycoef_multicomp(x, press, temp, Pc, Tc, acentric, BIP, eos_type);
catch
    u = z;
    converged = false;
    return;
end

if any(~isfinite(fugcoef_x)) || any(fugcoef_x <= 0)
    u = z;
    converged = false;
    return;
end

lnfugcoef_x = log(fugcoef_x);

% Initial residual
d = log(u) + lnfugcoef_x - log(z) - lnfugcoef_z;
err = norm(d);

% QNSS parameters
sigma = 1.0;
iter = 0;
converged = false;

while err > tol && iter < maxiter
    % Store old values
    d_old = d;
    
    % Update step
    dlnK = -sigma * d;
    
    % Limit step size
    max_step = 3.0;
    dlnK = sign(dlnK) .* min(abs(dlnK), max_step);
    
    % Update lnK
    lnK = lnK + dlnK;
    
    % Update composition
    u = exp(lnK) .* z;
    u = max(u, 1e-15);
    u = min(u, 1e10);
    
    x = u / sum(u);
    x = max(x, 1e-15);
    x = x / sum(x);
    
    % Calculate new fugacity coefficients
    try
        [fugcoef_x, ~] = fugacitycoef_multicomp(x, press, temp, Pc, Tc, acentric, BIP, eos_type);
    catch
        break;
    end
    
    if any(~isfinite(fugcoef_x)) || any(fugcoef_x <= 0)
        break;
    end
    
    lnfugcoef_x = log(fugcoef_x);
    
    % New residual
    d = log(u) + lnfugcoef_x - log(z) - lnfugcoef_z;
    
    % Update sigma using quasi-Newton formula
    dd = d - d_old;
    denominator = dlnK' * dd;
    if abs(denominator) > 1e-10
        sigma_new = abs(-(dlnK' * d_old) / denominator) * sigma;
        sigma = min(max(sigma_new, 0.1), 2.0);  % Limit sigma range
    end
    
    err = norm(d);
    iter = iter + 1;
    
    % Reset sigma periodically
    if mod(iter, 10) == 0
        sigma = 1.0;
    end
end

if iter < maxiter && err <= tol
    converged = true;
end

end

%% WILSON K-FACTOR CALCULATION
function lnK = wilson_K_factors(press, temp, Pc, Tc, acentric, factor)
% Calculate Wilson equation K-factors with scaling factor

ncomp = length(Pc);
lnK = zeros(ncomp, 1);

for i = 1:ncomp
    lnK(i) = factor * (5.373 * (1 + acentric(i)) * (1 - Tc(i)/temp) + log(Pc(i)/press));
end

% Limit to avoid overflow
lnK = max(lnK, -20);
lnK = min(lnK, 20);

end

%% TANGENT PLANE DISTANCE CALCULATION
function TPD = calculate_tpd(x_trial, z, fugcoef_trial, fugcoef_z)
% Calculate tangent plane distance

ncomp = length(z);
TPD = 0;

for i = 1:ncomp
    if x_trial(i) > 1e-15 && z(i) > 1e-15
        TPD = TPD + x_trial(i) * (log(x_trial(i)) + log(fugcoef_trial(i)) - ...
                                  log(z(i)) - log(fugcoef_z(i)));
    end
end

end

%% ESTIMATE COMPLEMENTARY PHASE
function x_comp = estimate_complementary_phase(x_trial, z)
% Estimate a complementary phase composition using lever rule approximation

ncomp = length(z);
x_comp = zeros(ncomp, 1);

% Simple complementary estimate
alpha = 0.5;  % Phase fraction guess
for i = 1:ncomp
    x_comp(i) = (z(i) - alpha * x_trial(i)) / (1 - alpha);
    x_comp(i) = max(x_comp(i), 1e-15);  % Ensure positive
end

% Normalize
x_comp = x_comp / sum(x_comp);

end