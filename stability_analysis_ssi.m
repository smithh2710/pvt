function [stability, x_trial1, x_trial2, TPD_min, converged] = stability_analysis_ssi(comp, press, temp, Pc, Tc, acentric, BIP, tol, maxiter, eos_type)
% PHASE STABILITY ANALYSIS USING SUCCESSIVE SUBSTITUTION ITERATION (SSI)
% Based on Michelsen (1982a,b) method
%
% Inputs:
%   comp     : Overall composition (mole fractions)
%   press    : Pressure [Pa]
%   temp     : Temperature [K]
%   Pc       : Critical pressure vector [Pa]
%   Tc       : Critical temperature vector [K]
%   acentric : Acentric factor vector
%   BIP      : Binary interaction parameter matrix
%   tol      : Convergence tolerance (default: 1e-9)
%   maxiter  : Maximum iterations (default: 100)
%   eos_type : 'PR' (default) or 'SRK'
%
% Outputs:
%   stability : 1 = stable (single phase), 2 = unstable (will split)
%   x_trial1  : First trial phase composition (if unstable)
%   x_trial2  : Second trial phase composition (if unstable)
%   TPD_min   : Minimum tangent plane distance value
%   converged : Convergence flag

%% Initialize parameters
if nargin < 8 || isempty(tol)
    tol = 1e-9;
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

% Calculate reference fugacity coefficients for feed composition z
[fugcoef_z, ~] = fugacitycoef_multicomp(comp, press, temp, Pc, Tc, acentric, BIP, eos_type);

% Check for valid fugacity coefficients
if any(~isfinite(fugcoef_z)) || any(fugcoef_z <= 0)
    % Cannot perform stability analysis
    stability = 1;
    x_trial1 = comp;
    x_trial2 = comp;
    TPD_min = 0;
    converged = false;
    return;
end

lnz_plus_lnphi_z = log(comp) + log(fugcoef_z);

%% Generate initial K-value guesses
% Wilson K-values as base
K_wilson = wilsoneq(press, temp, Pc, Tc, acentric);
K_wilson = max(K_wilson, 1e-10);
K_wilson = min(K_wilson, 1e10);

% Build initial K-value sets: 4 Wilson-based + Nc pure component guesses
n_initial = ncomp + 4;
K_initial_sets = zeros(ncomp, n_initial);

% Wilson-based guesses
K_initial_sets(:, 1) = 1 ./ K_wilson;       % Liquid-like
K_initial_sets(:, 2) = K_wilson;            % Vapor-like
K_initial_sets(:, 3) = (K_wilson).^(1/3);   % Intermediate
K_initial_sets(:, 4) = (K_wilson).^(3);     % Extreme

% Pure component guesses
for j = 1:ncomp
    K_pure = ones(ncomp, 1) * 0.1 / ((ncomp - 1) * max(min(comp), 1e-10));
    K_pure(j) = 0.9 / max(comp(j), 1e-10);
    K_initial_sets(:, 4 + j) = K_pure;
end

%% Main SSI loop for all trial phases
n_trials = size(K_initial_sets, 2);
TPD_values = inf(n_trials, 1);
X_solutions = zeros(ncomp, n_trials);
converged_flags = false(n_trials, 1);

for trial = 1:n_trials
    K_initial = K_initial_sets(:, trial);
    
    % Initialize trial composition
    X = comp .* K_initial;
    X = max(X, 1e-15);
    
    % SSI iterations
    converged_trial = false;
    
    for iter = 1:maxiter
        X_old = X;
        
        % Normalize X to get mole fractions
        sum_X = sum(X);
        if sum_X < 1e-15
            break;  % Near-zero, skip this trial
        end
        x = X / sum_X;
        
        % Protect against zero/negative compositions
        x = max(x, 1e-15);
        x = x / sum(x);
        
        % Calculate fugacity coefficients for trial composition
        try
            [fugcoef_X, ~] = fugacitycoef_multicomp(x, press, temp, Pc, Tc, acentric, BIP, eos_type);
        catch
            break;  % Fugacity calculation failed, skip trial
        end
        
        % Check for valid fugacity coefficients
        if any(~isfinite(fugcoef_X)) || any(fugcoef_X <= 0)
            break;
        end
        
        % SSI update: X_i = exp(ln z_i + ln φ_i(z) - ln φ_i(x))
        for i = 1:ncomp
            X(i) = exp(lnz_plus_lnphi_z(i) - log(fugcoef_X(i)));
        end
        
        % Protect against overflow/underflow
        X = max(X, 1e-15);
        X = min(X, 1e10);
        
        % Check for trivial solution (X ≈ z)
        sum_X_new = sum(X);
        x_new = X / sum_X_new;
        if max(abs(x_new - comp)) < 1e-6
            break;  % Trivial solution
        end
        
        % Check for divergence
        if any(~isfinite(X))
            break;
        end
        
        % Check convergence
        error_norm = norm(X - X_old) / max(norm(X_old), 1e-10);
        if error_norm <= tol
            converged_trial = true;
            break;
        end
    end
    
    % Store results if converged
    if converged_trial
        sum_X = sum(X);
        
        % Calculate TPD at solution point
        % TPD = Σ X_i * (ln X_i + ln φ_i(x) - ln z_i - ln φ_i(z))
        % At stationary point: TPD = 1 - Σ X_i  (modified form)
        % Or equivalently: TPD = -ln(Σ X_i)
        if sum_X > 1e-15
            TPD = 1 - sum_X;  % Modified TPD (negative = unstable)
        else
            TPD = inf;
        end
        
        X_solutions(:, trial) = X;
        TPD_values(trial) = TPD;
        converged_flags(trial) = true;
        
        % Early termination if clear instability detected
        if TPD < -0.01
            break;
        end
    end
end

%% Determine stability from all trials
valid_indices = converged_flags & isfinite(TPD_values);

if any(valid_indices)
    valid_TPD = TPD_values(valid_indices);
    valid_X = X_solutions(:, valid_indices);
    
    [TPD_min, min_idx] = min(valid_TPD);
    
    % Normalize to get composition
    x_trial1 = valid_X(:, min_idx) / sum(valid_X(:, min_idx));
    
    % Stability criterion: TPD < 0 indicates instability
    if TPD_min < -1e-8
        stability = 2;  % Unstable (will split)
        
        % Find second trial phase
        if size(valid_X, 2) > 1
            temp_TPD = valid_TPD;
            temp_TPD(min_idx) = inf;
            [~, idx_second] = min(temp_TPD);
            
            x_second = valid_X(:, idx_second) / sum(valid_X(:, idx_second));
            if norm(x_second - x_trial1) > 0.05
                x_trial2 = x_second;
            else
                x_trial2 = estimate_complementary_phase(x_trial1, comp);
            end
        else
            x_trial2 = estimate_complementary_phase(x_trial1, comp);
        end
    else
        stability = 1;  % Stable (single phase)
        x_trial2 = comp;
    end
    
    converged = true;
else
    % No converged solutions - assume stable
    stability = 1;
    x_trial1 = comp;
    x_trial2 = comp;
    TPD_min = 0;
    converged = false;
end

end

%% Helper function to estimate complementary phase
function x_comp = estimate_complementary_phase(x_trial, z)
% Estimate complementary phase using material balance

ncomp = length(z);
x_comp = zeros(ncomp, 1);

% Assume 50% phase split
alpha = 0.5;

for i = 1:ncomp
    % Material balance: z_i = alpha * x_trial_i + (1-alpha) * x_comp_i
    x_comp(i) = (z(i) - alpha * x_trial(i)) / (1 - alpha);
    x_comp(i) = max(x_comp(i), 1e-15);
end

% Normalize
x_comp = x_comp / sum(x_comp);

end