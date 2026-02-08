%% CALCULATE BUBBLE POINT PRESSURE BY NEWTON-RAPHSON METHOD
% -------------------------------------------------------------------------
% The Definition of Variables.
% comp_liq   : liquid composition
% pressb_ini : initial guess for bubble point pressure
% temp       : temperature
% pressc     : critical pressure
% tempc      : critical temperature
% acentric   : acentric factor
% BIP        : binary interaction parameter
% tol        : tolerance for convergence
% maxiter    : maximum iteration
% eos_type   : 'PR' (default) or 'SRK'
% pressb     : bubble point pressure
% comp_vap   : vapor composition at bubble point
% -------------------------------------------------------------------------

function [pressb, comp_vap] = pressbub_multicomp_newton(comp_liq, pressb_ini, temp, pressc, tempc, acentric, BIP, tol, maxiter, eos_type)

% Default to PR if eos_type not specified
if nargin < 10 || isempty(eos_type)
    eos_type = 'PR';
end

ncomp = size(comp_liq,1);

% Input initial values.
pressb = pressb_ini;
K = wilsoneq(pressb, temp, pressc, tempc, acentric);
K = updatekss(K, comp_liq, pressb, temp, pressc, tempc, acentric, BIP, eos_type);

% m : model parameters
fun = @(m) objfun(m(1:ncomp, 1), comp_liq, m(ncomp + 1, 1), temp, pressc, tempc, acentric, BIP, eos_type);
updatek = @(m) updatekss(m(1:ncomp, 1), comp_liq, m(ncomp + 1, 1), temp, pressc, tempc, acentric, BIP, eos_type);
x = [K; pressb];

for loop = 1:maxiter
    
    % Calculate update direction by using Newton-Raphson method.
    f = fun(x);
    J = jacob(f, x, fun);
    dx = -J\f;
    
    % Update x.
    x = x + dx;

    % Check convergence.
    eps = max(abs(f));
    if eps < tol
        break;
    end
    
end

% Echo a message if the loop did not converge.
if loop >= maxiter
    fprintf('The iteration in pressbub_multicomp_newton() did not converge: eps = %e\n', eps);
else
    fprintf('iter = %d, objfun = [ ', loop);
    for i = 1:ncomp+1
        fprintf('%1.3e ', f(i));
    end
    fprintf(']\n');
end

% Update K and pressb.
K = x(1:ncomp, 1);
pressb = x(ncomp + 1, 1);

% Calculate vapor composition.
comp_vap = calccompvap(K, comp_liq);

end

%% Calculate vapor composition from K-values
function comp_vap = calccompvap(K, comp_liq)

ncomp = size(comp_liq, 1);

comp_vap = zeros(ncomp, 1);
for i = 1:ncomp
    comp_vap(i) = K(i)*comp_liq(i);
end

end

%% Update K by using successive substitution
function K = updatekss(K, comp_liq, pressb, temp, pressc, tempc, acentric, BIP, eos_type)

ncomp = size(K,1);

% Calculate vapor composition.
comp_vap = calccompvap(K, comp_liq);

% Calculate fugacity coefficients in vapor and liquid phase.
[fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, pressb, temp, pressc, tempc, acentric, BIP, eos_type);
[fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, pressb, temp, pressc, tempc, acentric, BIP, eos_type);

% Calculate equilibrium constant.
K = zeros(ncomp, 1);
for i = 1:ncomp
    K(i) = fugcoef_liq(i)/fugcoef_vap(i);
end

end

%% Objective function
% f_i = K_i - phi_i^L / phi_i^V,  i = 1,...,Nc
% f_{Nc+1} = sum_i(z_i * K_i) - 1
function f = objfun(K, comp_liq, pressb, temp, pressc, tempc, acentric, BIP, eos_type)

ncomp = size(K,1);

% Calculate vapor composition.
comp_vap = calccompvap(K, comp_liq);

% Calculate fugacity coefficients in vapor and liquid phase.
[fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, pressb, temp, pressc, tempc, acentric, BIP, eos_type);
[fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, pressb, temp, pressc, tempc, acentric, BIP, eos_type);

% Calculate the objective function.
f = zeros(ncomp + 1, 1);
for i = 1:ncomp
    f(i) = K(i) - fugcoef_liq(i)/fugcoef_vap(i);
end
f(ncomp + 1) = sum(comp_vap) - 1;

end

%% Jacobian matrix
% J_{ij} = df_i / dx_j
function J = jacob(f0, x, fun)
N = size(x, 1);
J = zeros(N, N);
perturb_x = 1e-6;
for i = 1:N
    dx = zeros(N, 1);
    dx(i) = perturb_x;
    f1 = fun(x + dx);
    J(:, i) = (f1 - f0)/perturb_x;
end
end