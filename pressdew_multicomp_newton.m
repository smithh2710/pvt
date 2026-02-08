%% CALCULATE DEW POINT PRESSURE BY NEWTON-RAPHSON METHOD
% -------------------------------------------------------------------------
% The Definition of Variables.
% comp_vap   : vapor composition
% pressd_ini : initial guess for dew point pressure
% temp       : temperature
% pressc     : critical pressure
% tempc      : critical temperature
% acentric   : acentric factor
% BIP        : binary interaction parameter
% tol        : tolerance for convergence
% maxiter    : maximum iteration
% eos_type   : 'PR' (default) or 'SRK'
% pressd     : dew point pressure
% comp_liq   : liquid composition at dew point
% -------------------------------------------------------------------------

function [pressd, comp_liq] = pressdew_multicomp_newton(comp_vap, pressd_ini, temp, pressc, tempc, acentric, BIP, tol, maxiter, eos_type)

% Default to PR if eos_type not specified
if nargin < 10 || isempty(eos_type)
    eos_type = 'PR';
end

ncomp = size(comp_vap,1);

% Input initial values.
pressd = pressd_ini;
K = wilsoneq(pressd, temp, pressc, tempc, acentric);
K = updatek(K, comp_vap, pressd, temp, pressc, tempc, acentric, BIP, eos_type);

fun = @(x) objfun(x(1:ncomp), comp_vap, x(ncomp + 1), temp, pressc, tempc, acentric, BIP, eos_type);

% model parameters.
m = [K; pressd];

for loop = 1:maxiter
    
    f = fun(m);
    J = jacobfun(f, m, fun);
    dm = -J\f;
    
    % Update m.
    m = m + dm;
    
    % Check convergence.
    eps = max(abs(f));
    if eps < tol
        break;
    end

end

% Echo a message if the iteration did not converge.
if loop >= maxiter
    fprintf('The iteration in pressdew_multicomp_newton() did not converge: eps = %e\n', eps);
else
    fprintf('iter = %d, objfun = [ ', loop);
    for i = 1:ncomp+1
        fprintf('%1.3e ', f(i));
    end
    fprintf(']\n');
end

K = m(1:ncomp);
pressd = m(ncomp + 1);

% Calculate composition in liquid phase.
comp_liq = calccompliq(K, comp_vap);

end

%% Calculate liquid composition from K-values
function comp_liq = calccompliq(K, comp_vap)

ncomp = size(comp_vap, 1);

comp_liq = zeros(ncomp, 1);
for i = 1:ncomp
    comp_liq(i) = comp_vap(i)/K(i);
end

end

%% Update K by successive substitution
function Knew = updatek(K, comp_vap, pressd, temp, pressc, tempc, acentric, BIP, eos_type)

ncomp = size(comp_vap, 1);
Knew = zeros(ncomp, 1);

% Calculate composition in liquid phase.
comp_liq = calccompliq(K, comp_vap);

% Calculate fugacity coefficients in vapor and liquid phase.
[fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, pressd, temp, pressc, tempc, acentric, BIP, eos_type);
[fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, pressd, temp, pressc, tempc, acentric, BIP, eos_type);

for i = 1:ncomp
   Knew(i) = fugcoef_liq(i)/fugcoef_vap(i); 
end

end

%% Objective function
% f_i = K_i - phi_i^L / phi_i^V,  i = 1,...,Nc
% f_{Nc+1} = sum_i(y_i / K_i) - 1
function f = objfun(K, comp_vap, pressd, temp, pressc, tempc, acentric, BIP, eos_type)

ncomp = size(K,1);

% Calculate composition in liquid phase.
comp_liq = calccompliq(K, comp_vap);

% Calculate fugacity coefficients in vapor and liquid phase.
[fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, pressd, temp, pressc, tempc, acentric, BIP, eos_type);
[fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, pressd, temp, pressc, tempc, acentric, BIP, eos_type);

% Calculate objective function.
f = zeros(ncomp + 1, 1);

for i = 1:ncomp
    f(i) = K(i) - fugcoef_liq(i)/fugcoef_vap(i);
end
f(ncomp + 1) = sum(comp_liq) - 1;

end

%% Jacobian matrix
% J_{ij} = df_i / dx_j
function J = jacobfun(f0, x, fun)
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