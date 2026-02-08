%% CALCULATE BUBBLE POINT PRESSURE BY SUCCESSIVE SUBSTITUTION METHOD
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

function [pressb, comp_vap] = pressbub_multicomp_ss(comp_liq, pressb_ini, temp, pressc, tempc, acentric, BIP, tol, maxiter, eos_type)

% Default to PR if eos_type not specified
if nargin < 10 || isempty(eos_type)
    eos_type = 'PR';
end

ncomp = size(comp_liq,1);

% Input initial values.
pressb = pressb_ini;
K = wilsoneq(pressb, temp, pressc, tempc, acentric);

for loop = 1:maxiter
    
    [f, Knew, pressbnew] = objfun(K, comp_liq, pressb, temp, pressc, tempc, acentric, BIP, eos_type);
    
    % Update K and bubble point pressure.
    K = Knew;
    pressb = pressbnew;
    
    % Check convergence.
    eps = abs(f);
    if eps < tol
        break;
    end
    
end

% Echo a message if the loop did not converge.
if loop >= maxiter
    fprintf('The iteration in pressbub_multicomp_ss() did not converge. eps = %E\n', eps);
else
    fprintf('Iteration = %d\n', loop);
end

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

%% Objective function
function [f, Knew, pressbnew] = objfun(K, comp_liq, pressb, temp, pressc, tempc, acentric, BIP, eos_type)

ncomp = size(K,1);

% Calculate vapor composition.
comp_vap = calccompvap(K, comp_liq);

% Calculate fugacity coefficients in vapor and liquid phase.
[fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, pressb, temp, pressc, tempc, acentric, BIP, eos_type);
[fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, pressb, temp, pressc, tempc, acentric, BIP, eos_type);

% Calculate new K and bubble point pressure.
fug_vap = zeros(ncomp, 1);
fug_liq = zeros(ncomp, 1);
Knew = zeros(ncomp,1);
pressbnew = 0;
for i = 1:ncomp
    
    fug_vap(i) = comp_vap(i)*fugcoef_vap(i)*pressb;
    fug_liq(i) = comp_liq(i)*fugcoef_liq(i)*pressb;
    % Calculate new K by successive substitution method.
    Knew(i) = fugcoef_liq(i)/fugcoef_vap(i);
    % Calculate new bubble point pressure.
    pressbnew = pressbnew + fug_liq(i)/fugcoef_vap(i);
    
end

% Calculate objective function f.
f = comp_liq'*Knew - 1;

end