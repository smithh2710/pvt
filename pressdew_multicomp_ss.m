%% CALCULATE DEW POINT PRESSURE BY SUCCESSIVE SUBSTITUTION METHOD
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

function [pressd, comp_liq] = pressdew_multicomp_ss(comp_vap, pressd_ini, temp, pressc, tempc, acentric, BIP, tol, maxiter, eos_type)

% Default to PR if eos_type not specified
if nargin < 10 || isempty(eos_type)
    eos_type = 'PR';
end

ncomp = size(comp_vap,1);

% Input initial values.
pressd = pressd_ini;
K = wilsoneq(pressd, temp, pressc, tempc, acentric);

for loop = 1:maxiter
    
    [f, Knew, pressdnew] = objfun(K, comp_vap, pressd, temp, pressc, tempc, acentric, BIP, eos_type);
    
    % Update K and dew point pressure.
    K = Knew;
    pressd = pressdnew;
    
    % Check convergence.
    eps = abs(f);
    if eps < tol
        break;
    end
    
end

% Echo a message if the iteration did not converge.
if loop >= maxiter
    fprintf('The iteration in pressdew_multicomp_ss() did not converge. eps = %E\n', eps);
else
    fprintf('Iteration = %d\n', loop);
end

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
function [f, Knew, pressdnew] = objfun(K, comp_vap, pressd, temp, pressc, tempc, acentric, BIP, eos_type)

ncomp = size(K,1);

% Calculate composition in liquid phase.
comp_liq = calccompliq(K, comp_vap);

% Calculate fugacity coefficients in vapor and liquid phase.
[fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, pressd, temp, pressc, tempc, acentric, BIP, eos_type);
[fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, pressd, temp, pressc, tempc, acentric, BIP, eos_type);

% Calculate new K and dew point pressure.
fug_vap = zeros(ncomp, 1);
fug_liq = zeros(ncomp, 1);
Knew = zeros(ncomp,1);
f = 0;
pressdnew = 0;

for i = 1:ncomp
    
    fug_vap(i) = comp_vap(i)*fugcoef_vap(i)*pressd;
    fug_liq(i) = comp_liq(i)*fugcoef_liq(i)*pressd;
    % Calculate new K.
    Knew(i) = fugcoef_liq(i)/fugcoef_vap(i);
    % Calculate new dew point pressure.
    pressdnew = pressdnew + fug_vap(i)/fugcoef_liq(i);
    % Calculate the objective function.
    f = f + comp_vap(i)/Knew(i);
    
end

f = f - 1;

end