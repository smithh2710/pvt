%% VAPOR-LIQUID EQUILIBRIUM FLASH CALCULATION
% -------------------------------------------------------------------------
% The Definition of Variables.
% press        : Pressure
% % temp         : Temperature
% % comp_overall : Overall composition
% % pressc       : Critical pressure
% % tempc        : Critical temperature
% % acentric     : Acentric factor
% % BIP          : Binary interaction parameters
% % tol          : Tolerance for convergence
% % maxiter      : Maximum iterations
% % eos_type     : 'PR' (default) or 'SRK'
% % K            : Equilibrium K-values
% % comp_vap     : Vapor phase composition
% % comp_liq     : Liquid phase composition
% % phasefrac    : Vapor phase mole fraction
% % -------------------------------------------------------------------------
% 
% function [K, comp_vap, comp_liq, phasefrac] = vaporliquideq(press, temp, comp_overall, pressc, tempc, acentric, BIP, tol, maxiter, eos_type)
% 
% % Default to PR if eos_type not specified
% if nargin < 10 || isempty(eos_type)
%     eos_type = 'PR';
% end
% 
% % Initial estimate of equilibrium constant
% K = wilsoneq(press, temp, pressc, tempc, acentric);
% 
% fun = @(x) objfun(x, comp_overall, press, temp, pressc, tempc, acentric, BIP, eos_type);
% grad = @(x) jacobfun(x, comp_overall, press, temp, pressc, tempc, acentric, BIP, eos_type);
% 
% for i = 1:maxiter
% 
%     f = fun(K);
% 
%     eps = max(abs(f));
%     if eps < tol
%         break;
%     end
% 
%     dfdK = grad(K);
%     dK = -dfdK\f;
%     K = K + dK;
% 
% end
% 
% [phasefrac, comp] = phasefraction(K, comp_overall, tol, maxiter);
% comp_vap = comp(:, 1);
% comp_liq = comp(:, 2);
% 
% end
% 
% %% Objective function
% % f_i = K_i - phi_i^L / phi_i^V
% function f = objfun(K, comp_overall, press, temp, pressc, tempc, acentric, BIP, eos_type)
% 
% ncomp = size(comp_overall, 1);
% tol = 1e-5;
% maxiter = 20;
% [~, comp] = phasefraction(K, comp_overall, tol, maxiter);
% comp_vap = comp(:, 1);
% comp_liq = comp(:, 2);
% 
% [fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, press, temp, pressc, tempc, acentric, BIP, eos_type);
% [fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, press, temp, pressc, tempc, acentric, BIP, eos_type);
% 
% f = zeros(ncomp, 1);
% for i = 1:ncomp
%     f(i) = K(i) - fugcoef_liq(i)/fugcoef_vap(i);
% end
% 
% end
% 
% %% Jacobian matrix (numerical)
% function J = jacobfun(K, comp_overall, press, temp, pressc, tempc, acentric, BIP, eos_type)
% 
% ncomp = size(comp_overall, 1);
% 
% fun = @(x) objfun(x, comp_overall, press, temp, pressc, tempc, acentric, BIP, eos_type);
% 
% f0 = fun(K);
% 
% perturb_K = 1e-6;
% J = zeros(ncomp, ncomp);
% 
% for i = 1:ncomp
% 
%     dK = zeros(ncomp, 1);
%     dK(i) = perturb_K;
%     K1 = K + dK;
% 
%     f1 = fun(K1);
%     J(:, i) = (f1 - f0)/perturb_K;
% 
% end
% 
% end




 %% 



function [K, comp_vap, comp_liq, phasefrac] = vaporliquideq(press, temp, comp_overall, pressc, tempc, acentric, BIP, tol, maxiter, eos_type)

if nargin < 10 || isempty(eos_type)
    eos_type = 'PR';
end

K = wilsoneq(press, temp, pressc, tempc, acentric);
K = max(K, 1e-12);

ncomp = length(comp_overall);
eps_prev = Inf;
stall_count = 0;

for i = 1:maxiter

    f = objfun(K, comp_overall, press, temp, pressc, tempc, acentric, BIP, eos_type);

    if any(~isfinite(f))
        break;
    end

    eps_cur = max(abs(f));
    if eps_cur < tol
        break;
    end

    if eps_cur > 1e10
        break;
    end

    if eps_cur >= eps_prev
        stall_count = stall_count + 1;
        if stall_count > 5
            break;
        end
    else
        stall_count = 0;
    end
    eps_prev = eps_cur;

    dfdK = jacobfun(K, comp_overall, press, temp, pressc, tempc, acentric, BIP, eos_type);

    if rcond(dfdK) < 1e-15
        break;
    end

    dK = -dfdK\f;

    if any(~isfinite(dK))
        break;
    end

    max_step = 0.5 * max(abs(K));
    scale = max_step / max(max(abs(dK)), max_step);
    dK = dK * scale;

    K_new = K + dK;
    K_new = max(K_new, 1e-12);
    K = K_new;

end

[phasefrac, comp] = phasefraction(K, comp_overall, tol, maxiter);
comp_vap = comp(:, 1);
comp_liq = comp(:, 2);

end


function f = objfun(K, comp_overall, press, temp, pressc, tempc, acentric, BIP, eos_type)

ncomp = size(comp_overall, 1);
tol = 1e-5;
maxiter = 20;
[~, comp] = phasefraction(K, comp_overall, tol, maxiter);
comp_vap = comp(:, 1);
comp_liq = comp(:, 2);

[fugcoef_liq, ~] = fugacitycoef_multicomp_liquid(comp_liq, press, temp, pressc, tempc, acentric, BIP, eos_type);
[fugcoef_vap, ~] = fugacitycoef_multicomp_vapor(comp_vap, press, temp, pressc, tempc, acentric, BIP, eos_type);

f = zeros(ncomp, 1);
for i = 1:ncomp
    if fugcoef_vap(i) > 0 && isfinite(fugcoef_liq(i)) && isfinite(fugcoef_vap(i))
        f(i) = K(i) - fugcoef_liq(i)/fugcoef_vap(i);
    else
        f(i) = 0;
    end
end

end


function J = jacobfun(K, comp_overall, press, temp, pressc, tempc, acentric, BIP, eos_type)

ncomp = size(comp_overall, 1);

fun = @(x) objfun(x, comp_overall, press, temp, pressc, tempc, acentric, BIP, eos_type);

f0 = fun(K);

perturb_K = 1e-6;
J = zeros(ncomp, ncomp);

for i = 1:ncomp
    dK = zeros(ncomp, 1);
    dK(i) = perturb_K;
    K1 = K + dK;
    f1 = fun(K1);
    J(:, i) = (f1 - f0)/perturb_K;
end

end