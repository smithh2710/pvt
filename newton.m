%% Compute the minimum of a given function using Newton-Raphson iteration.
%
% fun        : Scalar function to be minimized
% grad       : Gradient function of 'fun'
% hessian    : Hessian function of 'fun'
% x0         : Initial vector
% tol        : Tolerance
% maxiter    : Maximum iteration
% maxstepfun : (optional) Function handle returning max step size alpha
%              given current x and Newton direction dx
function [x, converged] = newton(fun, grad, hessian, x0, tol, maxiter, maxstepfun)
  x = x0;
  step_norm = 1.0;
  iter = 0;
  while (step_norm > tol && iter < maxiter)
    dfdx = grad(x);
    H = hessian(x);
    dx = -H\dfdx;
    if nargin >= 7 && ~isempty(maxstepfun)
        alpha = maxstepfun(x, dx);
        alpha = min(alpha, 1.0) * 0.99;  % stay strictly inside feasible region
        dx = alpha * dx;
    end
    x = x + dx;
    step_norm = norm(dx);
    iter = iter + 1;
  end
  converged = (step_norm <= tol);
end
