%% Calculate phase mole fraction by solving Rachford-Rice equation
%
% phasefrac    : Phase fraction
% comp         : Phase composition
% converged    : Boolean flag indicating if Newton solver converged
% K            : K values or equilibrium ratios
% comp_overall : Overall composition
% tol          : Iteration tolerance
% maxiter      : Maximum number of iterations
function [phasefrac, comp, converged] = phasefraction(K, comp_overall, tol, maxiter)

% Calculate the initial estimate of phase mole fraction.
phasefrac = phasefracest(K, comp_overall);

if isempty(phasefrac)
    error('phasefraction:noFeasibleRegion', ...
          'No feasible region found for Rachford-Rice equation. Check K-values and composition.');
end

fun = @(x) minfun(K, comp_overall, x);
grad = @(x) gradfun(K, comp_overall, x);
hessian = @(x) hessianfun(K, comp_overall, x);
maxstepsize = @(x, dx) maxstepsizefun(K, comp_overall, x, dx);

[phasefrac, converged] = newton(fun, grad, hessian, phasefrac, tol, maxiter, maxstepsize);

if ~converged
    warning('phasefraction:notConverged', ...
            'Newton solver did not converge within %d iterations.', maxiter);
end

comp = calccomp(K, comp_overall, phasefrac);

end

%% Calculate variables.

function t = calct(K, phasefrac)
% K   : equilibrium constant
% phasefrac: phase mole fraction
nphase = size(phasefrac,1);  % the number of phases.
ncomp = size(K,1);         % the number of components.
Ktemp = ones(ncomp, nphase) - K;
t = ones(ncomp,1) - Ktemp*phasefrac;

end

function comp = calccomp(K, comp_overall, phasefrac)

ncomp = size(K, 1);
nphase = size(phasefrac, 1) + 1;
t = calct(K, phasefrac);

comp = zeros(ncomp, nphase);
for i = 1:ncomp
    comp(i, nphase) = comp_overall(i)/t(i);
    for j = 1:nphase-1
        comp(i, j) = K(i, j)*comp(i, nphase);
    end
end

end

%% Objective function to be minimized
% Based on Okuno et al., (2010)

function f = minfun(K, comp_overall, phasefrac)
ncomp = size(comp_overall, 1);
t = calct(K, phasefrac);
f = 0;
for i = 1:ncomp
    if t(i) <= 0
        f = inf;  % signal infeasibility to solver
        return;
    end
    f = f - comp_overall(i)*log(t(i));
end
end

function g = gradfun(K, comp_overall, phasefrac)

nphase = size(phasefrac, 1); % the number of phases - 1.
ncomp = size(K,1);              % the number of components.
t = calct(K, phasefrac);
temp = zeros(ncomp, 1);
for i = 1:ncomp
    if abs(t(i)) < 1e-15
        temp(i) = sign(t(i)) * comp_overall(i) / 1e-15;
    else
        temp(i) = comp_overall(i)/t(i);
    end
end
g = (ones(ncomp, nphase) - K)'*temp;
end

function H = hessianfun(K, z, beta)

nphase = size(beta,1);  % the number of phase - 1.
ncomp = size(K,1);      % the number of components
t = calct(K,beta);

% Calculate Hessian matrix.
H = zeros(nphase, nphase);
for j = 1:nphase
    for k = 1:nphase
        for i = 1:ncomp
            ti_sq = t(i)^2;
            if ti_sq < 1e-30
                ti_sq = 1e-30;
            end
            H(j, k) = H(j, k) + (1 - K(i, j))*(1 - K(i, k))*z(i)/ti_sq;
        end
    end
end

end

%% Phase mole fraction Calculation

function phasefrac = phasefracest(K, comp_overall)

% Calculate a and b.
a = calca(K);
b = calcb(K, comp_overall);

% Calculate feasible region.
fr = feasibleregion(a, b);

% center of min feasible region
nphase = size(fr,1);
ncomp = size(fr,2);
phasefrac = zeros(nphase,1);
for i = 1:nphase
    for j = 1:ncomp
        phasefrac(i) = phasefrac(i) + fr(i,j);
    end
end
phasefrac = phasefrac/ncomp;

end

%% Calculate a feasible region for the minimizing function
% Based on Okuno et al. (2010)

function a = calca(K)

ncomp = size(K,1);      % the number of components.
nphase = size(K,2);     % (the number of phases) - 1.

a = ones(ncomp, nphase) - K;

end

function b = calcb(K, comp_overall)

nphase = size(K,2); % the number of phases - 1.
ncomp  = size(K,1); % the number of components.

b = zeros(ncomp, 1);

for i = 1:ncomp
    x1 = 1 - comp_overall(i);
    x2 = 1 - K(i,1)*comp_overall(i);
    for j = 2:nphase
        temp = 1 - K(i,j)*comp_overall(i);
        if (x2 > temp)
            x2 = temp;
        end
    end
    if x1 < x2
        b(i) = x1;
    else
        b(i) = x2;
    end
end

end

function fr = feasibleregion(a, b)

ncomp = size(b, 1);
nphase = size(a, 2);
index = nchoosek(1:ncomp, nphase);
fr = [];
for i = 1:size(index,1)

    at = [];
    bt = [];

    for j = 1:nphase
        at = cat(1, at, a(index(i,j),:));
        bt = cat(1, bt, b(index(i,j)));
    end

    phasefrac = at\bt;
    abeta = a*phasefrac;
    flag = true;

    for j = 1:ncomp
        if abeta(j) > b(j)
            is_vertex = any(index(i,:) == j);
            if ~is_vertex
                flag = false;
                break;
            end
        end
    end

    if flag
        fr = cat(2, fr, phasefrac);
    end

end

end

%% Calculate the maximum step size for Newton iteration

function maxstepsize = maxstepsizefun(K, comp_overall, phasefrac, d)

ncomp = size(K, 1);         % the number of components

a = calca(K);
b = calcb(K, comp_overall);

temp1 = b - a*phasefrac;
temp2 = a*d;

maxstepsize = 1;
for i = 1:ncomp
    if temp2(i) > 0
        boundary = temp1(i)/temp2(i);
        if boundary < maxstepsize
            maxstepsize = boundary;
        end
    end
end
maxstepsize = max(maxstepsize, 1e-10);  % avoid zero step

end
