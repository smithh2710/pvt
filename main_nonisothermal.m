function [comp_h, press_h, temp_h, pressbub_h, pressdew_h] = main_nonisothermal(h, h_ref, comp_ref, press_ref, temp_ref, ...
    pressc, tempc, acentric, BIP, M_gmol, vt_method, vt_params, thermal_params)

R = 8.3144598;
g = 9.80665;
M = M_gmol / 1000;
M_g = M_gmol;

tol = 1e-10;
maxiter = 1500;
n = length(comp_ref);

comp_ref = comp_ref(:);
M = M(:);
M_g = M_g(:);
pressc = pressc(:);
tempc = tempc(:);
acentric = acentric(:);

if nargin < 12 || isempty(vt_params)
    vt_params = struct();
end
if nargin < 11 || isempty(vt_method)
    vt_method = 0;
end
if nargin < 13 || isempty(thermal_params)
    thermal_params = struct();
end

[dT_dh, Cp_coeffs, H_ig_ref_opt] = parse_thermal_params(thermal_params, n, M_gmol);

delta_h = h - h_ref;
temp_h = temp_ref + dT_dh * delta_h;
delta_T = temp_h - temp_ref;

[vt_type, vt_opts] = parse_vt_input(vt_method, vt_params, n, temp_ref, pressc, tempc, acentric, M_gmol);

c_ref = calc_volume_translation(comp_ref, press_ref, temp_ref, pressc, tempc, acentric, M_gmol, vt_type, vt_opts);

[fugcoef_ref_eos, ~] = fugacitycoef_multicomp(comp_ref, press_ref, temp_ref, pressc, tempc, acentric, BIP);
if vt_type > 0 || vt_type == -1
    fugcoef_ref = fugcoef_ref_eos .* exp(-c_ref * press_ref / (R * temp_ref));
else
    fugcoef_ref = fugcoef_ref_eos;
end

f_ref = fugcoef_ref .* comp_ref * press_ref;

[~, H_mix_ref, ~, ~, H_abs_spec_ref] = calculate_absolute_enthalpy(...
    temp_ref, press_ref, comp_ref, pressc, tempc, acentric, BIP, M_g, Cp_coeffs, H_ig_ref_opt);

M_avg_ref = sum(comp_ref .* M_g);
H_over_M_ref = H_mix_ref / M_avg_ref;

f_h = zeros(n, 1);
for i = 1:n
    gravity_term = (M(i) * g * delta_h) / (R * temp_h);
    thermal_driving_force = H_over_M_ref - H_abs_spec_ref(i);
    thermal_term = (M_g(i) * thermal_driving_force * delta_T) / (R * temp_h * temp_ref);
    
    ln_f_ref = log(f_ref(i));
    ln_f_h = (temp_ref / temp_h) * ln_f_ref + gravity_term - thermal_term;
    f_h(i) = exp(ln_f_h);
end

initial_guess = [comp_ref; press_ref];

vt_opts_h = vt_opts;
vt_opts_h.temp = temp_h;

fun = @(x) residual_fugacity_noniso(x(1:n), x(end), f_h, temp_h, pressc, tempc, acentric, BIP, M_gmol, vt_type, vt_opts_h, R);

options = optimoptions('fsolve', 'Display', 'none', 'FunctionTolerance', 1e-15, 'StepTolerance', 1e-15, 'MaxIterations', maxiter);
solution = fsolve(fun, initial_guess, options);

comp_h = solution(1:n);
press_h = solution(end);

comp_h = max(comp_h, 1e-15);
comp_h = comp_h / sum(comp_h);

try 
    pressbub_ini_h = pressbubest_multicomp(comp_h, temp_h, pressc, tempc, acentric);
    if isnan(pressbub_ini_h) || pressbub_ini_h <= 0
        pressbub_ini_h = press_h * 0.9;
    end
    [pressbub_h, ~] = pressbub_multicomp_newton(comp_h, pressbub_ini_h, temp_h, pressc, tempc, acentric, BIP, tol, maxiter);
    if ~isreal(pressbub_h) || pressbub_h <= 0 || ~isfinite(pressbub_h)
       pressbub_h = NaN;
    end
catch 
   pressbub_h = NaN; 
end   

try 
    pressdew_ini_h = pressdewest_multicomp(comp_h, temp_h, pressc, tempc, acentric);
    if isnan(pressdew_ini_h) || pressdew_ini_h <= 0
        pressdew_ini_h = press_h * 0.5;
    end
    [pressdew_h, ~] = pressdew_multicomp_newton(comp_h, pressdew_ini_h, temp_h, pressc, tempc, acentric, BIP, tol, maxiter);
    if ~isreal(pressdew_h) || pressdew_h <= 0 || ~isfinite(pressdew_h)
       pressdew_h = NaN;
    end
catch 
   pressdew_h = NaN; 
end 

end


function [dT_dh, Cp_coeffs, H_ig_ref_opt] = parse_thermal_params(thermal_params, n, M_gmol)

if isfield(thermal_params, 'dT_dh')
    dT_dh = thermal_params.dT_dh;
else
    dT_dh = 0.025;
end

if isfield(thermal_params, 'Cp_coeffs')
    Cp_coeffs = thermal_params.Cp_coeffs;
else
    Cp_coeffs = estimate_Cp_coeffs(n, M_gmol);
end

if isfield(thermal_params, 'H_ig_ref')
    H_ig_ref_opt = thermal_params.H_ig_ref;
else
    H_ig_ref_opt = [];
end

end


function Cp_coeffs = estimate_Cp_coeffs(n, M_gmol)

R = 8.314462618;
Cp_coeffs = zeros(n, 4);

for i = 1:n
    MW = M_gmol(i);
    A = (-0.33886 + 0.02827*MW - 1.6982e-5*MW^2) * R;
    B = (3.14e-3 + 2.697e-4*MW) * R;
    C = (-1.5e-6 - 7.97e-8*MW) * R;
    D = 2.0e-10 * R;
    Cp_coeffs(i, :) = [A, B, C, D];
end

end


function F = residual_fugacity_noniso(comp_h, press_h, f_h, temp_h, pressc, tempc, acentric, BIP, M_gmol, vt_type, vt_opts, R)

n = length(f_h);

[fugcoef_h_eos, ~] = fugacitycoef_multicomp(comp_h, press_h, temp_h, pressc, tempc, acentric, BIP);

if vt_type > 0 || vt_type == -1
    c_shift = calc_volume_translation(comp_h, press_h, temp_h, pressc, tempc, acentric, M_gmol, vt_type, vt_opts);
    fugcoef_h = fugcoef_h_eos .* exp(-c_shift * press_h / (R * temp_h));
else
    fugcoef_h = fugcoef_h_eos;
end

F = zeros(n+1, 1);
for i = 1:n
    f_i_calc = comp_h(i) * press_h * fugcoef_h(i);
    F(i) = f_i_calc - f_h(i);
end

F(n+1) = sum(comp_h) - 1;

end


function c = calc_volume_translation(comp, press, temp, Pc, Tc, acentric, M_gmol, vt_type, vt_opts)

n = length(comp);

switch vt_type
    case 0
        c = zeros(n, 1);
    case -1
        c = vt_opts.c_direct;
    case 1
        [~, c] = peneloux_volume_shift(Pc, Tc, acentric);
    case 2
        [~, c] = magoulas_tassios_volume_shift(temp, Pc, Tc, acentric);
    case 3
        [~, c] = ungerer_batut_volume_shift(temp, Pc, Tc, acentric, M_gmol);
    case 4
        components = vt_opts.components;
        [~, c] = baled_volume_shift(temp, Pc, Tc, acentric, M_gmol, components);
    case 5
        comp = comp(:);
        Vc = vt_opts.Vc;
        components = vt_opts.components;
        BIP_local = zeros(n);
        [~, c] = abudour_volume_shift(comp, press, temp, Pc, Tc, acentric, Vc, components, BIP_local, false);
    case 6
        c = vt_opts.c_direct;
    otherwise
        c = zeros(n, 1);
end

c = c(:);

end


function [vt_type, vt_opts] = parse_vt_input(vt_method, vt_params, n, temp, Pc, Tc, acentric, M_gmol)

vt_opts = struct();
R = 8.3144598;

if isfield(vt_params, 'Vc')
    vt_opts.Vc = vt_params.Vc(:);
else
    Zc_est = 0.2905 - 0.085 * acentric;
    vt_opts.Vc = Zc_est .* R .* Tc ./ Pc * 1e6;
end

if isfield(vt_params, 'components')
    vt_opts.components = vt_params.components;
else
    vt_opts.components = cell(n, 1);
end

vt_opts.Pc = Pc;
vt_opts.Tc = Tc;
vt_opts.acentric = acentric;
vt_opts.M_gmol = M_gmol;
vt_opts.temp = temp;
vt_opts.n = n;

if isempty(vt_method)
    vt_type = 0;
    return;
end

if isnumeric(vt_method) && length(vt_method) == n
    vt_type = -1;
    vt_opts.c_direct = vt_method(:) * 1e-6;
    return;
end

if ischar(vt_method) || isstring(vt_method)
    vt_type = get_vt_method_number(vt_method);
else
    vt_type = vt_method;
end

if vt_type == 6
    if isfield(vt_params, 'c_custom') && ~isempty(vt_params.c_custom)
        c_custom = vt_params.c_custom(:);
        vt_opts.c_direct = c_custom * 1e-6;
    else
        error('For vt_method = 6, provide vt_params.c_custom [cm3/mol]');
    end
end

if vt_type < 0 || vt_type > 6
    vt_type = 0;
end

end


function num = get_vt_method_number(name)

name = lower(char(name));

switch name
    case {'none', 'no', '0', ''}
        num = 0;
    case {'peneloux', 'pen', '1'}
        num = 1;
    case {'magoulas_tassios', 'magoulas', 'tassios', 'mt', '2'}
        num = 2;
    case {'ungerer_batut', 'ungerer', 'batut', 'ub', '3'}
        num = 3;
    case {'baled', 'bal', '4'}
        num = 4;
    case {'abudour', 'abu', '5'}
        num = 5;
    case {'custom', 'user', '6'}
        num = 6;
    otherwise
        num = 0;
end

end