%==========================================================================
% PFHX-5 pinch analysis: BASELINE + EJECTOR CASES, lowered T_HENE_in
%
% Reads ejector operating points from pinch_inputs.csv (one row per case)
% and runs pinch analysis for each. Baseline is unchanged, fully specified
% in run_baseline (T_HENE_in = 27.45 K, the Moro reference).
%
% For each ejector case the per-case T_HENE_in is read from the
% T_HENE_in_K column of pinch_inputs.csv. With the lowered T_HENE_in
% values, the baseline pinch target dT_min = 2.48 K is feasible at both
% ejector cases, so no fallback target is required.
%
% Output files:
%   pinch_curves_baseline.csv
%   pinch_curves_ejector_<label>.csv     one per row of pinch_inputs.csv
%
% Structure:
%   main script        : orchestrates baseline + N ejector cases
%   run_baseline()     : evaluates baseline (Riccardo / Moro) PFHX-5
%                        no solver; all streams fully specified.
%                        T_HENE_in fixed at 27.45 K (Moro).
%   run_ejector(spec, DT_TARGET, N_hot, N_cold, N_comp)
%                      : solves ejector case for (mdot_HENE, T_HENE_out)
%                        such that duty closure holds and dT_min = DT_TARGET.
%                        Stream conditions including T_HENE_in come from
%                        the `spec` struct (see run_ejector header).
%   run_ejector_with_fallback(spec, DT_primary, DT_fallback, ...)
%                      : retained for backward compatibility but not used
%                        in this run.
%
% Helpers (local functions): test_equilhyd, xp_equilibrium, build_HENE,
% pinch_residual, interp_clipped, write_curves_csv, compute_composite,
% pinch_metrics, build_hot
%==========================================================================

clear; clc; close all;

%% -------------------- Top-level parameters -----------------------------
% Each ejector case is run at the baseline pinch target dT_min = 2.48 K.
% This is feasible at both cases because the per-case T_HENE_in (read from
% pinch_inputs.csv) has been lowered from the Moro 27.45 K to a value
% just below the feasibility threshold for the 2.48 K target.
DT_TARGET          = 2.48;     % K, baseline pinch target (no fallback)
DT_FEASIBILITY_TOL = 0.01;     % K, tolerance on whether target is met

% Resolutions used in both baseline and ejector cases
N_hot  = 401;
N_cold = 401;
N_comp = 801;

%% -------------------- Run BASELINE -------------------------------------
fprintf('########################################################\n');
fprintf('# CASE 0: BASELINE (Riccardo / Moro)                    \n');
fprintf('########################################################\n\n');

base = run_baseline(N_hot, N_cold, N_comp);
write_curves_csv('pinch_curves_baseline.csv', ...
                 base.Q_hot_T, base.T_hot_grid, ...
                 base.Q_cold_comp, base.T_comp);
fprintf('Wrote pinch_curves_baseline.csv\n\n');

%% -------------------- Read ejector cases from CSV ----------------------
% Look for pinch_inputs.csv in this directory first; if not found, fall
% back to ../ejector/pinch_inputs.csv (the standard pipeline layout).
% If found in the sibling directory, copy it locally so the run is
% self-contained and reproducible.
inputs_path = 'pinch_inputs.csv';
if ~isfile(inputs_path)
    sibling_path = fullfile('..', 'ejector', 'pinch_inputs.csv');
    if isfile(sibling_path)
        fprintf(['pinch_inputs.csv not in current directory. ' ...
                 'Copying from %s ...\n'], sibling_path);
        copyfile(sibling_path, inputs_path);
    else
        error(['Could not find pinch_inputs.csv in either:\n' ...
               '  %s\n  %s\n' ...
               'Generate it first by running run_pinch_inputs.m in ' ...
               'the ejector directory.'], ...
              fullfile(pwd, inputs_path), ...
              fullfile(pwd, sibling_path));
    end
end
inputs_tbl = readtable(inputs_path);

required_cols = {'case_label','P_hot_bar','T_hot_in_K','T_hot_out_K', ...
                 'mdot_H_kgps','mdot_R_kgps','T_R_in_K','T_HENE_in_K', ...
                 'Q_PFHX5_kW'};
missing = setdiff(required_cols, inputs_tbl.Properties.VariableNames);
if ~isempty(missing)
    error('pinch_inputs.csv is missing required columns: %s', ...
          strjoin(missing, ', '));
end

n_cases = height(inputs_tbl);
fprintf('Loaded %d ejector case(s) from %s\n\n', n_cases, inputs_path);

%% -------------------- Run each ejector case ----------------------------
ejector_results = cell(n_cases, 1);

for k = 1:n_cases
    row = inputs_tbl(k, :);
    label = char(row.case_label);

    fprintf('########################################################\n');
    fprintf('# CASE %d/%d : EJECTOR  label = %s\n', k, n_cases, label);
    fprintf('#   P_hot     = %.3f bar\n',  row.P_hot_bar);
    fprintf('#   T_hot_in  = %.3f K\n',    row.T_hot_in_K);
    fprintf('#   T_hot_out = %.3f K\n',    row.T_hot_out_K);
    fprintf('#   mdot_H    = %.4f kg/s\n', row.mdot_H_kgps);
    fprintf('#   mdot_R    = %.4f kg/s\n', row.mdot_R_kgps);
    fprintf('#   T_R_in    = %.3f K\n',    row.T_R_in_K);
    fprintf('#   T_HENE_in = %.3f K\n',    row.T_HENE_in_K);
    fprintf('#   Q_PFHX5   = %.2f kW (ejector-side spec, for comparison)\n', ...
            row.Q_PFHX5_kW);
    fprintf('########################################################\n\n');

    spec = struct( ...
        'P_hot',        row.P_hot_bar, ...
        'T_hot_in',     row.T_hot_in_K, ...
        'T_hot_out',    row.T_hot_out_K, ...
        'mdot_H',       row.mdot_H_kgps, ...
        'mdot_R',       row.mdot_R_kgps, ...
        'T_R_in',       row.T_R_in_K, ...
        'T_HENE_in',    row.T_HENE_in_K, ...
        'Q_PFHX5_spec', row.Q_PFHX5_kW);

    fprintf('--- Solving at dT_min = %.2f K ---\n', DT_TARGET);
    ej_run = run_ejector(spec, DT_TARGET, N_hot, N_cold, N_comp);
    ej_run.target = DT_TARGET;
    ej_run.label  = label;

    % Sanity check that we actually hit the target
    if ej_run.dT_min < DT_TARGET - DT_FEASIBILITY_TOL
        warning(['Case %s: dT_min target = %.2f K NOT met ' ...
                 '(achieved %.3f K). Lower T_HENE_in further.'], ...
                label, DT_TARGET, ej_run.dT_min);
    else
        fprintf(['Target dT_min = %.2f K met (achieved %.3f K).\n\n'], ...
                DT_TARGET, ej_run.dT_min);
    end

    % Pack into the same shape as the old run_ejector_with_fallback so
    % downstream code (final summary, write_curves_csv) works unchanged.
    ej.primary               = ej_run;
    ej.primary_feasible      = (ej_run.dT_min >= DT_TARGET - DT_FEASIBILITY_TOL);
    ej.dT_min_max_achievable = ej_run.dT_min;
    ej.feasible              = ej_run;
    ej.label                 = label;
    ejector_results{k}       = ej;

    % Write the composite (target was met)
    out_csv = sprintf('pinch_curves_ejector_%s.csv', label);
    write_curves_csv(out_csv, ...
                     ej.feasible.Q_hot_T,    ej.feasible.T_hot_grid, ...
                     ej.feasible.Q_cold_comp, ej.feasible.T_comp);
    fprintf('Wrote %s\n\n', out_csv);
end

%% -------------------- Final summary ------------------------------------
fprintf('########################################################\n');
fprintf('# SUMMARY                                                \n');
fprintf('########################################################\n');

% Build column header
header = sprintf('  %-18s', 'quantity');
header = [header sprintf('  %12s', 'baseline')];
for k = 1:n_cases
    header = [header sprintf('  %12s', ejector_results{k}.label)];
end
fprintf('%s\n', header);

% Row helper
print_row = @(name, fmt, base_val, case_vals) ...
    fprintf(['  %-18s' repmat(['  ' fmt], 1, 1+numel(case_vals)) '\n'], ...
            name, base_val, case_vals);

% Numeric extraction helper
get_field_array = @(field) cellfun(@(s) s.feasible.(field), ejector_results);

print_row('Q_hot       [kW]',  '%12.2f', base.Q_hot_total,  get_field_array('Q_hot_total'));
print_row('Q_HENE      [kW]',  '%12.2f', base.Q_HENE_total, get_field_array('Q_HENE_total'));
print_row('Q_R         [kW]',  '%12.2f', base.Q_R_total,    get_field_array('Q_R_total'));
print_row('Q_cold      [kW]',  '%12.2f', base.Q_cold_total, get_field_array('Q_cold_total'));
print_row('Mismatch    [kW]',  '%+12.2f', ...
          base.Q_cold_total - base.Q_hot_total, ...
          get_field_array('Q_cold_total') - get_field_array('Q_hot_total'));
print_row('mdot_HENE   [kg/s]','%12.4f', base.mdot_HENE,    get_field_array('mdot_HENE'));
print_row('T_HENE_out  [K]',   '%12.3f', base.T_HENE_out,   get_field_array('T_HENE_out'));
print_row('dT_min      [K]',   '%12.3f', base.dT_min,       get_field_array('dT_min'));
print_row('Q_pinch     [kW]',  '%12.1f', base.Q_pinch,      get_field_array('Q_pinch'));
print_row('T_hot_p     [K]',   '%12.2f', base.T_hot_p,      get_field_array('T_hot_p'));
print_row('T_cold_p    [K]',   '%12.2f', base.T_cold_p,     get_field_array('T_cold_p'));

% Feasibility report (only meaningful for ejector cases; baseline N/A)
fprintf('\n  Feasibility against target dT_min = %.2f K:\n', DT_TARGET);
fprintf('  %-18s  %12s', 'case', '');
for k = 1:n_cases
    if ejector_results{k}.primary_feasible
        status = 'YES';
    else
        status = 'NO';
    end
    fprintf('  %12s', status);
end
fprintf('\n');
fprintf('  %-18s  %12s', 'max dT_min [K]', '');
for k = 1:n_cases
    fprintf('  %12.3f', ejector_results{k}.dT_min_max_achievable);
end
fprintf('\n');
fprintf('  %-18s  %12s', 'target [K]', '');
for k = 1:n_cases
    fprintf('  %12.3f', ejector_results{k}.feasible.target);
end
fprintf('\n');
fprintf('########################################################\n');


%% ======================================================================
%  WRAPPER: try primary target, fall back if infeasible
%  ======================================================================

function out = run_ejector_with_fallback(spec, DT_primary, DT_fallback, ...
                                          DT_tol, N_hot, N_cold, N_comp)
% RUN_EJECTOR_WITH_FALLBACK
%
%   Try DT_primary first. If the resulting dT_min is within DT_tol of
%   DT_primary, the primary target is feasible and we are done.
%   Otherwise, also run at DT_fallback to obtain a converged composite
%   for area sizing.
%
%   Returns:
%     out.primary             - struct returned by run_ejector at DT_primary
%     out.primary_feasible    - logical
%     out.dT_min_max_achievable - max dT_min reachable (= primary.dT_min,
%                                whether feasible or fminbnd-max)
%     out.fallback (optional) - struct returned by run_ejector at DT_fallback
%                               (only present if primary was infeasible)
%     out.feasible            - alias to whichever struct holds the
%                               "operating point" we'll use for area calc
%                               (primary if feasible, fallback otherwise)

    fprintf('--- Attempting primary target dT_min = %.2f K ---\n', DT_primary);
    primary = run_ejector(spec, DT_primary, N_hot, N_cold, N_comp);
    primary.target = DT_primary;

    % Did we actually meet the primary target?
    primary_feasible = (primary.dT_min >= DT_primary - DT_tol);

    out.primary               = primary;
    out.primary_feasible      = primary_feasible;
    out.dT_min_max_achievable = primary.dT_min;

    if primary_feasible
        fprintf(['Primary target dT_min = %.2f K is FEASIBLE ' ...
                 '(achieved %.3f K). No fallback needed.\n\n'], ...
                DT_primary, primary.dT_min);
        out.feasible = primary;
        return
    end

    fprintf(['Primary target dT_min = %.2f K is INFEASIBLE ' ...
             '(max achievable = %.3f K).\n'], ...
            DT_primary, primary.dT_min);
    fprintf('--- Re-running at fallback target dT_min = %.2f K ---\n', ...
            DT_fallback);
    fallback = run_ejector(spec, DT_fallback, N_hot, N_cold, N_comp);
    fallback.target = DT_fallback;

    out.fallback = fallback;
    out.feasible = fallback;
end


%% ======================================================================
%  CASE FUNCTIONS
%  ======================================================================

function out = run_baseline(N_hot, N_cold, N_comp)
% RUN_BASELINE  Evaluate PFHX-5 in baseline (Riccardo / Moro) configuration.
%
% All four stream conditions are fully specified. No solver. Computes
% Q_hot, Q_HENE, Q_R independently from REFPROP. Composite curves built
% as-is. dT_min reported on the common Q overlap (global pinch).

    %---- Stream specifications --------------------------------------------
    % Hot side: pressure drop ignored, P held constant at H6 value
    P_hot     = 75.0;             % bar
    T_hot_in  = 74.15;            % K (H6)
    T_hot_out = 27.95;            % K (H7)
    mdot_H    = 1.4847;           % kg/s

    % HENE (fully fixed)
    P_HENE    = 3.0;              % bar
    T_HENE_in = 27.45;            % K (HENE4)
    T_HENE_out= 74.05;            % K (HENE5)
    mdot_HENE = 6.9;              % kg/s
    x_He = 0.80; x_Ne = 0.20;

    % R stream: T_R_in raised above parahydrogen T_sat(1 bar) ~ 20.27 K
    P_R    = 1.0;                 % bar
    T_R_in = 20.40;               % K (H9, raised from 20.20)
    T_R_out= 46.65;               % K (H10)
    mdot_R = 0.4847;              % kg/s
    xp_R   = 0.999;

    %---- Hot-side enthalpy table ------------------------------------------
    use_equilhyd = test_equilhyd();
    [T_hot_grid, h_hot] = build_hot(T_hot_in, T_hot_out, P_hot, ...
                                    N_hot, use_equilhyd);
    Q_hot_T     = mdot_H * (h_hot - h_hot(1)) / 1000;
    Q_hot_total = Q_hot_T(end);
    fprintf('Q_hot  = %.2f kW\n', Q_hot_total);

    %---- HENE enthalpy table ----------------------------------------------
    T_HENE_grid = linspace(T_HENE_in, T_HENE_out, N_cold).';
    h_HENE = zeros(N_cold,1);
    for i = 1:N_cold
        h_HENE(i) = refpropm('H','T',T_HENE_grid(i),'P',P_HENE*100, ...
                             'helium','neon',[x_He, x_Ne]);
    end
    Q_HENE_T     = mdot_HENE * (h_HENE - h_HENE(1)) / 1000;
    Q_HENE_total = Q_HENE_T(end);
    fprintf('Q_HENE = %.2f kW\n', Q_HENE_total);

    %---- R stream enthalpy table ------------------------------------------
    T_R_grid = linspace(T_R_in, T_R_out, N_cold).';
    h_R = zeros(N_cold,1);
    for i = 1:N_cold
        h_R(i) = refpropm('H','T',T_R_grid(i),'P',P_R*100, ...
                          'orthohyd','parahyd',[1-xp_R, xp_R]);
    end
    Q_R_T     = mdot_R * (h_R - h_R(1)) / 1000;
    Q_R_total = Q_R_T(end);
    fprintf('Q_R    = %.2f kW\n', Q_R_total);

    Q_cold_total = Q_HENE_total + Q_R_total;
    fprintf('Q_cold = %.2f kW\n', Q_cold_total);
    duty_mismatch = Q_cold_total - Q_hot_total;
    fprintf('Duty mismatch (Q_cold - Q_hot) = %+.2f kW (%+.2f%%)\n\n', ...
            duty_mismatch, 100*duty_mismatch/Q_hot_total);

    %---- Cold composite (as-is, no rescaling) -----------------------------
    [T_comp, Q_cold_comp] = compute_composite(T_R_grid, Q_R_T, ...
                                              T_HENE_grid, Q_HENE_T, ...
                                              N_comp);

    %---- Pinch on common Q overlap ----------------------------------------
    [dT_min, Q_pinch, T_hot_p, T_cold_p] = pinch_metrics( ...
        Q_hot_T, T_hot_grid, Q_cold_comp, T_comp);

    dT_hot_end  = T_hot_in  - T_HENE_out;
    dT_cold_end = T_hot_out - min(T_R_in, T_HENE_in);

    %---- Report -----------------------------------------------------------
    fprintf('========================================================\n');
    fprintf('BASELINE PFHX-5 PINCH EVALUATION (no solver)\n');
    fprintf('========================================================\n');
    fprintf('Stream duties\n');
    fprintf('  Q_hot          = %.2f kW\n', Q_hot_total);
    fprintf('  Q_HENE         = %.2f kW\n', Q_HENE_total);
    fprintf('  Q_R            = %.2f kW\n', Q_R_total);
    fprintf('  Q_cold (total) = %.2f kW\n', Q_cold_total);
    fprintf('  Mismatch       = %+.2f kW  (%+.2f%%)\n', ...
            duty_mismatch, 100*duty_mismatch/Q_hot_total);
    fprintf('\nEndpoint temperature differences\n');
    fprintf('  dT_hot_end  (T_H6 - T_HENE5)   = %.3f K\n', dT_hot_end);
    fprintf('  dT_cold_end (T_H7 - T_cold_lo) = %.3f K\n', dT_cold_end);
    fprintf('\nGlobal pinch (common Q overlap)\n');
    fprintf('  dT_min        = %.3f K  at Q = %.1f kW\n', dT_min, Q_pinch);
    fprintf('    T_hot  pinch = %.2f K\n', T_hot_p);
    fprintf('    T_cold pinch = %.2f K\n', T_cold_p);
    fprintf('========================================================\n\n');

    %---- Pack output ------------------------------------------------------
    out.T_hot_grid   = T_hot_grid;
    out.Q_hot_T      = Q_hot_T;
    out.Q_hot_total  = Q_hot_total;
    out.Q_HENE_total = Q_HENE_total;
    out.Q_R_total    = Q_R_total;
    out.Q_cold_total = Q_cold_total;
    out.T_comp       = T_comp;
    out.Q_cold_comp  = Q_cold_comp;
    out.mdot_HENE    = mdot_HENE;
    out.T_HENE_out   = T_HENE_out;
    out.dT_min       = dT_min;
    out.Q_pinch      = Q_pinch;
    out.T_hot_p      = T_hot_p;
    out.T_cold_p     = T_cold_p;
end


function out = run_ejector(spec, DT_MIN_TARGET, N_hot, N_cold, N_comp)
% RUN_EJECTOR  Solve PFHX-5 in ejector-modified configuration.
%
% Free variables: mdot_HENE, T_HENE_out.
% Constraints:    Q_HENE + Q_R = Q_hot (duty closure),
%                 min(dT) on hot vs cold composite = DT_MIN_TARGET.
%
% spec struct fields (all hot- and R-side stream specs come from caller):
%   P_hot      [bar]    hot-side pressure (= ejector discharge)
%   T_hot_in   [K]      hot-side inlet temperature (= diffuser outlet)
%   T_hot_out  [K]      hot-side outlet temperature (= T_postPFHX5, fixed)
%   mdot_H     [kg/s]   hot-side mass flow (= mDIF)
%   mdot_R     [kg/s]   R-stream cold-side mass flow (= mR, recycle vapor)
%   T_R_in     [K]      R-stream inlet temperature (saturated vapor at Psep)
%   T_HENE_in  [K]      He-Ne inlet temperature (per-case, optional;
%                       defaults to 27.45 K Moro baseline if absent)
%   Q_PFHX5_spec [kW]   ejector-model duty estimate (printed for comparison)
%
% HENE pressure and He fraction are common across cases and set inside.

    %---- Hot-side and R-side stream specs (from spec) ---------------------
    P_hot      = spec.P_hot;
    T_hot_in   = spec.T_hot_in;
    T_hot_out  = spec.T_hot_out;
    mdot_H     = spec.mdot_H;
    Q_hot_spec = spec.Q_PFHX5_spec;

    P_R    = 1.0;                 % bar
    T_R_in = spec.T_R_in;
    T_R_out= 46.65;               % K, fixed across ejector cases
    mdot_R = spec.mdot_R;
    xp_R   = 0.999;

    %---- Defense-in-depth: guard against sub-saturation T_R_in -----------
    % At T = Tsat(P_R), REFPROP returns the saturated-LIQUID enthalpy for
    % parahydrogen. Sampling enthalpy at exactly Tsat then jumping to vapor
    % at the next grid point gives a phantom 445 kJ/kg latent-heat spike
    % that inflates Q_R by ~140 kW and makes the cold composite unphysical.
    % If the caller passes a T_R_in too close to Tsat (or below), bump it
    % to a safe vapor-side value matching the baseline (20.40 K at 1 bar).
    T_R_in_min_safe = 20.40;
    if T_R_in < T_R_in_min_safe
        fprintf(['WARNING: T_R_in = %.4f K is below the safe vapor-side ' ...
                 'value (%.2f K).\n         Bumping T_R_in to %.2f K to ' ...
                 'avoid REFPROP latent-heat spike.\n'], ...
                T_R_in, T_R_in_min_safe, T_R_in_min_safe);
        T_R_in = T_R_in_min_safe;
    end

    % HENE inlet (free: mdot, T_out)
    %
    % T_HENE_in is read from the per-case spec (column T_HENE_in_K in
    % pinch_inputs.csv). The lowered values were chosen to make
    % dT_min = 2.48 K feasible at each case while still being only ~1 K
    % below the Moro 27.45 K reference.
    P_HENE = 3.0;            % bar
    if isfield(spec, 'T_HENE_in') && ~isempty(spec.T_HENE_in)
        T_HENE_in = spec.T_HENE_in;
    else
        T_HENE_in = 27.45;   % K, fall back to Moro default if not provided
    end
    x_He = 0.80; x_Ne = 0.20;

    %---- Hot-side enthalpy table ------------------------------------------
    use_equilhyd = test_equilhyd();
    [T_hot_grid, h_hot] = build_hot(T_hot_in, T_hot_out, P_hot, ...
                                    N_hot, use_equilhyd);
    Q_hot_T     = mdot_H * (h_hot - h_hot(1)) / 1000;
    Q_hot_total = Q_hot_T(end);
    fprintf('Hot duty = %.2f kW (spec %.2f kW, diff %.2f%%)\n', ...
            Q_hot_total, Q_hot_spec, ...
            100*(Q_hot_total - Q_hot_spec)/Q_hot_spec);

    %---- R stream enthalpy table ------------------------------------------
    T_R_grid = linspace(T_R_in, T_R_out, N_cold).';
    h_R = zeros(N_cold,1);
    for i = 1:N_cold
        h_R(i) = refpropm('H','T',T_R_grid(i),'P',P_R*100, ...
                          'orthohyd','parahyd',[1-xp_R, xp_R]);
    end
    Q_R_T     = mdot_R * (h_R - h_R(1)) / 1000;
    Q_R_total = Q_R_T(end);
    fprintf('R duty   = %.2f kW\n', Q_R_total);

    Q_HENE_required = Q_hot_total - Q_R_total;
    fprintf('HENE duty required = %.2f kW\n\n', Q_HENE_required);

    %---- Bisection on T_HENE_out ------------------------------------------
    T_HENE_out_min = T_HENE_in + 0.5;
    T_HENE_out_max = T_hot_in - DT_MIN_TARGET;

    fprintf('Bisecting T_HENE_out in [%.2f, %.2f] K for dT_min = %.2f K\n', ...
            T_HENE_out_min, T_HENE_out_max, DT_MIN_TARGET);

    f_lo = pinch_residual(T_HENE_out_min, DT_MIN_TARGET, ...
                          Q_HENE_required, T_HENE_in, P_HENE, x_He, x_Ne, ...
                          N_cold, T_R_grid, Q_R_T, T_R_in, N_comp, ...
                          T_hot_grid, Q_hot_T, Q_hot_total);
    f_hi = pinch_residual(T_HENE_out_max, DT_MIN_TARGET, ...
                          Q_HENE_required, T_HENE_in, P_HENE, x_He, x_Ne, ...
                          N_cold, T_R_grid, Q_R_T, T_R_in, N_comp, ...
                          T_hot_grid, Q_hot_T, Q_hot_total);

    fprintf('  at T_out = %.2f K: f = %+.3f K (min_dT = %.3f K)\n', ...
            T_HENE_out_min, f_lo, f_lo + DT_MIN_TARGET);
    fprintf('  at T_out = %.2f K: f = %+.3f K (min_dT = %.3f K)\n', ...
            T_HENE_out_max, f_hi, f_hi + DT_MIN_TARGET);

    if f_lo * f_hi > 0
        if f_lo < 0 && f_hi < 0
            fprintf(['Target %.2f K infeasible. Searching for T_HENE_out ' ...
                     'that maximizes pinch...\n'], DT_MIN_TARGET);
            objective = @(T) -(pinch_residual(T, 0, Q_HENE_required, ...
                T_HENE_in, P_HENE, x_He, x_Ne, N_cold, T_R_grid, Q_R_T, ...
                T_R_in, N_comp, T_hot_grid, Q_hot_T, Q_hot_total));
            T_HENE_out_solved = fminbnd(objective, T_HENE_out_min, ...
                T_HENE_out_max, optimset('TolX',1e-3,'Display','off'));
            fprintf('Max-pinch T_HENE_out = %.3f K\n', T_HENE_out_solved);
        else
            T_HENE_out_solved = T_HENE_out_max;
            fprintf('Both endpoints feasible. Hot-end cap binds.\n');
        end
    else
        a = T_HENE_out_min; b = T_HENE_out_max;
        for it = 1:50
            c = 0.5*(a+b);
            fc = pinch_residual(c, DT_MIN_TARGET, ...
                                Q_HENE_required, T_HENE_in, P_HENE, ...
                                x_He, x_Ne, N_cold, T_R_grid, Q_R_T, ...
                                T_R_in, N_comp, T_hot_grid, Q_hot_T, ...
                                Q_hot_total);
            if f_lo * fc < 0
                b = c; f_hi = fc; %#ok<NASGU>
            else
                a = c; f_lo = fc;
            end
            if abs(b-a) < 1e-3, break; end
        end
        T_HENE_out_solved = 0.5*(a+b);
        fprintf('Converged in %d iterations: T_HENE_out = %.3f K\n', ...
                it, T_HENE_out_solved);
    end

    %---- Rebuild HENE at solved T_HENE_out --------------------------------
    [T_HENE_grid, Q_HENE_T, mdot_HENE_solved] = build_HENE( ...
        T_HENE_out_solved, Q_HENE_required, T_HENE_in, P_HENE, ...
        x_He, x_Ne, N_cold);
    Q_HENE_total = Q_HENE_T(end);

    %---- Cold composite ---------------------------------------------------
    [T_comp, Q_cold_comp] = compute_composite(T_R_grid, Q_R_T, ...
                                              T_HENE_grid, Q_HENE_T, ...
                                              N_comp);
    Q_cold_total = Q_cold_comp(end);

    %---- Pinch on common Q overlap ----------------------------------------
    [dT_min, Q_pinch, T_hot_p, T_cold_p] = pinch_metrics( ...
        Q_hot_T, T_hot_grid, Q_cold_comp, T_comp);

    dT_hot_end  = T_hot_in  - T_HENE_out_solved;
    dT_cold_end = T_hot_out - min(T_R_in, T_HENE_in);

    %---- Report -----------------------------------------------------------
    fprintf('\n========================================================\n');
    fprintf('SOLUTION (target dT_min = %.2f K)\n', DT_MIN_TARGET);
    fprintf('========================================================\n');
    fprintf('  T_HENE_out     = %.3f K\n', T_HENE_out_solved);
    fprintf('  mdot_HENE      = %.4f kg/s  (baseline 6.9, ratio %.2f)\n', ...
            mdot_HENE_solved, mdot_HENE_solved/6.9);
    fprintf('  Q_HENE         = %.2f kW\n', Q_HENE_total);
    fprintf('  Q_R            = %.2f kW\n', Q_R_total);
    fprintf('  Q_cold total   = %.2f kW\n', Q_HENE_total + Q_R_total);
    fprintf('  Q_hot          = %.2f kW\n', Q_hot_total);
    fprintf('  Duty closure   = %+.2f kW  (%+.2f%%)\n', ...
            Q_cold_total - Q_hot_total, ...
            100*(Q_cold_total - Q_hot_total)/Q_hot_total);
    fprintf('  dT_hot_end     = %.3f K\n', dT_hot_end);
    fprintf('  dT_cold_end    = %.3f K\n', dT_cold_end);
    fprintf('  dT_min (pinch) = %.3f K  at Q = %.1f kW\n', dT_min, Q_pinch);
    fprintf('    T_hot  pinch = %.2f K\n', T_hot_p);
    fprintf('    T_cold pinch = %.2f K\n', T_cold_p);
    fprintf('========================================================\n\n');

    %---- Pack output ------------------------------------------------------
    out.T_hot_grid   = T_hot_grid;
    out.Q_hot_T      = Q_hot_T;
    out.Q_hot_total  = Q_hot_total;
    out.Q_HENE_total = Q_HENE_total;
    out.Q_R_total    = Q_R_total;
    out.Q_cold_total = Q_cold_total;
    out.T_comp       = T_comp;
    out.Q_cold_comp  = Q_cold_comp;
    out.mdot_HENE    = mdot_HENE_solved;
    out.T_HENE_out   = T_HENE_out_solved;
    out.dT_min       = dT_min;
    out.Q_pinch      = Q_pinch;
    out.T_hot_p      = T_hot_p;
    out.T_cold_p     = T_cold_p;
end

%% ======================================================================
%  HELPER FUNCTIONS
%  ======================================================================

function tf = test_equilhyd()
% TEST_EQUILHYD  Return true if REFPROP supports the 'equilhyd' fluid.
    try
        h = refpropm('H','T',40,'P',10*100,'equilhyd'); %#ok<NASGU>
        tf = true;
    catch
        tf = false;
    end
end

function xp = xp_equilibrium(T)
% XP_EQUILIBRIUM  Para-fraction at thermodynamic equilibrium (Farkas /
% Leachman 2009 tabulated, pchip-interpolated). Below 20 K: xp = 1.
% Above 300 K: xp = 0.25 (high-T limit).
    T_tab = [ 20  25  30  35  40  45  50  60  70  80  90 100 ...
             120 150 200 250 300 ].';
    x_tab = [0.9989 0.9911 0.9702 0.9316 0.8797 0.8193 0.7558 ...
             0.6309 0.5271 0.4487 0.3913 0.3496 0.2984 0.2638 ...
             0.2507 0.2500 0.2500 ].';
    if T < T_tab(1)
        xp = 1.0;
    elseif T > T_tab(end)
        xp = 0.25;
    else
        xp = interp1(T_tab, x_tab, T, 'pchip');
    end
end

function [T_grid, h_hot] = build_hot(T_in, T_out, P_bar, N, use_equilhyd)
% BUILD_HOT  Hot-side enthalpy table at local equilibrium para-fraction.
% Grid runs from T_out (cold end) to T_in (hot end).
    T_grid = linspace(T_out, T_in, N).';
    h_hot  = zeros(N,1);
    fprintf('Building hot-side enthalpy table (%d points)...\n', N);
    if use_equilhyd
        for i = 1:N
            h_hot(i) = refpropm('H','T',T_grid(i),'P',P_bar*100,'equilhyd');
        end
    else
        for i = 1:N
            xp_i = xp_equilibrium(T_grid(i));
            h_hot(i) = refpropm('H','T',T_grid(i),'P',P_bar*100, ...
                                'orthohyd','parahyd',[1-xp_i, xp_i]);
        end
    end
end

function [T_grid, Q_T, mdot] = build_HENE(T_out, Q_required, T_in, P, ...
                                          xHe, xNe, N)
% BUILD_HENE  HENE enthalpy table sized to deliver exactly Q_required
% between T_in and T_out. Returns the implied mass flow rate.
    T_grid = linspace(T_in, T_out, N).';
    h = zeros(N,1);
    for i = 1:N
        h(i) = refpropm('H','T',T_grid(i),'P',P*100, ...
                        'helium','neon',[xHe,xNe]);
    end
    dh   = h(end) - h(1);
    mdot = Q_required*1000 / dh;
    Q_T  = mdot * (h - h(1)) / 1000;
end

function f = pinch_residual(T_HENE_out, DT_target, Q_HENE_req, T_HENE_in, ...
                            P_HENE, xHe, xNe, N_cold, T_R_grid, Q_R_T, ...
                            T_R_in, N_comp, T_hot_grid, Q_hot_T, ...
                            Q_hot_total)
% PINCH_RESIDUAL  f = min_dT(T_HENE_out) - DT_target. Used by the bisection
% / fminbnd solver in run_ejector. Duty closure is enforced inside
% build_HENE which sizes mdot_HENE to absorb exactly Q_HENE_req.
    [T_HENE_grid, Q_HENE_T, ~] = build_HENE(T_HENE_out, Q_HENE_req, ...
        T_HENE_in, P_HENE, xHe, xNe, N_cold);

    T_cold_lo = min(T_R_in, T_HENE_in);
    T_cold_hi = max(T_R_grid(end), T_HENE_out);
    T_comp = linspace(T_cold_lo, T_cold_hi, N_comp).';

    Q_R_on    = interp_clipped(T_R_grid,    Q_R_T,    T_comp);
    Q_HENE_on = interp_clipped(T_HENE_grid, Q_HENE_T, T_comp);
    Q_cold    = cummax(Q_R_on + Q_HENE_on);

    Q_max = min(Q_hot_total, Q_cold(end));
    Qg    = linspace(0, Q_max, 2001).';
    T_h   = interp1(Q_hot_T, T_hot_grid, Qg, 'pchip');
    T_c   = interp1(Q_cold,  T_comp,     Qg, 'pchip');

    min_dT = min(T_h - T_c);
    f      = min_dT - DT_target;
end

function Qout = interp_clipped(T_in, Q_in, T_query)
% INTERP_CLIPPED  Interpolate Q(T) for query points inside the source
% range; clip to 0 below the source and to Q_in(end) above. Used to
% combine streams that span different T ranges into one composite.
    Qout = zeros(size(T_query));
    inside = T_query >= T_in(1) & T_query <= T_in(end);
    Qout(inside) = interp1(T_in, Q_in, T_query(inside), 'pchip');
    Qout(T_query > T_in(end)) = Q_in(end);
end

function [T_comp, Q_cold_comp] = compute_composite(T_R_grid, Q_R_T, ...
                                                   T_HENE_grid, Q_HENE_T, ...
                                                   N_comp)
% COMPUTE_COMPOSITE  Build the cold composite curve on a common T grid
% spanning both cold streams. cummax enforces monotonicity in case
% interpolation produces a tiny non-monotonic step.
    T_cold_lo = min(T_R_grid(1),    T_HENE_grid(1));
    T_cold_hi = max(T_R_grid(end),  T_HENE_grid(end));
    T_comp    = linspace(T_cold_lo, T_cold_hi, N_comp).';

    Q_R_on_comp    = interp_clipped(T_R_grid,    Q_R_T,    T_comp);
    Q_HENE_on_comp = interp_clipped(T_HENE_grid, Q_HENE_T, T_comp);
    Q_cold_comp    = cummax(Q_R_on_comp + Q_HENE_on_comp);
end

function [dT_min, Q_pinch, T_hot_p, T_cold_p] = pinch_metrics( ...
    Q_hot_T, T_hot_grid, Q_cold_comp, T_comp)
% PINCH_METRICS  Global pinch on the common Q overlap of hot and cold
% composites. Returns dT_min and the (Q, T_hot, T_cold) at the pinch.
    Q_max_common = min(Q_hot_T(end), Q_cold_comp(end));
    Q_grid       = linspace(0, Q_max_common, 4001).';
    T_hot_of_Q   = interp1(Q_hot_T,     T_hot_grid, Q_grid, 'pchip');
    T_cold_of_Q  = interp1(Q_cold_comp, T_comp,     Q_grid, 'pchip');
    dT_profile   = T_hot_of_Q - T_cold_of_Q;
    [dT_min, i_pinch] = min(dT_profile);
    Q_pinch  = Q_grid(i_pinch);
    T_hot_p  = T_hot_of_Q(i_pinch);
    T_cold_p = T_cold_of_Q(i_pinch);
end

function write_curves_csv(filename, Q_hot, T_hot, Q_cold, T_cold)
% WRITE_CURVES_CSV  Pad the shorter vector with NaN and write the
% rectangular table that plot_pinch_PFHX5_final.py expects.
    n_hot  = length(Q_hot);
    n_cold = length(Q_cold);
    n_max  = max(n_hot, n_cold);

    Q_hot_pad  = [Q_hot;  NaN(n_max - n_hot,  1)];
    T_hot_pad  = [T_hot;  NaN(n_max - n_hot,  1)];
    Q_cold_pad = [Q_cold; NaN(n_max - n_cold, 1)];
    T_cold_pad = [T_cold; NaN(n_max - n_cold, 1)];

    T = table(Q_hot_pad, T_hot_pad, Q_cold_pad, T_cold_pad, ...
        'VariableNames', {'Q_hot_kW','T_hot_K','Q_cold_kW','T_cold_K'});
    writetable(T, filename);
end
