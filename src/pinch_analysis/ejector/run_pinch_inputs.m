%% Generate PFHX-5 pinch-analysis inputs from ejector_solve
%
%  Calls ejector_solve at two operating points:
%     case 1:  f = 0.20, eta = 0.20  (existing pinch low-eta point)
%     case 2:  f = 0.10, eta = 0.30  (new high-eta point, max accessible Pd)
%
%  For each case the seven stream values needed by pinch_PFHX5.m are
%  packed into a CSV (one row per case, header on top):
%
%      case_label, P_hot_bar, T_hot_in_K, T_hot_out_K, mdot_H_kgps,
%      mdot_R_kgps, T_R_in_K, Q_PFHX5_kW
%
%  T_hot_out is fixed at the model's T_postPFHX5 input (27.95 K).
%  T_R_in is computed as Tsat at Psep so its provenance is explicit;
%  it should reproduce ~20.23 K for parahydrogen at 1 bar.
%
%  Verification: at f=0.20, eta=0.20 the v13 self-consistent loop will
%  converge to different values than the v12 single-pass script
%  (v12 reported x_refinement_error = 0.0857 at this point). v13 is
%  the corrected result. Compare:
%      Reference v12 at f=0.20, eta=0.20 (UNREFINED):
%        Pd = 26.386 bar, Tdiffusor = 64.86 K, mDIF = 1.5897 kg/s,
%        mR = 0.3078 kg/s, x = 0.242029, mL = 1.2050 kg/s,
%        Q_PFHX5 = 1166.30 kW

clear
clc

% =====================================================================
%  Operating points to evaluate
%
%  T_HENE_in_K is set per-case at the threshold value where dT_min = 2.48 K
%  becomes feasible (from the T_HENE_in sweep done separately):
%     low_eta:  T_HENE_in lowered by 1.25 K from the Moro 27.45 K -> 26.20 K
%     high_eta: T_HENE_in lowered by 0.75 K -> 26.70 K
%  Both values sit on the threshold of feasibility for the baseline-matching
%  dT_min = 2.48 K target. Lower values would give more pinch margin but
%  cost more in the upstream He-Ne Brayton stage.
% =====================================================================
cases = struct( ...
    'label',     {'low_eta',  'high_eta'}, ...
    'f',         {0.20,        0.10}, ...
    'eta',       {0.20,        0.30}, ...
    'T_HENE_in', {26.20,       26.70});

% Solver options shared by both runs (defaults match Moro baseline,
% see ejector_solve.m header for the full list)
opts = struct('verbose', true);

n_cases = numel(cases);

% =====================================================================
%  Run both operating points
%  Note: do NOT preallocate results as struct(); MATLAB will reject
%  later assignments because the empty struct schema differs from
%  the populated one returned by ejector_solve. Grow on first assign.
% =====================================================================
for k = 1:n_cases
    fprintf('\n############################################################\n');
    fprintf('# Case %d/%d : %s   (f = %.3f, eta = %.3f)\n', ...
            k, n_cases, cases(k).label, cases(k).f, cases(k).eta);
    fprintf('############################################################\n');

    r = ejector_solve(cases(k).f, cases(k).eta, opts);

    if ~r.converged
        error('Case %s did not converge: %s', ...
              cases(k).label, r.failure_reason);
    end

    if k == 1
        results = r;
    else
        results(k) = r;
    end
end

% =====================================================================
%  T_R_in: must be ABOVE saturated vapor temperature at separator
%  pressure. At Tsat itself (~20.23 K for pH2 at 1 bar) REFPROP returns
%  the saturated-liquid enthalpy, not the vapor enthalpy, which causes a
%  phantom 445 kJ/kg latent-heat jump when the next grid point lands in
%  the vapor region. The pinch baseline raises T_R_in to 20.40 K for
%  exactly this reason; we apply the same offset here.
%
%  T_R_in = Tsat + T_VAPOR_OFFSET, with T_VAPOR_OFFSET set to ~0.17 K
%  so the value matches the baseline 20.40 K.
% =====================================================================
T_VAPOR_OFFSET = 0.17;                  % K, vapor-side margin above Tsat

Psep_kPa = results(1).Psep_kPa;         % all cases use the same default
T_sat    = refpropm('T', 'P', Psep_kPa, 'Q', 1, 'parahydrogen');
T_R_in   = T_sat + T_VAPOR_OFFSET;
fprintf(['\nT_sat (P=%.1f kPa, Q=1, parahydrogen) = %.4f K\n' ...
         'T_R_in = T_sat + %.2f K = %.4f K  (vapor-side margin to ' ...
         'avoid latent-heat spike)\n'], ...
        Psep_kPa, T_sat, T_VAPOR_OFFSET, T_R_in);

% =====================================================================
%  Build the seven-column row per case and assemble table
% =====================================================================
T_hot_out_K = results(1).T_postPFHX5;   % 27.95 K, identical for both cases

case_label    = cell(n_cases, 1);
P_hot_bar     = zeros(n_cases, 1);
T_hot_in_K    = zeros(n_cases, 1);
T_hot_out_col = zeros(n_cases, 1);
mdot_H_kgps   = zeros(n_cases, 1);
mdot_R_kgps   = zeros(n_cases, 1);
T_R_in_K      = zeros(n_cases, 1);
T_HENE_in_K   = zeros(n_cases, 1);
Q_PFHX5_kW    = zeros(n_cases, 1);

for k = 1:n_cases
    case_label{k}     = cases(k).label;
    P_hot_bar(k)      = results(k).Pd_bar;
    T_hot_in_K(k)     = results(k).Tdiffusor;
    T_hot_out_col(k)  = T_hot_out_K;
    mdot_H_kgps(k)    = results(k).mDIF;
    mdot_R_kgps(k)    = results(k).mR;
    T_R_in_K(k)       = T_R_in;
    T_HENE_in_K(k)    = cases(k).T_HENE_in;
    Q_PFHX5_kW(k)     = results(k).Q_PFHX5_kW;
end

T_out = table( ...
    case_label, P_hot_bar, T_hot_in_K, T_hot_out_col, ...
    mdot_H_kgps, mdot_R_kgps, T_R_in_K, T_HENE_in_K, Q_PFHX5_kW, ...
    'VariableNames', {'case_label', 'P_hot_bar', 'T_hot_in_K', ...
                      'T_hot_out_K', 'mdot_H_kgps', 'mdot_R_kgps', ...
                      'T_R_in_K', 'T_HENE_in_K', 'Q_PFHX5_kW'});

writetable(T_out, 'pinch_inputs.csv');
fprintf('\nWrote pinch_inputs.csv\n');

% =====================================================================
%  Side-by-side console summary
% =====================================================================
fprintf('\n############################################################\n');
fprintf('# PINCH ANALYSIS INPUTS\n');
fprintf('############################################################\n');
fprintf('%-22s  %12s  %12s\n', 'quantity', cases(1).label, cases(2).label);
fprintf('%-22s  %12.4f  %12.4f\n', 'f',           cases(1).f,           cases(2).f);
fprintf('%-22s  %12.4f  %12.4f\n', 'eta',         cases(1).eta,         cases(2).eta);
fprintf('%-22s  %12.4f  %12.4f\n', 'P_hot [bar]', P_hot_bar(1),         P_hot_bar(2));
fprintf('%-22s  %12.4f  %12.4f\n', 'T_hot_in [K]', T_hot_in_K(1),       T_hot_in_K(2));
fprintf('%-22s  %12.4f  %12.4f\n', 'T_hot_out [K]', T_hot_out_col(1),   T_hot_out_col(2));
fprintf('%-22s  %12.4f  %12.4f\n', 'mdot_H [kg/s]', mdot_H_kgps(1),     mdot_H_kgps(2));
fprintf('%-22s  %12.4f  %12.4f\n', 'mdot_R [kg/s]', mdot_R_kgps(1),     mdot_R_kgps(2));
fprintf('%-22s  %12.4f  %12.4f\n', 'T_R_in [K]',   T_R_in_K(1),         T_R_in_K(2));
fprintf('%-22s  %12.4f  %12.4f\n', 'T_HENE_in [K]', T_HENE_in_K(1),     T_HENE_in_K(2));
fprintf('%-22s  %12.4f  %12.4f\n', 'Q_PFHX5 [kW]', Q_PFHX5_kW(1),       Q_PFHX5_kW(2));
fprintf('############################################################\n\n');

% =====================================================================
%  Verification reminder for the f=0.20, eta=0.20 case
% =====================================================================
fprintf('Reference v12 values at f=0.20, eta=0.20 (UNREFINED single-pass):\n');
fprintf('  Pd        = 26.386 bar    (v13 result: %.3f bar)\n',   P_hot_bar(1));
fprintf('  Tdiffusor = 64.86  K      (v13 result: %.3f K)\n',     T_hot_in_K(1));
fprintf('  mDIF      = 1.5897 kg/s   (v13 result: %.4f kg/s)\n',  mdot_H_kgps(1));
fprintf('  mR        = 0.3078 kg/s   (v13 result: %.4f kg/s)\n',  mdot_R_kgps(1));
fprintf('  Q_PFHX5   = 1166.30 kW    (v13 result: %.2f kW)\n',    Q_PFHX5_kW(1));
fprintf(['\nNon-trivial differences are expected because v12 reported ' ...
         'x_refinement_error = 0.0857 at this point.\nv13 is the ' ...
         'self-consistent fixed-point result and is what feeds the ' ...
         'pinch analysis going forward.\n\n']);
