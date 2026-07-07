% pinch_nelium90_comparison.m
% -------------------------------------------------------------------
% Compares PFHX-5 composite curves for two He-Ne compositions:
%
%   Case 1 (baseline):  80 % He / 20 % Ne at 3 bar
%     -- Fully specified Moro et al. streams (mdot_HENE = 6.9 kg/s,
%        T_HENE_out = 74.05 K), same as run_baseline in pinch_PFHX5.m.
%        dT_min is an output (~2.5 K).
%
%   Case 2 (Nelium-90): 10 % He / 90 % Ne at 10 bar
%     -- T_HENE_in = 27.45 K fixed. Bisection solver on T_HENE_out
%        to hit dT_min = 2.5 K (same thermal margin as baseline).
%        mdot_HENE from duty closure at each trial T_HENE_out.
%
% Same hot side (H2 at 75 bar) and R stream in both cases.
%
% Output: nelium90_composite_baseline.csv
%         nelium90_composite_nelium90.csv
%         nelium90_summary.csv
% -------------------------------------------------------------------

clear; clc; close all;

%% ---- Configuration ----------------------------------------------------
N_hot  = 401;
N_cold = 401;
N_comp = 801;
DT_TARGET = 2.48;   % K, match the existing baseline pinch target

%% ---- Hot side (H2, common to both cases) -------------------------------
P_hot     = 75.0;       % bar
T_hot_in  = 74.15;      % K (H6)
T_hot_out = 27.95;      % K (H7)
mdot_H    = 1.4847;     % kg/s

use_equilhyd = test_equilhyd();
[T_hot_grid, h_hot] = build_hot(T_hot_in, T_hot_out, P_hot, ...
                                N_hot, use_equilhyd);
Q_hot_T     = mdot_H * (h_hot - h_hot(1)) / 1000;
Q_hot_total = Q_hot_T(end);
fprintf('Hot side: Q_hot = %.2f kW\n', Q_hot_total);

%% ---- R stream (common to both cases) -----------------------------------
P_R    = 1.0;           % bar
T_R_in = 20.40;         % K
T_R_out= 46.65;         % K
mdot_R = 0.4847;        % kg/s
xp_R   = 0.999;

T_R_grid = linspace(T_R_in, T_R_out, N_cold).';
h_R = zeros(N_cold, 1);
for i = 1:N_cold
    h_R(i) = refpropm('H','T',T_R_grid(i),'P',P_R*100, ...
                      'orthohyd','parahyd',[1-xp_R, xp_R]);
end
Q_R_T     = mdot_R * (h_R - h_R(1)) / 1000;
Q_R_total = Q_R_T(end);
fprintf('R stream: Q_R = %.2f kW\n\n', Q_R_total);

%% ========================================================================
%  CASE 1: BASELINE  (80 % He / 20 % Ne, 3 bar, fully specified)
%  ========================================================================
fprintf('============================================================\n');
fprintf('  CASE 1: Baseline (80%% He, 20%% Ne, P = 3 bar)\n');
fprintf('============================================================\n');

P_HENE_1     = 3.0;
T_HENE_in_1  = 27.45;
T_HENE_out_1 = 74.05;
mdot_HENE_1  = 6.9;       % kg/s, fixed Moro value
xHe_1 = 0.80;  xNe_1 = 0.20;

% Build HENE enthalpy table at fixed mdot (not from duty closure)
T_HENE_grid_1 = linspace(T_HENE_in_1, T_HENE_out_1, N_cold).';
h_HENE_1 = zeros(N_cold, 1);
for i = 1:N_cold
    h_HENE_1(i) = refpropm('H','T',T_HENE_grid_1(i),'P',P_HENE_1*100, ...
                           'helium','neon',[xHe_1, xNe_1]);
end
Q_HENE_T_1     = mdot_HENE_1 * (h_HENE_1 - h_HENE_1(1)) / 1000;
Q_HENE_total_1 = Q_HENE_T_1(end);

Q_cold_total_1 = Q_HENE_total_1 + Q_R_total;
mismatch_1 = Q_cold_total_1 - Q_hot_total;
fprintf('  Q_HENE = %.2f kW,  Q_cold = %.2f kW\n', Q_HENE_total_1, Q_cold_total_1);
fprintf('  Duty mismatch = %+.2f kW (%+.2f%%)\n', ...
        mismatch_1, 100*mismatch_1/Q_hot_total);

% Composite and pinch
[T_comp_1, Q_comp_1] = compute_composite( ...
    T_R_grid, Q_R_T, T_HENE_grid_1, Q_HENE_T_1, N_comp);
[dTmin_1, Qp_1, Thp_1, Tcp_1] = pinch_metrics( ...
    Q_hot_T, T_hot_grid, Q_comp_1, T_comp_1);
fprintf('  dT_min = %.3f K  at Q = %.1f kW\n', dTmin_1, Qp_1);
fprintf('  T_hot_p = %.2f K,  T_cold_p = %.2f K\n\n', Thp_1, Tcp_1);

%% ========================================================================
%  CASE 2: NELIUM-90  (10 % He / 90 % Ne, 10 bar, solver)
%  ========================================================================
fprintf('============================================================\n');
fprintf('  CASE 2: Nelium-90 (10%% He, 90%% Ne, P = 10 bar)\n');
fprintf('============================================================\n');

P_HENE_2    = 10.0;
T_HENE_in_2 = 27.45;
xHe_2 = 0.10;  xNe_2 = 0.90;

% --- VLE boundaries ---
T_bub = NaN;  T_dew = NaN;  h_bub = NaN;  h_dew = NaN;
try
    T_bub = refpropm('T','P',P_HENE_2*100,'Q',0, ...
                     'helium','neon',[xHe_2,xNe_2]);
    T_dew = refpropm('T','P',P_HENE_2*100,'Q',1, ...
                     'helium','neon',[xHe_2,xNe_2]);
    h_bub = refpropm('H','P',P_HENE_2*100,'Q',0, ...
                     'helium','neon',[xHe_2,xNe_2]);
    h_dew = refpropm('H','P',P_HENE_2*100,'Q',1, ...
                     'helium','neon',[xHe_2,xNe_2]);
    fprintf('  VLE at %.0f bar:\n', P_HENE_2);
    fprintf('    T_bubble = %.2f K   h_bubble = %.1f J/kg\n', T_bub, h_bub);
    fprintf('    T_dew    = %.2f K   h_dew    = %.1f J/kg\n', T_dew, h_dew);
    fprintf('    Glide    = %.2f K\n\n', T_dew - T_bub);
catch ME
    fprintf('  WARNING: VLE boundaries failed: %s\n\n', ME.message);
end

% --- Duty required for He-Ne ---
Q_HENE_req = Q_hot_total - Q_R_total;
fprintf('  Q_HENE required = %.2f kW\n', Q_HENE_req);

% --- Bisection on T_HENE_out for dT_min = DT_TARGET ---
T_out_lo = T_HENE_in_2 + 0.5;
T_out_hi = T_hot_in - DT_TARGET;

fprintf('  Bisecting T_HENE_out in [%.2f, %.2f] K for dT_min = %.2f K\n', ...
        T_out_lo, T_out_hi, DT_TARGET);

f_lo = nel_pinch_residual(T_out_lo, DT_TARGET, Q_HENE_req, ...
    T_HENE_in_2, P_HENE_2, xHe_2, xNe_2, N_cold, ...
    T_R_grid, Q_R_T, T_R_in, N_comp, T_hot_grid, Q_hot_T, Q_hot_total, ...
    T_bub, T_dew, h_bub, h_dew);
f_hi = nel_pinch_residual(T_out_hi, DT_TARGET, Q_HENE_req, ...
    T_HENE_in_2, P_HENE_2, xHe_2, xNe_2, N_cold, ...
    T_R_grid, Q_R_T, T_R_in, N_comp, T_hot_grid, Q_hot_T, Q_hot_total, ...
    T_bub, T_dew, h_bub, h_dew);

fprintf('  f(%.2f) = %+.3f,  f(%.2f) = %+.3f\n', T_out_lo, f_lo, T_out_hi, f_hi);

if f_lo * f_hi > 0
    if f_lo < 0 && f_hi < 0
        fprintf('  Target %.2f K infeasible. Searching for max-pinch T_HENE_out...\n', ...
                DT_TARGET);
        objective = @(T) -(nel_pinch_residual(T, 0, Q_HENE_req, ...
            T_HENE_in_2, P_HENE_2, xHe_2, xNe_2, N_cold, ...
            T_R_grid, Q_R_T, T_R_in, N_comp, T_hot_grid, Q_hot_T, ...
            Q_hot_total, T_bub, T_dew, h_bub, h_dew));
        T_HENE_out_2 = fminbnd(objective, T_out_lo, T_out_hi, ...
                               optimset('TolX',1e-3,'Display','off'));
        fprintf('  Max-pinch T_HENE_out = %.3f K\n', T_HENE_out_2);
    else
        T_HENE_out_2 = T_out_hi;
        fprintf('  Both endpoints feasible. Hot-end cap binds.\n');
    end
else
    a = T_out_lo;  b = T_out_hi;
    for it = 1:60
        c = 0.5*(a+b);
        fc = nel_pinch_residual(c, DT_TARGET, Q_HENE_req, ...
            T_HENE_in_2, P_HENE_2, xHe_2, xNe_2, N_cold, ...
            T_R_grid, Q_R_T, T_R_in, N_comp, T_hot_grid, Q_hot_T, ...
            Q_hot_total, T_bub, T_dew, h_bub, h_dew);
        if f_lo * fc < 0
            b = c;
        else
            a = c;  f_lo = fc;
        end
        if abs(b-a) < 1e-3, break; end
    end
    T_HENE_out_2 = 0.5*(a+b);
    fprintf('  Converged in %d iterations: T_HENE_out = %.3f K\n', it, T_HENE_out_2);
end

% --- Rebuild at solved T_HENE_out ---
[T_HENE_grid_2, Q_HENE_T_2, mdot_HENE_2, nfail_2] = build_HENE_safe( ...
    T_HENE_in_2, T_HENE_out_2, P_HENE_2, xHe_2, xNe_2, ...
    Q_HENE_req, N_cold, T_bub, T_dew, h_bub, h_dew);

[T_comp_2, Q_comp_2] = compute_composite( ...
    T_R_grid, Q_R_T, T_HENE_grid_2, Q_HENE_T_2, N_comp);
[dTmin_2, Qp_2, Thp_2, Tcp_2] = pinch_metrics( ...
    Q_hot_T, T_hot_grid, Q_comp_2, T_comp_2);

fprintf('\n  SOLUTION:\n');
fprintf('  T_HENE_out = %.3f K\n', T_HENE_out_2);
fprintf('  mdot_HENE  = %.4f kg/s  (baseline 6.9, ratio %.2f)\n', ...
        mdot_HENE_2, mdot_HENE_2/6.9);
fprintf('  dT_min     = %.3f K  at Q = %.1f kW\n', dTmin_2, Qp_2);
fprintf('  T_hot_p    = %.2f K,  T_cold_p = %.2f K\n', Thp_2, Tcp_2);
fprintf('  REFPROP failures: %d / %d\n\n', nfail_2, N_cold);

%% ---- Comparison table ---------------------------------------------------
fprintf('============================================================\n');
fprintf('  COMPARISON\n');
fprintf('============================================================\n');
fprintf('  %-24s  %14s  %14s\n', '', 'Baseline', 'Nelium-90');
fprintf('  %-24s  %14s  %14s\n', '', '(80/20, 3bar)', '(10/90, 10bar)');
fprintf('  %-24s  %14.4f  %14.4f\n', 'mdot_HENE [kg/s]', mdot_HENE_1, mdot_HENE_2);
fprintf('  %-24s  %14.2f  %14.2f\n', 'mdot ratio [-]', 1.00, mdot_HENE_2/mdot_HENE_1);
fprintf('  %-24s  %14.2f  %14.2f\n', 'T_HENE_in [K]', T_HENE_in_1, T_HENE_in_2);
fprintf('  %-24s  %14.2f  %14.2f\n', 'T_HENE_out [K]', T_HENE_out_1, T_HENE_out_2);
fprintf('  %-24s  %14.3f  %14.3f\n', 'dT_min [K]', dTmin_1, dTmin_2);
fprintf('  %-24s  %14.1f  %14.1f\n', 'Q_pinch [kW]', Qp_1, Qp_2);
if ~isnan(T_bub)
    fprintf('  %-24s  %14s  %14.2f\n', 'T_bubble [K]', 'N/A', T_bub);
    fprintf('  %-24s  %14s  %14.2f\n', 'T_dew [K]',    'N/A', T_dew);
end
fprintf('============================================================\n');

%% ---- Save ---------------------------------------------------------------
write_curves('nelium90_composite_baseline.csv', ...
             Q_hot_T, T_hot_grid, Q_comp_1, T_comp_1);
write_curves('nelium90_composite_nelium90.csv', ...
             Q_hot_T, T_hot_grid, Q_comp_2, T_comp_2);

fid = fopen('nelium90_summary.csv', 'w');
fprintf(fid, 'case,xHe,xNe,P_bar,mdot_HENE_kgps,T_HENE_out_K,dTmin_K,Qpinch_kW,T_bub_K,T_dew_K\n');
fprintf(fid, 'baseline,0.80,0.20,3.0,%.6f,%.4f,%.4f,%.2f,NaN,NaN\n', ...
        mdot_HENE_1, T_HENE_out_1, dTmin_1, Qp_1);
fprintf(fid, 'nelium90,0.10,0.90,10.0,%.6f,%.4f,%.4f,%.2f,%.4f,%.4f\n', ...
        mdot_HENE_2, T_HENE_out_2, dTmin_2, Qp_2, T_bub, T_dew);
fclose(fid);

fprintf('\nSaved: nelium90_composite_baseline.csv\n');
fprintf('       nelium90_composite_nelium90.csv\n');
fprintf('       nelium90_summary.csv\n');
fprintf('Done.\n');


%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================

function tf = test_equilhyd()
    try
        refpropm('H','T',40,'P',10*100,'equilhyd'); %#ok<NASGU>
        tf = true;
    catch
        tf = false;
    end
end

function xp = xp_equilibrium(T)
    T_tab = [ 20  25  30  35  40  45  50  60  70  80  90 100 ...
             120 150 200 250 300 ].';
    x_tab = [0.9989 0.9911 0.9702 0.9316 0.8797 0.8193 0.7558 ...
             0.6309 0.5271 0.4487 0.3913 0.3496 0.2984 0.2638 ...
             0.2507 0.2500 0.2500 ].';
    if T < T_tab(1),      xp = 1.0;
    elseif T > T_tab(end), xp = 0.25;
    else,                  xp = interp1(T_tab, x_tab, T, 'pchip');
    end
end

function [T_grid, h_hot] = build_hot(T_in, T_out, P_bar, N, use_equilhyd)
    T_grid = linspace(T_out, T_in, N).';
    h_hot  = zeros(N,1);
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

function [T_grid, Q_T, mdot, n_fail] = build_HENE_safe( ...
    T_in, T_out, P_bar, xHe, xNe, Q_req_kW, N, ...
    T_bub, T_dew, h_bub, h_dew)
% BUILD_HENE_SAFE  Build He-Ne enthalpy table sized for duty closure.
%
%   For each grid point, attempt REFPROP T,P flash. On failure:
%     - If T is between T_bub and T_dew, interpolate linearly
%       between h_bub and h_dew (first-order two-phase fallback).
%     - Otherwise, mark NaN and fill from neighbors later.

    if nargin < 8,  T_bub = NaN; end
    if nargin < 9,  T_dew = NaN; end
    if nargin < 10, h_bub = NaN; end
    if nargin < 11, h_dew = NaN; end

    T_grid = linspace(T_in, T_out, N).';
    h      = NaN(N, 1);
    P_kPa  = P_bar * 100;
    n_fail = 0;
    have_vle = ~isnan(T_bub) && ~isnan(T_dew) && ...
               ~isnan(h_bub) && ~isnan(h_dew);

    for k = 1:N
        try
            h(k) = refpropm('H','T',T_grid(k),'P',P_kPa, ...
                            'helium','neon',[xHe, xNe]);
        catch
            n_fail = n_fail + 1;
            if have_vle && T_grid(k) >= T_bub && T_grid(k) <= T_dew
                frac = (T_grid(k) - T_bub) / (T_dew - T_bub);
                h(k) = h_bub + frac * (h_dew - h_bub);
            end
        end
    end

    good = ~isnan(h);
    if any(~good) && sum(good) >= 2
        h = interp1(T_grid(good), h(good), T_grid, 'pchip', 'extrap');
    end

    dh   = h(end) - h(1);
    mdot = Q_req_kW * 1000 / dh;
    Q_T  = mdot * (h - h(1)) / 1000;
end

function f = nel_pinch_residual(T_HENE_out, DT_target, Q_HENE_req, ...
    T_HENE_in, P_HENE, xHe, xNe, N_cold, T_R_grid, Q_R_T, T_R_in, ...
    N_comp, T_hot_grid, Q_hot_T, Q_hot_total, ...
    T_bub, T_dew, h_bub, h_dew)
% NEL_PINCH_RESIDUAL  f = min_dT - DT_target, used by bisection solver.
%   Builds HENE at trial T_HENE_out with two-phase fallback, computes
%   cold composite, and returns the residual for the pinch target.

    [T_HENE_grid, Q_HENE_T, ~, ~] = build_HENE_safe( ...
        T_HENE_in, T_HENE_out, P_HENE, xHe, xNe, Q_HENE_req, N_cold, ...
        T_bub, T_dew, h_bub, h_dew);

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

function [T_comp, Q_comp] = compute_composite(T_R, Q_R, T_HENE, Q_HENE, N)
    T_lo = min(T_R(1),   T_HENE(1));
    T_hi = max(T_R(end), T_HENE(end));
    T_comp = linspace(T_lo, T_hi, N).';
    Q_R_on    = interp_clipped(T_R,    Q_R,    T_comp);
    Q_HENE_on = interp_clipped(T_HENE, Q_HENE, T_comp);
    Q_comp    = cummax(Q_R_on + Q_HENE_on);
end

function Qout = interp_clipped(T_in, Q_in, T_query)
    Qout = zeros(size(T_query));
    inside = T_query >= T_in(1) & T_query <= T_in(end);
    Qout(inside) = interp1(T_in, Q_in, T_query(inside), 'pchip');
    Qout(T_query > T_in(end)) = Q_in(end);
end

function [dTmin, Qp, Thp, Tcp] = pinch_metrics(Q_hot, T_hot, Q_cold, T_cold)
    Q_max = min(Q_hot(end), Q_cold(end));
    Qg    = linspace(0, Q_max, 4001).';
    Th    = interp1(Q_hot,  T_hot,  Qg, 'pchip');
    Tc    = interp1(Q_cold, T_cold, Qg, 'pchip');
    dT    = Th - Tc;
    [dTmin, idx] = min(dT);
    Qp  = Qg(idx);
    Thp = Th(idx);
    Tcp = Tc(idx);
end

function write_curves(filename, Q_hot, T_hot, Q_cold, T_cold)
    n_hot  = length(Q_hot);
    n_cold = length(Q_cold);
    n_max  = max(n_hot, n_cold);
    Q_h = [Q_hot;  NaN(n_max-n_hot,  1)];
    T_h = [T_hot;  NaN(n_max-n_hot,  1)];
    Q_c = [Q_cold; NaN(n_max-n_cold, 1)];
    T_c = [T_cold; NaN(n_max-n_cold, 1)];
    T = table(Q_h, T_h, Q_c, T_c, ...
        'VariableNames', {'Q_hot_kW','T_hot_K','Q_cold_kW','T_cold_K'});
    writetable(T, filename);
end