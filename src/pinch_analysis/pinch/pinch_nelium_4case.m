% pinch_nelium_4case.m
% -------------------------------------------------------------------
% Compares PFHX-5 composite curves across four cases:
%
%   Case A: H2 at 75 bar + He-Ne 80/20 at 3 bar  (Moro baseline)
%   Case B: H2 at 75 bar + He-Ne 10/90 at 10 bar (Nelium-90, Wilhelmsen)
%   Case C: H2 at 21 bar + He-Ne 80/20 at 3 bar  (ejector, Pp=40)
%   Case D: H2 at 21 bar + He-Ne 10/90 at 10 bar (ejector + Nelium)
%
%   Case A uses fixed Moro streams (mdot_HENE=6.9, T_HENE_out=74.05).
%   Cases B-D use bisection on T_HENE_out for dT_min = 2.48 K.
%
% Output: nelium_case_A.csv through nelium_case_D.csv
%         nelium_4case_summary.csv
% -------------------------------------------------------------------

clear; clc; close all;

%% ---- Configuration ----------------------------------------------------
N_hot  = 401;
N_cold = 401;
N_comp = 801;
DT_TARGET = 2.48;
P_R = 1.0;
T_R_in  = 20.40;
T_R_out = 46.65;
xp_R = 0.999;

use_equilhyd = test_equilhyd();

%% ---- H2/R-stream configurations ----------------------------------------
cfg = struct();

cfg(1).label     = '75 bar';
cfg(1).P_hot     = 75.0;
cfg(1).T_hot_in  = 74.15;
cfg(1).T_hot_out = 27.95;
cfg(1).mdot_H    = 1.4847;
cfg(1).mdot_R    = 0.4847;
cfg(1).T_HENE_in = 27.45;

cfg(2).label     = '21 bar';
cfg(2).P_hot     = 21.0;
cfg(2).T_hot_in  = 68.12;
cfg(2).T_hot_out = 27.95;
cfg(2).mdot_H    = 1.5878;
cfg(2).mdot_R    = 0.2994;
cfg(2).T_HENE_in = 26.20;

%% ---- Precompute R-stream enthalpy --------------------------------------
fprintf('Computing R-stream enthalpy...\n');
T_R_grid = linspace(T_R_in, T_R_out, N_cold).';
h_R = zeros(N_cold, 1);
for i = 1:N_cold
    h_R(i) = refpropm('H','T',T_R_grid(i),'P',P_R*100, ...
                      'orthohyd','parahyd',[1-xp_R, xp_R]);
end

%% ---- VLE boundaries for Nelium-90 at 10 bar ----------------------------
P_HENE_NEL = 10.0;   % Nelium cases at 10 bar (Wilhelmsen conditions)
xHe_nel = 0.10;  xNe_nel = 0.90;
T_bub = NaN;  T_dew = NaN;  h_bub = NaN;  h_dew = NaN;
try
    T_bub = refpropm('T','P',P_HENE_NEL*100,'Q',0,'helium','neon',[xHe_nel,xNe_nel]);
    T_dew = refpropm('T','P',P_HENE_NEL*100,'Q',1,'helium','neon',[xHe_nel,xNe_nel]);
    h_bub = refpropm('H','P',P_HENE_NEL*100,'Q',0,'helium','neon',[xHe_nel,xNe_nel]);
    h_dew = refpropm('H','P',P_HENE_NEL*100,'Q',1,'helium','neon',[xHe_nel,xNe_nel]);
    fprintf('VLE for Nelium-90 at %.0f bar:\n', P_HENE_NEL);
    fprintf('  T_bubble = %.2f K,  T_dew = %.2f K,  Glide = %.2f K\n\n', ...
            T_bub, T_dew, T_dew - T_bub);
catch ME
    fprintf('WARNING: VLE failed: %s\n\n', ME.message);
end

%% ---- Case definitions ---------------------------------------------------
cases = struct();

% Case A: 75 bar H2, 80/20 at 3 bar, fixed baseline
cases(1).tag    = 'A';
cases(1).label  = 'H2 75 bar, He-Ne 80/20 (3 bar)';
cases(1).cfg    = 1;
cases(1).xHe    = 0.80;
cases(1).xNe    = 0.20;
cases(1).P_HENE = 3.0;
cases(1).mode   = 'fixed';

% Case B: 75 bar H2, 10/90 at 10 bar, solver
cases(2).tag    = 'B';
cases(2).label  = 'H2 75 bar, He-Ne 10/90 (10 bar)';
cases(2).cfg    = 1;
cases(2).xHe    = 0.10;
cases(2).xNe    = 0.90;
cases(2).P_HENE = P_HENE_NEL;
cases(2).mode   = 'solver';

% Case C: 21 bar H2, 80/20 at 3 bar, solver
cases(3).tag    = 'C';
cases(3).label  = 'H2 21 bar, He-Ne 80/20 (3 bar)';
cases(3).cfg    = 2;
cases(3).xHe    = 0.80;
cases(3).xNe    = 0.20;
cases(3).P_HENE = 3.0;
cases(3).mode   = 'solver';

% Case D: 21 bar H2, 10/90 at 10 bar, solver
cases(4).tag    = 'D';
cases(4).label  = 'H2 21 bar, He-Ne 10/90 (10 bar)';
cases(4).cfg    = 2;
cases(4).xHe    = 0.10;
cases(4).xNe    = 0.90;
cases(4).P_HENE = P_HENE_NEL;
cases(4).mode   = 'solver';

%% ---- Run each case -----------------------------------------------------
Ncases = length(cases);
for ic = 1:Ncases
    c  = cases(ic);
    cf = cfg(c.cfg);
    P_HENE = c.P_HENE;

    fprintf('============================================================\n');
    fprintf('  CASE %s: %s\n', c.tag, c.label);
    fprintf('============================================================\n');

    % --- Hot side ---
    [T_hot, h_hot] = build_hot(cf.T_hot_in, cf.T_hot_out, cf.P_hot, ...
                               N_hot, use_equilhyd);
    Q_hot = cf.mdot_H * (h_hot - h_hot(1)) / 1000;
    Q_hot_total = Q_hot(end);

    % --- R stream ---
    Q_R = cf.mdot_R * (h_R - h_R(1)) / 1000;
    Q_R_total = Q_R(end);

    Q_HENE_req = Q_hot_total - Q_R_total;
    fprintf('  Q_hot = %.2f kW,  Q_R = %.2f kW,  Q_HENE_req = %.2f kW\n', ...
            Q_hot_total, Q_R_total, Q_HENE_req);

    % --- VLE info (only for Nelium cases at 10 bar) ---
    if c.xNe > 0.5
        vle = struct('T_bub',T_bub,'T_dew',T_dew,'h_bub',h_bub,'h_dew',h_dew);
    else
        vle = struct('T_bub',NaN,'T_dew',NaN,'h_bub',NaN,'h_dew',NaN);
    end

    % --- Build HENE and composite ---
    if strcmp(c.mode, 'fixed')
        % Case A: fixed Moro baseline
        T_HENE_out = 74.05;
        mdot_HENE  = 6.9;
        T_HENE_grid = linspace(cf.T_HENE_in, T_HENE_out, N_cold).';
        h_HENE = zeros(N_cold, 1);
        nfail = 0;
        for k = 1:N_cold
            try
                h_HENE(k) = refpropm('H','T',T_HENE_grid(k),'P',P_HENE*100, ...
                                     'helium','neon',[c.xHe, c.xNe]);
            catch
                h_HENE(k) = NaN;
                nfail = nfail + 1;
            end
        end
        good = ~isnan(h_HENE);
        if any(~good) && sum(good) >= 2
            h_HENE = interp1(T_HENE_grid(good), h_HENE(good), T_HENE_grid, 'pchip', 'extrap');
        end
        Q_HENE = mdot_HENE * (h_HENE - h_HENE(1)) / 1000;
        fprintf('  Mode: fixed  mdot=%.2f kg/s  T_out=%.2f K\n', mdot_HENE, T_HENE_out);

    else
        % Solver mode: bisect on T_HENE_out for dT_min = DT_TARGET
        T_out_lo = cf.T_HENE_in + 0.5;
        T_out_hi = cf.T_hot_in - DT_TARGET;

        f_lo = pinch_residual(T_out_lo, DT_TARGET, Q_HENE_req, ...
            cf.T_HENE_in, P_HENE, c.xHe, c.xNe, N_cold, ...
            T_R_grid, Q_R, T_R_in, N_comp, T_hot, Q_hot, Q_hot_total, vle);
        f_hi = pinch_residual(T_out_hi, DT_TARGET, Q_HENE_req, ...
            cf.T_HENE_in, P_HENE, c.xHe, c.xNe, N_cold, ...
            T_R_grid, Q_R, T_R_in, N_comp, T_hot, Q_hot, Q_hot_total, vle);

        fprintf('  Bisecting T_HENE_out in [%.2f, %.2f] K\n', T_out_lo, T_out_hi);
        fprintf('  f(lo) = %+.3f,  f(hi) = %+.3f\n', f_lo, f_hi);

        if f_lo * f_hi > 0
            if f_lo < 0 && f_hi < 0
                fprintf('  Target infeasible. Searching max-pinch...\n');
                obj = @(T) -(pinch_residual(T, 0, Q_HENE_req, ...
                    cf.T_HENE_in, P_HENE, c.xHe, c.xNe, N_cold, ...
                    T_R_grid, Q_R, T_R_in, N_comp, T_hot, Q_hot, ...
                    Q_hot_total, vle));
                T_HENE_out = fminbnd(obj, T_out_lo, T_out_hi, ...
                                     optimset('TolX',1e-3,'Display','off'));
            else
                T_HENE_out = T_out_hi;
            end
        else
            a = T_out_lo;  b = T_out_hi;
            for it = 1:60
                mid = 0.5*(a+b);
                fm = pinch_residual(mid, DT_TARGET, Q_HENE_req, ...
                    cf.T_HENE_in, P_HENE, c.xHe, c.xNe, N_cold, ...
                    T_R_grid, Q_R, T_R_in, N_comp, T_hot, Q_hot, ...
                    Q_hot_total, vle);
                if f_lo * fm < 0
                    b = mid;
                else
                    a = mid; f_lo = fm;
                end
                if abs(b-a) < 1e-3, break; end
            end
            T_HENE_out = 0.5*(a+b);
            fprintf('  Converged in %d iterations\n', it);
        end

        [T_HENE_grid, Q_HENE, mdot_HENE, nfail] = build_HENE_safe( ...
            cf.T_HENE_in, T_HENE_out, P_HENE, c.xHe, c.xNe, ...
            Q_HENE_req, N_cold, vle);
        fprintf('  T_HENE_out = %.3f K  mdot = %.4f kg/s\n', T_HENE_out, mdot_HENE);
    end

    % --- Composite and pinch ---
    [T_comp, Q_comp] = compute_composite(T_R_grid, Q_R, T_HENE_grid, Q_HENE, N_comp);
    [dTmin, Qp, Thp, Tcp] = pinch_metrics(Q_hot, T_hot, Q_comp, T_comp);

    fprintf('  dT_min = %.3f K  at Q = %.1f kW\n', dTmin, Qp);
    fprintf('  REFPROP failures: %d / %d\n\n', nfail, N_cold);

    % --- Store results ---
    cases(ic).Q_hot      = Q_hot;
    cases(ic).T_hot      = T_hot;
    cases(ic).Q_comp     = Q_comp;
    cases(ic).T_comp     = T_comp;
    cases(ic).dTmin      = dTmin;
    cases(ic).Qpinch     = Qp;
    cases(ic).mdot_HENE  = mdot_HENE;
    cases(ic).T_HENE_out = T_HENE_out;
    cases(ic).P_hot_val  = cf.P_hot;
    cases(ic).T_HENE_in  = cf.T_HENE_in;
    cases(ic).nfail      = nfail;

    % --- Save CSV ---
    fname = sprintf('nelium_case_%s.csv', c.tag);
    write_curves(fname, Q_hot, T_hot, Q_comp, T_comp);
    fprintf('  Saved %s\n\n', fname);
end

%% ---- Summary table ------------------------------------------------------
fprintf('============================================================\n');
fprintf('  SUMMARY\n');
fprintf('============================================================\n');
fprintf('  %-6s  %-8s  %-6s  %6s  %8s  %8s  %8s  %8s\n', ...
        'Case', 'P_H2', 'xHe', 'P_HN', 'mdot_HN', 'T_HN_out', 'dTmin', 'Q_pinch');
for ic = 1:Ncases
    fprintf('  %-6s  %5.0f bar  %5.2f  %4.0f bar  %8.3f  %8.2f  %8.3f  %8.1f\n', ...
        cases(ic).tag, cases(ic).P_hot_val, cases(ic).xHe, cases(ic).P_HENE, ...
        cases(ic).mdot_HENE, cases(ic).T_HENE_out, ...
        cases(ic).dTmin, cases(ic).Qpinch);
end
fprintf('============================================================\n');

%% ---- Summary CSV --------------------------------------------------------
fid = fopen('nelium_4case_summary.csv', 'w');
fprintf(fid, 'case,P_H2_bar,xHe,xNe,P_HENE_bar,T_HENE_in_K,T_HENE_out_K,mdot_HENE_kgps,dTmin_K,Qpinch_kW,T_bub_K,T_dew_K,nfail\n');
for ic = 1:Ncases
    c = cases(ic);
    if c.xNe > 0.5
        tb = T_bub; td = T_dew;
    else
        tb = NaN; td = NaN;
    end
    fprintf(fid, '%s,%.1f,%.2f,%.2f,%.1f,%.4f,%.4f,%.6f,%.4f,%.2f,%.4f,%.4f,%d\n', ...
        c.tag, c.P_hot_val, c.xHe, c.xNe, c.P_HENE, c.T_HENE_in, ...
        c.T_HENE_out, c.mdot_HENE, c.dTmin, c.Qpinch, tb, td, c.nfail);
end
fclose(fid);
fprintf('\nSaved nelium_4case_summary.csv\n');
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
    T_tab = [20 25 30 35 40 45 50 60 70 80 90 100 120 150 200 250 300].';
    x_tab = [0.9989 0.9911 0.9702 0.9316 0.8797 0.8193 0.7558 ...
             0.6309 0.5271 0.4487 0.3913 0.3496 0.2984 0.2638 ...
             0.2507 0.2500 0.2500].';
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
    T_in, T_out, P_bar, xHe, xNe, Q_req_kW, N, vle)

    T_grid = linspace(T_in, T_out, N).';
    h      = NaN(N, 1);
    P_kPa  = P_bar * 100;
    n_fail = 0;
    have_vle = ~isnan(vle.T_bub) && ~isnan(vle.T_dew) && ...
               ~isnan(vle.h_bub) && ~isnan(vle.h_dew);

    for k = 1:N
        try
            h(k) = refpropm('H','T',T_grid(k),'P',P_kPa, ...
                            'helium','neon',[xHe, xNe]);
        catch
            n_fail = n_fail + 1;
            if have_vle && T_grid(k) >= vle.T_bub && T_grid(k) <= vle.T_dew
                frac = (T_grid(k) - vle.T_bub) / (vle.T_dew - vle.T_bub);
                h(k) = vle.h_bub + frac * (vle.h_dew - vle.h_bub);
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

function f = pinch_residual(T_HENE_out, DT_target, Q_HENE_req, ...
    T_HENE_in, P_HENE, xHe, xNe, N_cold, T_R_grid, Q_R_T, T_R_in, ...
    N_comp, T_hot_grid, Q_hot_T, Q_hot_total, vle)

    [T_HN, Q_HN, ~, ~] = build_HENE_safe( ...
        T_HENE_in, T_HENE_out, P_HENE, xHe, xNe, Q_HENE_req, N_cold, vle);

    T_lo = min(T_R_in, T_HENE_in);
    T_hi = max(T_R_grid(end), T_HENE_out);
    T_comp = linspace(T_lo, T_hi, N_comp).';
    Q_R_on  = interp_clipped(T_R_grid, Q_R_T, T_comp);
    Q_HN_on = interp_clipped(T_HN,    Q_HN,  T_comp);
    Q_cold  = cummax(Q_R_on + Q_HN_on);

    Q_max = min(Q_hot_total, Q_cold(end));
    Qg  = linspace(0, Q_max, 2001).';
    T_h = interp1(Q_hot_T, T_hot_grid, Qg, 'pchip');
    T_c = interp1(Q_cold,  T_comp,     Qg, 'pchip');
    f   = min(T_h - T_c) - DT_target;
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