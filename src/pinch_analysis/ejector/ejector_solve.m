function results = ejector_solve(f, eta, opts)
% EJECTOR_SOLVE  Solve the ejector model for given f and eta.
%
%   results = EJECTOR_SOLVE(f, eta) solves with default operating conditions.
%   results = EJECTOR_SOLVE(f, eta, opts) overrides defaults via opts struct.
%
%   Inputs:
%     f   - vapor recirculation fraction [-]
%     eta - ejector efficiency [-]
%     opts (optional) struct with any of the following fields:
%       Pprimary      - primary stream pressure [kPa]   (default 7500)
%       Tprimary      - primary stream temperature [K]  (default 74.15)
%       mp            - product mass flow [kg/s]        (default 1.4847)
%       TBOG          - BOG temperature [K]             (default 21.7)
%       PBOG_init     - BOG pressure before throttle    (default 150 kPa)
%       BOR           - boil-off ratio mBOG/mL [-]      (default 0.0233)
%       Psep          - separator pressure [kPa]        (default 100)
%       T_postPFHX5   - PFHX-5 outlet temp [K]          (default 27.95)
%       fluid         - REFPROP fluid name              (default 'parahydrogen')
%       tol_inner     - Pd bisection tolerance on eta   (default 1e-5)
%       maxiter_inner - max bisection iterations        (default 200)
%       tol_outer     - x fixed-point tolerance         (default 1e-4)
%       maxiter_outer - max outer iterations            (default 10)
%       verbose       - print results to console        (default true)
%
%   Output:
%     results - struct with all inputs echoed plus:
%       converged          - logical
%       failure_reason     - string ('none' on success, otherwise diagnostic)
%       x                  - separator quality [-]
%       D                  - D-factor [-]
%       entrainment_ratio  - ms/mp [-]
%       Pd_kPa, Pd_bar     - diffuser pressure
%       pressure_penalty_bar
%       mp, mDIF, ms, mBOG, mE, mv, mR, mL  - mass flows [kg/s]
%       Tdiffusor          - diffuser outlet temperature [K]
%       xdiffusor          - diffuser outlet quality (-998 if single-phase)
%       hdiffusor          - diffuser outlet enthalpy [J/kg]
%       Q_PFHX5_kW         - PFHX-5 duty [kW]
%       iter_inner         - cumulative inner bisection iterations
%       iter_outer         - outer fixed-point iterations
%
%   On any failure, all numeric output fields remain NaN, converged is
%   false, and failure_reason takes one of:
%     'refprop_nan', 'x_out_of_dome', 'D_nonphysical',
%     'bisection_no_converge', 'fixed_point_no_converge'

    % =================================================================
    %  Apply default options
    % =================================================================
    if nargin < 3, opts = struct(); end
    defaults = struct( ...
        'Pprimary',      75 * 100, ...
        'Tprimary',      74.15, ...
        'mp',            1.4847, ...
        'TBOG',          21.7, ...
        'PBOG_init',     1.5 * 100, ...
        'BOR',           0.0233, ...
        'Psep',          100, ...
        'T_postPFHX5',   27.95, ...
        'fluid',         'parahydrogen', ...
        'tol_inner',     1e-5, ...
        'maxiter_inner', 200, ...
        'tol_outer',     1e-4, ...
        'maxiter_outer', 10, ...
        'verbose',       true);
    fns = fieldnames(defaults);
    for k = 1:numel(fns)
        if ~isfield(opts, fns{k})
            opts.(fns{k}) = defaults.(fns{k});
        end
    end

    % =================================================================
    %  Initialize NaN-filled results struct
    % =================================================================
    results = init_results(f, eta, opts);

    % =================================================================
    %  Precompute fixed stream properties
    % =================================================================
    try
        hprimary = refpropm('h', 'T', opts.Tprimary, 'P', opts.Pprimary, opts.fluid);
        sprimary = refpropm('s', 'T', opts.Tprimary, 'P', opts.Pprimary, opts.fluid);
        hBOG     = refpropm('h', 'T', opts.TBOG,     'P', opts.PBOG_init, opts.fluid);
        hL_sat   = refpropm('h', 'P', opts.Psep, 'Q', 0, opts.fluid);
        hV_sat   = refpropm('h', 'P', opts.Psep, 'Q', 1, opts.fluid);
        TV_sat   = refpropm('T', 'P', opts.Psep, 'Q', 1, opts.fluid);
    catch ME
        results = fail(results, 'refprop_nan', opts, ME.message);
        return
    end

    if any(isnan([hprimary, sprimary, hBOG, hL_sat, hV_sat, TV_sat]))
        results = fail(results, 'refprop_nan', opts, 'NaN in fixed properties');
        return
    end

    hrecirc  = hV_sat;
    Psuction = opts.Psep;

    % =================================================================
    %  Initial x estimate using hPFHX5 at reference (motive) pressure
    % =================================================================
    try
        hPFHX5_ref = refpropm('h', 'T', opts.T_postPFHX5, 'P', opts.Pprimary, opts.fluid);
    catch ME
        results = fail(results, 'refprop_nan', opts, ME.message);
        return
    end
    x = (hPFHX5_ref - hL_sat) / (hV_sat - hL_sat);

    if isnan(x) || x < 0 || x > 1
        results.x = x;
        results = fail(results, 'x_out_of_dome', opts, ...
                       sprintf('Initial x = %.4f outside [0,1]', x));
        return
    end

    % =================================================================
    %  Outer fixed-point loop:  x  <->  Pd
    % =================================================================
    pdiffusor        = NaN;
    iter_inner_total = 0;
    fp_converged     = false;
    iter_outer       = 0;
    hPFHX5           = hPFHX5_ref;
    mb               = struct();
    hsecondary       = NaN;
    x_old            = x;     % defensive init for post-loop diagnostic

    for iter_outer = 1:opts.maxiter_outer
        x_old = x;

        % --- Mass balance from current x ---
        D = 1 / ((1 - opts.BOR) + x*(opts.BOR - f));
        if isnan(D) || isinf(D) || D < 0
            results.D = D;
            results = fail(results, 'D_nonphysical', opts, ...
                           sprintf('D = %.4f at x = %.4f, f = %.4f', D, x, f));
            return
        end

        mb = mass_balance(x, f, opts);

        % --- Mixed secondary properties at suction ---
        if mb.mE < 1e-12
            hsecondary = hBOG;
        else
            hsecondary = (mb.mBOG * hBOG + mb.mE * hrecirc) / mb.ms;
        end

        try
            ssecondary = refpropm('s', 'P', Psuction, 'H', hsecondary, opts.fluid);
        catch ME
            results = fail(results, 'refprop_nan', opts, ME.message);
            return
        end
        if isnan(hsecondary) || isnan(ssecondary)
            results = fail(results, 'refprop_nan', opts, 'NaN in secondary props');
            return
        end

        % --- Bisection on Pd to satisfy efficiency constraint ---
        massratio = mb.ms / opts.mp;
        pmin      = Psuction;
        pmax      = opts.Pprimary;
        iter_in   = 0;
        bis_conv  = false;
        etaCalc   = NaN;     % so post-loop diagnostic prints cleanly even if
                             % every iteration hits a continue guard

        while iter_in < opts.maxiter_inner
            iter_in   = iter_in + 1;
            pdiffusor = (pmin + pmax) / 2;

            try
                h_iss = refpropm('h', 'P', pdiffusor, 'S', ssecondary, opts.fluid);
                h_isp = refpropm('h', 'P', pdiffusor, 'S', sprimary,   opts.fluid);
            catch
                pmin = pmin + 0.01 * (pmax - pmin);
                continue
            end

            denom = hprimary - h_isp;
            if abs(denom) < 1e-10 || isnan(h_iss) || isnan(h_isp)
                pmin = pmin + 0.01 * (pmax - pmin);
                continue
            end

            etaCalc   = massratio * (h_iss - hsecondary) / denom;
            err_inner = etaCalc - eta;

            if abs(err_inner) <= opts.tol_inner
                bis_conv = true;
                break
            end

            if err_inner > 0
                pmax = pdiffusor;
            else
                pmin = pdiffusor;
            end

            if (pmax - pmin) < 1e-8
                break
            end
        end

        iter_inner_total = iter_inner_total + iter_in;

        if ~bis_conv
            results = fail(results, 'bisection_no_converge', opts, ...
                           sprintf('eta_calc = %.6f, eta = %.6f, Pd = %.2f kPa', ...
                                   etaCalc, eta, pdiffusor));
            return
        end

        % --- Update x from actual hPFHX5 at solved Pd ---
        try
            hPFHX5 = refpropm('h', 'T', opts.T_postPFHX5, 'P', pdiffusor, opts.fluid);
        catch ME
            results = fail(results, 'refprop_nan', opts, ME.message);
            return
        end
        x_new = (hPFHX5 - hL_sat) / (hV_sat - hL_sat);

        if isnan(x_new) || x_new < 0 || x_new > 1
            results.x = x_new;
            results = fail(results, 'x_out_of_dome', opts, ...
                           sprintf('Refined x = %.4f outside [0,1]', x_new));
            return
        end

        x = x_new;

        if abs(x - x_old) < opts.tol_outer
            fp_converged = true;
            break
        end
    end

    if ~fp_converged
        results = fail(results, 'fixed_point_no_converge', opts, ...
                       sprintf('|dx| = %.2e after %d iters', ...
                               abs(x - x_old), iter_outer));
        return
    end

    % =================================================================
    %  Final self-consistent state at converged x
    % =================================================================
    mb = mass_balance(x, f, opts);

    if mb.mE < 1e-12
        hsecondary = hBOG;
    else
        hsecondary = (mb.mBOG * hBOG + mb.mE * hrecirc) / mb.ms;
    end

    hdiffusor = (opts.mp * hprimary + mb.ms * hsecondary) / mb.mDIF;

    try
        Tdiffusor = refpropm('T', 'P', pdiffusor, 'H', hdiffusor, opts.fluid);
        xdiffusor = refpropm('q', 'P', pdiffusor, 'H', hdiffusor, opts.fluid);
    catch ME
        results = fail(results, 'refprop_nan', opts, ME.message);
        return
    end

    Q_PFHX5 = mb.mDIF * (hdiffusor - hPFHX5);

    % =================================================================
    %  Pack results
    % =================================================================
    results.converged            = true;
    results.failure_reason       = 'none';
    results.x                    = x;
    results.D                    = mb.D;
    results.entrainment_ratio    = mb.ms / opts.mp;
    results.Pd_kPa               = pdiffusor;
    results.Pd_bar               = pdiffusor / 100;
    results.pressure_penalty_bar = (opts.Pprimary - pdiffusor) / 100;
    results.mp                   = opts.mp;
    results.mDIF                 = mb.mDIF;
    results.ms                   = mb.ms;
    results.mBOG                 = mb.mBOG;
    results.mE                   = mb.mE;
    results.mv                   = mb.mv;
    results.mR                   = mb.mR;
    results.mL                   = mb.mL;
    results.Tdiffusor            = Tdiffusor;
    results.xdiffusor            = xdiffusor;
    results.hdiffusor            = hdiffusor;
    results.Q_PFHX5_kW           = Q_PFHX5 / 1000;
    results.iter_inner           = iter_inner_total;
    results.iter_outer           = iter_outer;

    if opts.verbose
        print_results(results, opts);
    end
end


% =====================================================================
%  Local helpers
% =====================================================================

function results = init_results(f, eta, opts)
% Build a NaN-filled results struct with consistent schema.
    results = struct( ...
        'f',                    f, ...
        'eta',                  eta, ...
        'BOR',                  opts.BOR, ...
        'mp_input',             opts.mp, ...
        'T_postPFHX5',          opts.T_postPFHX5, ...
        'Pprimary_kPa',         opts.Pprimary, ...
        'Tprimary',             opts.Tprimary, ...
        'Psep_kPa',             opts.Psep, ...
        'converged',            false, ...
        'failure_reason',       'unknown', ...
        'x',                    NaN, ...
        'D',                    NaN, ...
        'entrainment_ratio',    NaN, ...
        'Pd_kPa',               NaN, ...
        'Pd_bar',               NaN, ...
        'pressure_penalty_bar', NaN, ...
        'mp',                   NaN, ...
        'mDIF',                 NaN, ...
        'ms',                   NaN, ...
        'mBOG',                 NaN, ...
        'mE',                   NaN, ...
        'mv',                   NaN, ...
        'mR',                   NaN, ...
        'mL',                   NaN, ...
        'Tdiffusor',            NaN, ...
        'xdiffusor',            NaN, ...
        'hdiffusor',            NaN, ...
        'Q_PFHX5_kW',           NaN, ...
        'iter_inner',           0, ...
        'iter_outer',           0);
end


function mb = mass_balance(x, f, opts)
% Closed-form D-factor mass balance.
    D = 1 / ((1 - opts.BOR) + x*(opts.BOR - f));
    mb.D    = D;
    mb.mDIF = D * opts.mp;
    mb.mv   = x * D * opts.mp;
    mb.mR   = (1-f) * x * D * opts.mp;
    mb.mE   = f     * x * D * opts.mp;
    mb.mBOG = opts.BOR * (1-x) * D * opts.mp;
    mb.mL   = (1-x) * D * opts.mp;
    mb.ms   = (opts.BOR*(1-x) + f*x) * D * opts.mp;
end


function r = fail(r, reason, opts, detail)
% Mark results struct as failed and optionally print a verbose message.
    r.failure_reason = reason;
    r.converged      = false;
    if opts.verbose
        if nargin < 4 || isempty(detail)
            fprintf('Solver failed: %s\n', reason);
        else
            fprintf('Solver failed: %s  (%s)\n', reason, detail);
        end
    end
end


function print_results(r, opts)
% Verbose console report (success path only).
    fprintf('=== Ejector Results ===\n');
    fprintf('  Status                  : CONVERGED\n');
    fprintf('  BOR = %.4f, f = %.3f, eta = %.3f\n', opts.BOR, r.f, r.eta);
    fprintf('  Outer fixed-point iters : %d\n', r.iter_outer);
    fprintf('  Inner bisection iters   : %d (cumulative)\n', r.iter_inner);
    fprintf('  Separator quality x     : %.6f\n', r.x);
    fprintf('  D-factor                : %.6f\n', r.D);
    fprintf('\n');

    fprintf('  --- Ejector ---\n');
    fprintf('  Diffuser pressure       : %.2f kPa  (%.3f bar)\n', r.Pd_kPa, r.Pd_bar);
    fprintf('  Pressure penalty        : %.3f bar\n', r.pressure_penalty_bar);
    fprintf('  Diffuser temperature    : %.2f K\n', r.Tdiffusor);
    if r.xdiffusor >= 0 && r.xdiffusor <= 1
        fprintf('  Diffuser quality        : %.4f\n', r.xdiffusor);
    else
        fprintf('  Diffuser quality        : single-phase\n');
    end
    fprintf('  Entrainment ratio       : %.6f  (ms/mp)\n', r.entrainment_ratio);
    fprintf('\n');

    fprintf('  --- PFHX-5 ---\n');
    fprintf('  Inlet temperature       : %.2f K  (ejector discharge)\n', r.Tdiffusor);
    fprintf('  Outlet temperature      : %.2f K  (fixed input)\n', opts.T_postPFHX5);
    fprintf('  Pressure                : %.3f bar\n', r.Pd_bar);
    fprintf('  Duty                    : %.2f kW\n', r.Q_PFHX5_kW);
    fprintf('\n');

    fprintf('  --- Mass flows ---\n');
    fprintf('  Product (mp)            : %.4f kg/s\n', r.mp);
    fprintf('  Diffuser (mDIF)         : %.4f kg/s\n', r.mDIF);
    fprintf('  Secondary (ms)          : %.4f kg/s\n', r.ms);
    fprintf('    BOG (mBOG)            : %.4f kg/s  (BOR * mL)\n', r.mBOG);
    fprintf('    Recirc (mE)           : %.4f kg/s  (f * mv)\n', r.mE);
    fprintf('  Vapor (mv)              : %.4f kg/s\n', r.mv);
    fprintf('    To ejector (mE)       : %.4f kg/s\n', r.mE);
    fprintf('    To compressor (mR)    : %.4f kg/s\n', r.mR);
    fprintf('  Liquid (mL)             : %.4f kg/s\n', r.mL);
    fprintf('\n');

    fprintf('  --- Sanity checks ---\n');
    fprintf('  mBOG/mL                 : %.4f  (should be BOR = %.4f)\n', ...
            r.mBOG/r.mL, opts.BOR);
    fprintf('  mL - mBOG + mR          : %.4f  (should be mp = %.4f)\n', ...
            r.mL - r.mBOG + r.mR, r.mp);
    fprintf('  mp + ms                 : %.4f  (should be mDIF = %.4f)\n', ...
            r.mp + r.ms, r.mDIF);
    fprintf('\n');
end
