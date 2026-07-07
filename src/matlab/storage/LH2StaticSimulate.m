function data = LH2StaticSimulate(name)
% LH2StaticSimulate  Single-vessel LH2 static boil-off simulation.
%
%   data = LH2StaticSimulate(name)
%
%   KEY DIFFERENCE from Petitpas original:
%   Vapor state variable is TEMPERATURE, not internal energy. This
%   eliminates all (D,U) and (P,U) REFPROP flashes, which are fragile
%   near the saturation dome. Every REFPROP call is now (T,D), (T,P),
%   or (T,Q) -- all robust.
%
%   Physics are identical to Petitpas (2018):
%     - Non-equilibrium condensation-evaporation at the saturated film
%     - Wall thermal mass with temperature-dependent Cp (SS304)
%     - Real gas EOS via REFPROP (para-hydrogen)
%     - Vent valve (PRD) with hysteresis and event detection
%
%   Rayleigh number bug from original code is FIXED (nu -> nu^2).

% =========================================================================
%  1. Load parameters
% =========================================================================
try
    P = evalin('base','LH2Model');
catch
    evalin('base','inputs_StaticBoiloff');
    P = evalin('base','LH2Model');
end

if nargin < 1, name = 'Static boil-off'; end

P.waitbar = waitbar(0, ['Simulating ' name '...']);

% =========================================================================
%  2. Build initial state vector
% =========================================================================
% State vector layout:
%   x(1)                     = mL     (liquid mass)          [kg]
%   x(2 : 1+nL)             = uL     (liquid internal energies) [J/kg]
%   x(2+nL)                 = mv     (vapor mass)           [kg]
%   x(3+nL : 2+nL+nV)      = Tv     (vapor TEMPERATURES)   [K]  <-- KEY CHANGE
%   x(3+nL+nV)              = Ts     (film temperature)     [K]
%   x(4+nL+nV)              = Tw     (wall temperature)     [K]
%   x(5+nL+nV : 12+nL+nV)  = aux    (diagnostic outputs)

nL = P.nL;
nV = P.nV;

% Liquid initial internal energies (uniform)
UL0 = refpropm('U','T',P.TL0,'Q',0,'PARAHYD') * ones(nL,1);

% Vapor initial TEMPERATURES (uniform)
Tv0_vec = P.Tv0 * ones(nV,1);

x0 = [P.mL0;        % 1
      UL0;           % 2 : 1+nL
      P.mv0;         % 2+nL
      Tv0_vec;       % 3+nL : 2+nL+nV   (TEMPERATURES, not energies)
      P.Ts0;         % 3+nL+nV
      P.Tw0;         % 4+nL+nV
      zeros(8,1)];   % auxiliary

% =========================================================================
%  3. Global vent valve state
% =========================================================================
global VentStateGlobal;
VentStateGlobal = P.VentState;

% =========================================================================
%  4. ODE integration with event detection (vent valve toggling)
% =========================================================================
tstart = 0;
tout   = tstart;
xout   = [x0', VentStateGlobal];
teout  = [];
tfinal = P.tFinal;

% ---- Nested ODE right-hand-side ----------------------------------------
function dxdt = odefun(t, x)

    % Throttled waitbar
    persistent lastDay;
    if isempty(lastDay), lastDay = -1; end
    day = floor(t / 86400);
    if day > lastDay
        lastDay = day;
        waitbar(t/P.tFinal, P.waitbar, ...
            sprintf('Day %.1f / %.0f', t/86400, P.tFinal/86400));
    end

    % --- Unpack state ----------------------------------------------------
    mL  = x(1);
    uL  = x(2 : 1+nL);
    mv  = max(x(2+nL), 0.01);          % vapor mass [kg]
    Tv  = x(3+nL : 2+nL+nV);           % vapor TEMPERATURES [K]
    Ts  = max(x(3+nL+nV), 14);         % film temperature [K]
    Tw  = x(4+nL+nV);                  % wall temperature [K]

    % Clamp vapor temperatures to physical range
    for iv = 1:nV
        Tv(iv) = max(Tv(iv), 14);
        Tv(iv) = min(Tv(iv), 100);     % para-H2 shouldn't exceed this
    end

    % --- Liquid properties (polynomials, no REFPROP) ---------------------
    rho_L = polyval_rhoL(uL(nL));
    VL      = mL / rho_L;
    hL      = VL / P.A;
    Vullage = max(P.VTotal - VL, 0.001);
    rhov    = max(mv / Vullage, 0.01);

    TL = zeros(1, nL);
    for iL = 1:nL
        TL(iL) = polyval_TL(uL(iL));
        TL(iL) = max(13.804, min(32.93, TL(iL)));
    end

    % --- Vapor pressure and properties ------------------------------------
    % Check if (T,D) state is two-phase. If so, use saturated vapor
    % properties via (T,Q=1). Transport properties are undefined for
    % two-phase mixtures in REFPROP.
    quality_v = refpropm('q','T',Tv(nV),'D',rhov,'PARAHYD');
    is_twophase = (quality_v >= 0 && quality_v < 1);

    if is_twophase
        % State is inside the dome -- use saturated vapor properties
        pv      = refpropm('P','T',Tv(nV),'Q',1,'PARAHYD') * 1e3;
        Pr_v    = refpropm('^','T',Tv(nV),'Q',1,'PARAHYD');
        kappa_v = refpropm('L','T',Tv(nV),'Q',1,'PARAHYD');
        mu_v    = refpropm('V','T',Tv(nV),'Q',1,'PARAHYD');
        cv_v    = refpropm('O','T',Tv(nV),'Q',1,'PARAHYD');
        cp_v    = refpropm('C','T',Tv(nV),'Q',1,'PARAHYD');
        beta_v  = refpropm('B','T',Tv(nV),'Q',1,'PARAHYD');
    else
        % Single-phase vapor -- use (T,D) directly
        pv      = refpropm('P','T',Tv(nV),'D',rhov,'PARAHYD') * 1e3;
        Pr_v    = refpropm('^','T',Tv(nV),'D',rhov,'PARAHYD');
        kappa_v = refpropm('L','T',Tv(nV),'D',rhov,'PARAHYD');
        mu_v    = refpropm('V','T',Tv(nV),'D',rhov,'PARAHYD');
        cv_v    = refpropm('O','T',Tv(nV),'D',rhov,'PARAHYD');
        cp_v    = refpropm('C','T',Tv(nV),'D',rhov,'PARAHYD');
        beta_v  = refpropm('B','T',Tv(nV),'D',rhov,'PARAHYD');
    end
    pv = max(pv, 1000);

    % --- Film temperature dynamics ---------------------------------------
    Ts0_eq = P.T_c * (pv / P.p_c)^(1/P.lambda);
    dTsdt  = (Ts0_eq - Ts) / P.tminL;

    % Enthalpy of vaporization (polynomial, no REFPROP)
    qh = 1000 * (-0.002445451720487*Ts^6 ...
                 + 0.3629946692976*Ts^5 ...
                 - 22.28028769483*Ts^4 ...
                 + 723.6541112107*Ts^3 ...
                 - 13116.31006512*Ts^2 ...
                 + 125780.2915522*Ts ...
                 - 498095.5392318);

    % --- Vent valve flow -------------------------------------------------
    Jvvalve = VentStateGlobal * gasFlow(P.S_valve, P.gamma_, rhov, pv, P.p_atm);

    % =====================================================================
    %  Wall-to-fluid heat transfer (Nusselt correlations)
    % =====================================================================
    % Vapor side (FIXED Rayleigh: nu^2 in denominator)
    nu_v  = mu_v / rhov;
    H_vap = max(P.H - hL, 0.001);
    Ra_v  = abs(P.g * beta_v * (Tw - Tv(nV)) * H_vap^3 * Pr_v / nu_v^2);
    Psi_v = (1 + (0.492/Pr_v)^(9/16))^(-16/9);
    Nu_v  = 0.68 + 0.503 * (Ra_v * Psi_v)^(1/4);

    % Liquid side transport properties via (T, Q=0) -- always robust
    Pr_L    = refpropm('^','T',TL(nL),'Q',0,'PARAHYD');
    kappa_L = refpropm('L','T',TL(nL),'Q',0,'PARAHYD');
    mu_L    = refpropm('V','T',TL(nL),'Q',0,'PARAHYD');
    cv_L    = refpropm('O','T',TL(nL),'Q',0,'PARAHYD');
    cp_L    = refpropm('C','T',TL(nL),'Q',0,'PARAHYD');
    beta_L  = refpropm('B','T',TL(nL),'Q',0,'PARAHYD');

    % Liquid side Rayleigh (FIXED: nu^2)
    nu_L        = mu_L / rho_L;
    hL_safe     = max(hL, 0.001);
    Ra_L_side   = abs(P.g * beta_L * (Tw - TL(nL)) * hL_safe^3 * Pr_L / nu_L^2);
    Psi_L_side  = (1 + (0.492/Pr_L)^(9/16))^(-16/9);
    Nu_L_side   = 0.68 + 0.503 * (Ra_L_side * Psi_L_side)^(1/4);

    Ra_L_bottom = abs(P.g * beta_L * (Tw - TL(nL)) * (P.R/2)^3 * Pr_L / nu_L^2);
    Nu_L_bottom = 0.27 * Ra_L_bottom^(1/4);

    % Heat transfer coefficients
    hWV        = Nu_v * kappa_v / H_vap;
    hWL_side   = Nu_L_side * kappa_L / hL_safe;
    hWL_bottom = Nu_L_bottom * kappa_L / (P.R/2);

    % Wall-to-fluid heat flows
    QdotWL = (Tw - TL(nL)) * (hWL_bottom * P.A + hWL_side * 2*pi*P.R*hL_safe);
    QdotWV = hWV * (Tw - Tv(nV)) * (P.A + 2*pi*P.R*H_vap);

    % =====================================================================
    %  Film-to-fluid heat transfer
    % =====================================================================
    % Boundary layer grids (liquid)
    lmin_L   = sqrt(kappa_L * P.tminL / cv_L / rho_L);
    l_L(1)   = lmin_L / (1 + exp(pi/2/sqrt(nL)));
    l12_L(1) = lmin_L;
    for i = 2:nL
        l12_L(i) = l12_L(i-1) * exp(pi/sqrt(nL));
        l_L(i)   = sqrt(l12_L(i-1) * l12_L(i));
    end

    % Boundary layer grids (vapor)
    lmin_V   = sqrt(kappa_v * P.tminV / cv_v / rhov);
    l_V(1)   = lmin_V / (1 + exp(pi/2/sqrt(nV)));
    l12_V(1) = lmin_V;
    for i = 2:nV
        l12_V(i) = l12_V(i-1) * exp(pi/sqrt(nV));
        l_V(i)   = sqrt(l12_V(i-1) * l12_V(i));
    end

    % Film heat transfer coefficients
    hVS_cond = kappa_v / l12_V(1);
    hVS_conv = kappa_v * 0.156 * (P.g * beta_v * cp_v * rhov^2 ...
               * max(Ts - Tv(nV), 0) / kappa_v / mu_v)^(1/3);
    hLS_cond = kappa_L / l12_L(1);
    hLS_conv = kappa_L * 0.156 * (P.g * beta_L * cp_L * rho_L^2 ...
               * abs(TL(nL) - Ts) / kappa_L / mu_L)^(1/3);

    % Film-liquid heat flow
    QdotLS_cond = hLS_cond * P.A * (TL(1) - Ts) - l_L(1)*cp_L*rho_L*dTsdt;
    QdotLS_conv = hLS_conv * P.A * (TL(1) - Ts) * (TL(nL) > Ts);
    if QdotLS_conv > 0
        QdotLS = max(QdotLS_conv, QdotLS_cond);
    else
        QdotLS = QdotLS_cond;
    end

    % Film-vapor heat flow
    QdotVS_conv = hVS_conv * P.A * (Tv(1) - Ts) * (Ts > Tv(nV));
    QdotVS_cond = hVS_cond * P.A * (Tv(1) - Ts) - l_V(1)*cv_v*rhov*dTsdt;
    if QdotVS_conv < 0
        QdotVS = min(QdotVS_conv, QdotVS_cond);
    else
        QdotVS = QdotVS_cond;
    end

    % =====================================================================
    %  Condensation
    % =====================================================================
    if qh <= 0
        Jcd = 0;
    else
        Jcd = -(QdotLS + QdotVS) / qh;
    end

    % Guard: limit condensation when vapor mass is low
    if Jcd > 0 && mv < 0.1
        Jcd = Jcd * max(0, (mv - 0.01) / (0.1 - 0.01));
    end

    % =====================================================================
    %  Mass balances
    % =====================================================================
    JL = Jcd;
    Jv = -Jcd - Jvvalve;

    % =====================================================================
    %  pdV work
    % =====================================================================
    pdV = -pv * (JL / rho_L);

    % =====================================================================
    %  Enthalpy terms (all via (T,Q) or (T,D) -- robust)
    % =====================================================================
    if Ts > 32
        hcd = P.c_p * Ts;
    else
        hcd = refpropm('H','T',Ts,'Q',1,'PARAHYD');
    end

    if is_twophase
        hvalve = refpropm('H','T',Tv(nV),'Q',1,'PARAHYD');
    else
        hvalve = refpropm('H','T',Tv(nV),'D',rhov,'PARAHYD');
    end
    vv = Jvvalve / max(P.S_valve * rhov, 1e-6);

    % =====================================================================
    %  Total heat flows to each phase
    % =====================================================================
    QdotV = QdotWV - QdotVS - pdV ...
          - Jvvalve * (hvalve + 0.5*vv^2) ...
          - Jcd * hcd;

    QdotL = QdotWL - QdotLS + pdV ...
          + Jcd * hcd;

    % =====================================================================
    %  LIQUID internal energy derivatives (unchanged from Petitpas)
    % =====================================================================
    duLdt = zeros(nL, 1);
    for i = 1:nL-1
        if i == 1, TLim1 = Ts; else, TLim1 = TL(i-1); end
        rho_Li = refpropm('D','T',TL(i),'Q',0,'PARAHYD');
        duLdt(i) = ((TL(i+1)-TL(i))/l12_L(i+1) - (TL(i)-TLim1)/l12_L(i)) ...
                   * kappa_L / (l_L(i) * rho_Li);
    end
    duLdt(nL) = (QdotL - JL * refpropm('U','T',TL(nL),'Q',0,'PARAHYD')) / mL;

    % =====================================================================
    %  VAPOR TEMPERATURE derivatives (KEY CHANGE: T instead of u)
    % =====================================================================
    % Boundary layers: heat equation dT/dt = alpha * d2T/dx2
    % where alpha = kappa / (rho * cv)
    dTvdt = zeros(nV, 1);
    for i = 1:nV-1
        if i == 1, Tvim1 = Ts; else, Tvim1 = Tv(i-1); end
        % Check if this layer is two-phase at current pressure
        q_layer = refpropm('q','T',Tv(i),'P',pv/1000,'PARAHYD');
        if q_layer >= 0 && q_layer < 1
            rhovi = refpropm('D','T',Tv(i),'Q',1,'PARAHYD');
            cv_vi = refpropm('O','T',Tv(i),'Q',1,'PARAHYD');
        else
            rhovi = refpropm('D','T',Tv(i),'P',pv/1000,'PARAHYD');
            cv_vi = refpropm('O','T',Tv(i),'P',pv/1000,'PARAHYD');
        end
        dTvdt(i) = ((Tv(i+1)-Tv(i))/l12_V(i+1) - (Tv(i)-Tvim1)/l12_V(i)) ...
                   * kappa_v / (l_V(i) * rhovi * cv_vi);
    end

    % Bulk vapor energy balance, converted to dT/dt:
    %   mv * cv * dTv/dt = QdotV - Jv * u(Tv, rhov)
    if is_twophase
        uv_bulk = refpropm('U','T',Tv(nV),'Q',1,'PARAHYD');
    else
        uv_bulk = refpropm('U','T',Tv(nV),'D',rhov,'PARAHYD');
    end
    dTvdt(nV) = (QdotV - Jv * uv_bulk) / (mv * cv_v);

    % =====================================================================
    %  Wall temperature
    % =====================================================================
    cw    = 2.516173240451E-11*Tw^6 - 2.695483209737E-08*Tw^5 ...
          + 0.00001122596286143*Tw^4 - 0.002261465800734*Tw^3 ...
          + 0.214810433559*Tw^2 - 5.41715155529*Tw + 51.75489930095;
    dcwdT = 6*2.516173240451E-11*Tw^5 - 5*2.695483209737E-08*Tw^4 ...
          + 4*0.00001122596286143*Tw^3 - 3*0.002261465800734*Tw^2 ...
          + 2*0.214810433559*Tw - 5.41715155529;

    if strcmp(P.heatLeakModel, 'petitpas')
        QdotEW = -7.462776654302E-02*VL^2 + 4.445867251697E+00*VL ...
               + 3.108170556297E+01;
    else
        QdotEW = P.QdotEW_const;
    end

    dTwdt = (QdotEW - QdotWL - QdotWV) / (P.mw * (cw + Tw * dcwdT));

    % =====================================================================
    %  Assemble derivative vector
    % =====================================================================
    dxdt = zeros(size(x));

    dxdt(1)                = JL;
    dxdt(2 : 1+nL)        = duLdt;
    dxdt(2+nL)             = Jv;
    dxdt(3+nL : 2+nL+nV)  = dTvdt;     % TEMPERATURES, not energies
    dxdt(3+nL+nV)          = dTsdt;
    dxdt(4+nL+nV)          = dTwdt;

    % Auxiliary (first-order filtered for smooth output)
    base = 4 + nL + nV;
    dxdt(base+1) = (Jcd     - x(base+1));
    dxdt(base+2) = (QdotWV  - x(base+2));
    dxdt(base+3) = (QdotWL  - x(base+3));
    dxdt(base+4) = (QdotVS  - x(base+4));
    dxdt(base+5) = (QdotLS  - x(base+5));
    dxdt(base+6) = (pdV     - x(base+6));
    dxdt(base+7) = (QdotEW  - x(base+7));
    dxdt(base+8) = (Jvvalve - x(base+8));

end % end of odefun

% =========================================================================
%  5. Integration loop (restarts at each vent event)
% =========================================================================
while tout(end) < P.tFinal

    evtFcn  = @(t,x) ventEvents(x, P, VentStateGlobal, nL, nV);
    rhsFcn  = @(t,x) odefun(t, x);

    options = odeset('RelTol', P.relTol, 'Events', evtFcn, ...
                     'Refine', 4, 'MaxStep', 3600);
    if tstart > 0
        options = odeset(options, 'InitialStep', 10);
    end

    [t, x, te, xe, ie] = ode15s(rhsFcn, [tstart, tfinal], x0, options);

    % Accumulate output
    nt = length(t);
    ventcol = VentStateGlobal * ones(nt, 1);
    tout = [tout; t(2:nt)];
    xout = [xout; x(2:nt,:), ventcol(2:nt)];

    % Prepare for next segment
    x0     = x(end, :)';
    tstart = t(nt);

    % Check which event fired
    if ~isempty(ie) && any(ie == 3)
        % Liquid depleted -- stop simulation
        fprintf('Liquid depleted at day %.1f. Stopping.\n', tstart/86400);
        break;
    end

    % Toggle vent valve (events 1 or 2)
    VentStateGlobal = abs(VentStateGlobal - 1);
end

close(P.waitbar);
fprintf('ODE integration complete (%.1f days simulated).\n', tout(end)/86400);

% =========================================================================
%  6. Package output data
% =========================================================================
data.name = name;
data.t    = tout;

data.mL  = xout(:, 1);
data.uL  = xout(:, 2 : 1+nL);
data.mv  = xout(:, 2+nL);
data.Tv  = xout(:, 3+nL : 2+nL+nV);     % TEMPERATURES directly
data.Ts  = xout(:, 3+nL+nV);
data.Tw  = xout(:, 4+nL+nV);

base = 4 + nL + nV;
data.Jcd     = xout(:, base+1);
data.QdotWV  = xout(:, base+2);
data.QdotWL  = xout(:, base+3);
data.QdotVS  = xout(:, base+4);
data.QdotLS  = xout(:, base+5);
data.pdV     = xout(:, base+6);
data.QdotEW  = xout(:, base+7);
data.Jvvalve = xout(:, base+8);

data.VentState = xout(:, base+9);

end % end of LH2StaticSimulate


% =========================================================================
%  Local functions
% =========================================================================

function rho = polyval_rhoL(uL_val)
% Liquid density from internal energy (REFPROP v9.1 polynomial fit)
    u = uL_val / 1000;
    rho = -5.12074746E-07*u^3 - 1.56628367E-05*u^2 ...
          - 1.18436797E-01*u + 7.06218354E+01;
    rho = max(rho, 50);
end

function T = polyval_TL(uL_val)
% Liquid temperature from internal energy (REFPROP v9.1 polynomial fit)
    u = uL_val / 1000;
    T = 1.44867559E-07*u^3 - 2.53438808E-04*u^2 ...
        + 1.05449468E-01*u + 2.03423757E+01;
end

function mdot = gasFlow(CA, gamma, rho, P1, P2)
% Choked / non-choked compressible flow through an orifice.
    if P1 < P2
        mdot = -gasFlow(CA, gamma, rho, P2, P1);
    else
        thresh = ((gamma+1)/2)^(gamma/(gamma-1));
        if P1/P2 >= thresh
            mdot = CA * sqrt(gamma*rho*P1 * (2/(gamma+1))^((gamma+1)/(gamma-1)));
        else
            mdot = CA * sqrt(2*rho*P1*(gamma/(gamma-1)) ...
                   * ((P2/P1)^(2/gamma) - (P2/P1)^((gamma+1)/gamma)));
        end
    end
end

function [value, isterminal, direction] = ventEvents(x, P, ventstate, nL, nV)
% Event function for vent valve hysteresis + liquid depletion.
    rho_L_ev   = polyval_rhoL(x(1+nL));
    VL_ev      = x(1) / rho_L_ev;
    Vullage_ev = max(P.VTotal - VL_ev, 0.001);
    mv_ev      = max(x(2+nL), 0.01);
    rhov_ev    = mv_ev / Vullage_ev;
    Tv_ev      = max(x(2+nL+nV), 14);

    % Vapor pressure -- handle two-phase gracefully
    q_ev = refpropm('q','T',Tv_ev,'D',rhov_ev,'PARAHYD');
    if q_ev >= 0 && q_ev < 1
        pv_ev = refpropm('P','T',Tv_ev,'Q',1,'PARAHYD') * 1e3;
    else
        pv_ev = refpropm('P','T',Tv_ev,'D',rhov_ev,'PARAHYD') * 1e3;
    end

    value = [pv_ev - P.p_low;     % 1: vent close threshold
             pv_ev - P.p_high;    % 2: vent open threshold
             x(1) - 5.0];         % 3: liquid depletion (mL < 5 kg)

    if ventstate > 0
        isterminal = [1; 1; 1];
    else
        isterminal = [0; 1; 1];
    end
    direction = [-1; +1; -1];     % liquid depletion: decreasing
end
