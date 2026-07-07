% Defines LH2 model parameters for LOADING AT LIQUID TERMINAL
% Adapted from Petitpas (2018) inputs_TrailerToDewar.m
%
% Configuration:
%   (ST) or 1 = 5 x 350 m3 manifolded horizontal storage (FEEDING vessel)
%   (ET) or 2 = 60 m3 horizontal LH2 trailer (RECEIVING vessel)
%
% Scenario: ONE representative truck loading event at an 86 TPD plant.
%   The 5 x 350 m3 Air Liquide horizontal tanks are manifolded (common
%   liquid + vapor headers) and operate as a single 1750 m3 buffer. The
%   liquefier feeds continuously (~1 kg/s), trucks draw periodically
%   (~22 trucks/day, ~4000 kg each, ~65 min apart).
%
%   STEADY-STATE ASSUMPTION: The manifold is ALREADY at its operating
%   pressure (2.5 bar). The vaporizer only maintains pressure as liquid
%   is withdrawn -- there is NO pressurization startup transient and NO
%   depressurization between trucks. This represents the steady-state
%   condition during normal plant operation. Per-truck depletion is only
%   4000 kg out of ~111,000 kg (3.6%), so source conditions barely budge.
%
%   The result of this simulation is multiplied by 22 for daily losses.
%
% References:
%   - PRESLHY D2.3 (Jallais, 2018): Air Liquide Waziers tank dimensions
%   - Cryolor LH2 brochure: European trailer specifications
%   - Petitpas (2018): Transfer model and heat leak estimates
%   - Linde/Decker (2019): ~1 hour loading time for LH2 trailers
%   - BoilFAST steady-state results: 648 W per 350 m3 tank

clear LH2Model;

%% ===== Unit Conversions =====
psiToPa = 6894.75729;        % psi to Pascals
galTom3 = 0.00378541;        % gallons to cubic meters
inTom   = 0.0254;            % inches to meters
barToPa = 1e5;               % bar to Pascals

%% ===== Physical Constants =====
LH2Model.p_atm = 1.01325e5;  % [Pa] atmospheric pressure
LH2Model.g     = 9.8;        % [m/s^2] gravitational acceleration

%% ===== Tank Geometry =====
% --- ST: 5 x 350 m3 Manifolded Industrial Horizontal Storage ---
% Based on Air Liquide Waziers design (PRESLHY D2.3)
% Inner diameter = 4.0 m, evacuated perlite insulation 500 mm
% Modeled as single equivalent horizontal cylinder (R = 2.0 m, L = 139.3 m)
% NOTE: Interface area is wrong for a single long cylinder vs 5 short ones,
% but source-side dynamics are negligible (3.6% depletion per truck).
LH2Model.VTotal1 = 1750;                                  % [m^3] manifolded volume (5 x 350)
LH2Model.R1      = 2.0;                                   % [m] inner radius (physical)
LH2Model.A1      = pi * LH2Model.R1^2;                    % [m^2] cross-section area
LH2Model.Lcyl    = LH2Model.VTotal1 / LH2Model.A1;       % [m] equivalent length (~139.3 m)

% --- ET: 60 m3 European LH2 Trailer ---
% Based on Cryolor/Linde European trailer specifications
% Payload ~3,800 kg at 90% fill, MAWP ~10 bar
LH2Model.VTotal2 = 60;                                    % [m^3] total volume
LH2Model.R2      = 1.15;                                  % [m] inner radius (ID ~2.3 m)
LH2Model.A2      = pi * LH2Model.R2^2;                    % [m^2] cross-section area
LH2Model.Lcyl2   = LH2Model.VTotal2 / LH2Model.A2;       % [m] cylinder length (~14.4 m)
LH2Model.H       = 2 * LH2Model.R2;                       % [m] diameter (used as reference height)

% Flag to indicate ET is a horizontal cylinder (used in LH2Simulate)
LH2Model.ET_horizontal = 1;

%% ===== LH2 Fluid Properties =====
% Para-hydrogen properties from REFPROP v9.1
LH2Model.T_c     = 32.938;                 % [K] critical temperature
LH2Model.p_c     = 186.49 * psiToPa;       % [Pa] critical pressure
LH2Model.lambda  = 5;                      % [-] exponent for film temperature
LH2Model.rho_L   = 70.9;                   % [kg/m^3] liquid density @ 1 bar
LH2Model.c_L     = 9702.5;                 % [J/kg/K] liquid specific heat @ 1 bar
LH2Model.kappa_L = 0.10061;                % [W/mK] liquid thermal conductivity @ 1 bar
LH2Model.mu_L    = 13.54e-6;               % [Pa*s] liquid dynamic viscosity @ 1 bar

%% ===== GH2 (Vapor) Properties =====
LH2Model.R_v     = 4124;                   % [J/kg/K] specific gas constant
LH2Model.c_v     = 6490;                   % [J/kg/K] Cv (rotational DOF frozen)
LH2Model.c_p     = LH2Model.c_v + LH2Model.R_v;
LH2Model.gamma_  = 5/3;                    % [-] ratio of specific heats
LH2Model.Gamma_  = ((LH2Model.gamma_+1)/2)^((LH2Model.gamma_+1)/2/(LH2Model.gamma_-1));
LH2Model.mu_v    = 0.98e-6;                % [Pa*s] vapor dynamic viscosity @ 1 bar
LH2Model.kappa_v = 0.0166;                 % [W/mK] vapor thermal conductivity @ 1 bar

%% ===== Grid Parameters =====
% Boundary layer grid for condensation/evaporation model
LH2Model.nL1    = 3;      LH2Model.tminL1 = 0.1;   % ST liquid grid
LH2Model.nL2    = 3;      LH2Model.tminL2 = 0.1;   % ET liquid grid
LH2Model.nV1    = 3;      LH2Model.tminV1 = 0.1;   % ST vapor grid
LH2Model.nV2    = 4;      LH2Model.tminV2 = 0.1;   % ET vapor grid

%% ===== Initial Conditions: ST (Manifolded Storage, STEADY-STATE) =====
% The manifold is ALREADY at its continuous operating pressure of 2.5 bar.
% No pressurization startup. The vaporizer just maintains pressure as
% liquid is withdrawn (~4000 kg out of ~111,000 kg = 3.6% depletion).
LH2Model.p10        = 2.5 * barToPa;       % [Pa] already at operating pressure
LH2Model.TL10       = 21.0;                % [K] initial liquid temperature
LH2Model.totalmass10 = 0.90 * 1750 * 70.9; % [kg] 90% fill at rho_L ~70.9 kg/m3 (~111,667 kg)

% Compute saturation temperature of vapor from pressure (REFPROP polynomial)
LH2Model.Tv10 = 0.1 + (-1.603941638811E-11*(LH2Model.p10/psiToPa)^6 ...
    + 7.830478134841E-09*(LH2Model.p10/psiToPa)^5 ...
    - 1.549372675881E-06*(LH2Model.p10/psiToPa)^4 ...
    + 1.614567978153E-04*(LH2Model.p10/psiToPa)^3 ...
    - 9.861776990784E-03*(LH2Model.p10/psiToPa)^2 ...
    + 4.314905904166E-01*(LH2Model.p10/psiToPa)^1 ...
    + 1.559843335080E+01);

LH2Model.Ts10    = LH2Model.T_c * (LH2Model.p10/LH2Model.p_c)^(1/LH2Model.lambda);
LH2Model.rhov10  = refpropm('D','T',LH2Model.Tv10,'P',LH2Model.p10/1000,'PARAHYD');
LH2Model.rhoL10  = refpropm('D','T',LH2Model.TL10,'Q',0,'PARAHYD');
LH2Model.Vullage10 = (LH2Model.totalmass10 - LH2Model.rhoL10*LH2Model.VTotal1) / ...
                      (LH2Model.rhov10 - LH2Model.rhoL10);
LH2Model.mL10    = LH2Model.rhoL10 * (LH2Model.VTotal1 - LH2Model.Vullage10);
LH2Model.mv10    = LH2Model.totalmass10 - LH2Model.mL10;

%% ===== Initial Conditions: ET (Trailer) =====
% Freshly purged trailer: low pressure, small residual liquid heel.
% Returned from delivery, depressurized for road transport.
LH2Model.p20        = 1.4 * barToPa;       % [Pa] initial pressure (1.4 bar, ~20.3 psia)
LH2Model.TL20       = 20.4;                % [K] initial liquid temperature
LH2Model.Tw20       = 21.0;                % [K] initial wall temperature (cold trailer)
LH2Model.pct_hL20   = 0.10;                % [-] initial fill fraction (10% heel)

LH2Model.Tv20 = 0.1 + (-1.603941638811E-11*(LH2Model.p20/psiToPa)^6 ...
    + 7.830478134841E-09*(LH2Model.p20/psiToPa)^5 ...
    - 1.549372675881E-06*(LH2Model.p20/psiToPa)^4 ...
    + 1.614567978153E-04*(LH2Model.p20/psiToPa)^3 ...
    - 9.861776990784E-03*(LH2Model.p20/psiToPa)^2 ...
    + 4.314905904166E-01*(LH2Model.p20/psiToPa)^1 ...
    + 1.559843335080E+01);

LH2Model.Ts20    = LH2Model.T_c * (LH2Model.p20/LH2Model.p_c)^(1/LH2Model.lambda);

% For horizontal trailer, initial liquid volume from fill fraction
LH2Model.VL20    = LH2Model.pct_hL20 * LH2Model.VTotal2;              % [m^3]
LH2Model.rhoL20  = refpropm('D','T',LH2Model.TL20,'Q',0,'PARAHYD');
LH2Model.mL20    = LH2Model.rhoL20 * LH2Model.VL20;                   % [kg]
LH2Model.rhov20  = refpropm('D','T',LH2Model.Tv20,'P',LH2Model.p20/1000,'PARAHYD');
LH2Model.mv20    = LH2Model.rhov20 * (LH2Model.VTotal2 - LH2Model.VL20);

%% ===== Initial Flows =====
LH2Model.Jboil0 = 0;       % [kg/s] initial vaporizer boiling flow
LH2Model.Jtr0   = 0;       % [kg/s] initial transfer line flow

%% ===== Heat Transfer: ST (Manifolded Storage) =====
% Total heat leak: 5 x 648 W = 3240 W (from BoilFAST, 350 m3 tanks with
% U = 0.005 W/(m2K), evacuated perlite). Split 83/17 liquid/vapor by
% wetted area at 90% fill (following Petitpas 200/40 W ratio).
LH2Model.QdotEL1 = 2690;   % [W] environment to liquid (83% of 3240 W)
LH2Model.QdotEV1 = 550;    % [W] environment to vapor  (17% of 3240 W)

%% ===== Heat Transfer: ET (Trailer) =====
% Heat leak ~270 W (Petitpas 2018 estimate for near-full trailer)
LH2Model.QdotEW2 = 270;    % [W] total environment to trailer wall (constant)
LH2Model.mw2     = 4000;   % [kg] estimated mass of trailer inner vessel
                            % Scaled from Petitpas Dewar (2529 kg for 12.5 m3)
                            % to 60 m3 trailer, with lighter construction

%% ===== Vaporizer Parameters (on ST) =====
% In steady-state, the vaporizer only MAINTAINS pressure at 2.5 bar as
% liquid is withdrawn. With a 1750 m3 manifold losing only 57 m3 of
% liquid per truck, the ullage expansion is small and the pressure drop
% is modest. Coefficient scaled up from original (3e-4 for 350 m3) to
% account for the larger system, though the vaporizer does less work
% per event now because the manifold acts as a pressure buffer.
LH2Model.mVap0      = 0;       % [kg] initial liquid mass in vaporizer
LH2Model.Tboil      = 22;      % [K] temperature of vapor bubbles
LH2Model.tau_vap    = 2.0;     % [s] vaporizer time constant
LH2Model.c_vap      = 5e-4;    % [-] vaporizer valve flow coefficient (scaled for 1750 m3)
LH2Model.VapValveState = 0;    % initial valve state

%% ===== Transmission Line Parameters =====
% Industrial cryogenic transfer line between storage tank and trailer
% loading bay. 1.25-inch vacuum-jacketed piping is representative of
% fixed loading arm installations at European LH2 terminals.
%
% CHILL-DOWN ASSUMPTION: At an 86 TPD plant loading ~22 trucks/day
% (one truck every ~65 min), the transfer line stays cold between
% consecutive loadings. Pipeline chill-down losses are therefore
% negligible and not modelled. This assumption would NOT hold for
% infrequent deliveries where the line warms back to ambient.
LH2Model.DPipe  = 1.25 * inTom;   % [m] pipe diameter (1.25 inch, vacuum-jacketed)
LH2Model.LPipe  = 15;             % [m] transmission line length
LH2Model.drPipe = 1e-6;           % [-] pipe roughness
LH2Model.f      = 1.3 / log(LH2Model.DPipe/2/LH2Model.drPipe)^2;
LH2Model.tau_tr = 15;             % [s] transmission line delay constant
LH2Model.dE     = 1.25 * inTom;   % [m] transfer line valve diameter (matched to pipe)
LH2Model.kE     = 8;              % [-] transfer line valve loss coefficient

%% ===== Vent Valves =====
LH2Model.S_valve1    = 0.005;  % [m^2] ST vent valve orifice area (scaled for 1750 m3 manifold)
LH2Model.STVentState = 0;      % initial ST vent state (closed)
LH2Model.S_valve2    = 0.0005; % [m^2] ET vent valve orifice area
LH2Model.ETVentState = 0;      % initial ET vent state (closed)

%% ===== Pressure Settings =====
% STEADY-STATE: The manifold stays at 2.5 bar continuously.
% p_ST_final = p_ST_fast: NO depressurization between trucks.
% The next truck connects and loading continues immediately.
%
% Delivery pressure of 2.5 bar gives margin above the ET vent
% setpoint (2.5 bar). With the 1750 m3 manifold acting as a pressure
% buffer, the driving pressure stays stable throughout filling.
LH2Model.p_ST_slow  = 2.5 * barToPa;  % [Pa] vaporizer target for slow fill
LH2Model.p_ST_fast  = 2.5 * barToPa;  % [Pa] vaporizer target for fast fill
LH2Model.p_ST_final = 2.5 * barToPa;  % [Pa] stays at operating pressure (NO depressurization)

% Vent valve thresholds for ET (trailer) -- PRD hysteresis.
% As liquid enters the trailer, ullage compression raises the pressure.
% Venting captures the BOG and returns it to the plant feed compressor.
% Wide hysteresis band (0.5 bar) reduces cycling frequency, preventing
% solver stalling from rapid open/close oscillation.
LH2Model.p_ET_low   = 1.5 * barToPa;  % [Pa] ET vent close pressure (1.5 bar)
LH2Model.p_ET_high  = 1.8 * barToPa;  % [Pa] ET vent open pressure (1.8 bar)

%% ===== Fill Control =====
LH2Model.TopET = 0.90;                 % [-] maximum fill fraction for ET (90%)
LH2Model.ratio_top_bottom = 0.0;       % [-] fraction of liquid going to top fill (0 = all bottom)

%% ===== Solver Options =====
% No pressurization transient needed (ST already at 2.5 bar). Allow
% 90 minutes: ~50-60 min fill + ~10 min post-fill ET settling.
% Note: ST does NOT depressurize, so no post-fill ST vent phase.
LH2Model.tFinal = 90 * 60;             % [s] simulation time (90 minutes)
LH2Model.relTol = 5e-4;                % relative tolerance for ODE solver

%% ===== Display Summary =====
fprintf('\n=== Manifolded Storage to Trailer Loading (Steady-State) ===\n');
fprintf('ST (Manifold): %.0f m3 (5 x 350), R=%.2f m, L_eq=%.1f m\n', ...
    LH2Model.VTotal1, LH2Model.R1, LH2Model.Lcyl);
fprintf('  Initial: %.0f kg total, %.1f bar (steady-state), %.1f K\n', ...
    LH2Model.totalmass10, LH2Model.p10/barToPa, LH2Model.TL10);
fprintf('  Fill: %.0f%%, Heat leak: %.0f W (5 tanks)\n', ...
    (1 - LH2Model.Vullage10/LH2Model.VTotal1)*100, LH2Model.QdotEL1+LH2Model.QdotEV1);
fprintf('  p_ST_final = p_ST_fast = %.1f bar (no depressurization)\n', ...
    LH2Model.p_ST_final/barToPa);
fprintf('ET (Trailer):  %.0f m3, R=%.2f m, L=%.1f m\n', ...
    LH2Model.VTotal2, LH2Model.R2, LH2Model.Lcyl2);
fprintf('  Initial: %.1f kg liquid + %.1f kg vapor, %.1f bar, %.1f K\n', ...
    LH2Model.mL20, LH2Model.mv20, LH2Model.p20/barToPa, LH2Model.TL20);
fprintf('  Target fill: %.0f%% (%.0f kg LH2)\n', ...
    LH2Model.TopET*100, LH2Model.TopET*LH2Model.VTotal2*LH2Model.rho_L*0.99);
fprintf('Delivery pressure: %.1f bar\n', LH2Model.p_ST_slow/barToPa);
fprintf('Transfer line: %.2f inch ID, %.0f m length\n', ...
    LH2Model.DPipe/inTom, LH2Model.LPipe);
fprintf('Trailer vent: %.1f / %.1f bar (low/high)\n', ...
    LH2Model.p_ET_low/barToPa, LH2Model.p_ET_high/barToPa);
fprintf('=============================================================\n\n');
