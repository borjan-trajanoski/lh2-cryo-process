% inputs_StaticBoiloff.m
% Parameters for single-vessel LH2 static boil-off simulation.
%
% DEFAULT: 3,300-gallon Dewar at LLNL, matching Petitpas (2018) Figs 9-11.

clear LH2Model;

%% Unit conversions
psiToPa = 6894.75729;
galTom3 = 0.00378541;

%% Physical constants
LH2Model.p_atm = 1.01325e5;   % [Pa]
LH2Model.g     = 9.8;         % [m/s^2]

%% Tank geometry (vertical cylinder)
LH2Model.VTotal = 3300 * galTom3;              % [m^3] total volume
LH2Model.R      = 1.0;                         % [m] inner radius
LH2Model.A      = pi * LH2Model.R^2;           % [m^2] cross-section
LH2Model.H      = LH2Model.VTotal / LH2Model.A; % [m] height

%% Critical properties (para-hydrogen)
LH2Model.T_c    = 32.938;                      % [K]
LH2Model.p_c    = 186.49 * psiToPa;            % [Pa]
LH2Model.lambda = 5;                           % exponent for film T

%% Reference vapor properties at ~1 bar (for gasFlow choked-flow calc)
LH2Model.gamma_ = 5/3;
LH2Model.c_p    = 6490 + 4124;                 % [J/kg/K] cp = cv + R
LH2Model.Gamma_ = ((LH2Model.gamma_+1)/2)^((LH2Model.gamma_+1)/2/(LH2Model.gamma_-1));

%% Boundary layer grid parameters
LH2Model.nL    = 3;       % liquid grid size
LH2Model.tminL = 0.1;     % [s] liquid grid time constant
LH2Model.nV    = 4;       % vapor grid size
LH2Model.tminV = 0.1;     % [s] vapor grid time constant

%% Initial conditions
LH2Model.p0  = 20 * psiToPa;    % [Pa] initial vapor pressure
LH2Model.TL0 = 21;              % [K] initial liquid temperature
LH2Model.Tw0 = 21;              % [K] initial wall temperature

% Saturation temperature at initial pressure (polynomial from REFPROP v9.1)
% Add 0.5 K superheat so initial state is clearly single-phase vapor.
p0_psi = LH2Model.p0 / psiToPa;
Tsat0 = -1.603941638811E-11*p0_psi^6 ...
    + 7.830478134841E-09*p0_psi^5 ...
    - 1.549372675881E-06*p0_psi^4 ...
    + 1.614567978153E-04*p0_psi^3 ...
    - 9.861776990784E-03*p0_psi^2 ...
    + 4.314905904166E-01*p0_psi ...
    + 1.559843335080E+01;
LH2Model.Tv0 = Tsat0 + 0.5;   % [K] slight superheat

% Film temperature (Osipov 2008 model)
LH2Model.Ts0 = LH2Model.T_c * (LH2Model.p0 / LH2Model.p_c)^(1/LH2Model.lambda);

% Initial fill: ~70% (2300 of 3300 gallons)
VL0 = 2300 * galTom3;                                                     % [m^3]
LH2Model.rhoL0 = refpropm('D','T',LH2Model.TL0,'Q',0,'PARAHYD');         % [kg/m^3]
LH2Model.mL0   = LH2Model.rhoL0 * VL0;                                   % [kg]
LH2Model.rhov0 = refpropm('D','T',LH2Model.Tv0,'P',LH2Model.p0/1000,'PARAHYD');
LH2Model.mv0   = LH2Model.rhov0 * (LH2Model.VTotal - VL0);              % [kg]

%% Wall properties (SA240 T304 stainless steel inner vessel)
LH2Model.mw = 2529;   % [kg] mass of inner vessel

%% Heat leak model
% 'petitpas': fill-level polynomial from Fig 8 (Summer 2015, 3300-gal)
% 'constant': fixed value in QdotEW_const
LH2Model.heatLeakModel = 'petitpas';
LH2Model.QdotEW_const  = 50;    % [W] only used if model = 'constant'

%% Vent valve (PRD)
LH2Model.S_valve   = pi * (0.4/100)^2;   % [m^2] orifice area (4mm radius)
LH2Model.VentState = 0;                   % initial state (closed)
LH2Model.p_low     = 43 * psiToPa;        % [Pa] vent close pressure
LH2Model.p_high    = 45 * psiToPa;        % [Pa] vent open pressure

%% Solver options
LH2Model.tFinal = 50 * 24 * 3600;    % [s] 50 days (Petitpas validation)
LH2Model.relTol = 1e-4;
