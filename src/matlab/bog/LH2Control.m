function U = LH2Control(fill_frac,p1,p2,ET_fill_complete,ST_vent_complete)
% U = LH2Control(fill_frac,p1,p2,ET_fill_complete,ST_vent_complete)
%   Determines control inputs for given ullage pressures and ET fill fraction.
%
%   (ST) or 1 = Storage tank at liquefaction terminal  (FEEDING vessel)
%   (ET) or 2 = LH2 trailer                           (RECEIVING vessel)
%
%   fill_frac = VL2/VTotal2, the volume fill fraction of the trailer (0..1)
%
% E = Transfer Line Valve
% V = Vaporizer (on ST)

% obtain model parameters structure
P = evalin('base','LH2Model');

% Fill regime thresholds (volume-fraction based, for horizontal ET)
%   0.00 - 0.15: slow fill   (half-open valve, stabilise flow)
%   0.15 - 0.50: fast fill   (full-open valve, full vaporizer)
%   0.50 - 0.80: reduced fast (full-open valve, reduced vaporizer to limit overpressure)
%   0.80 - TopET: topping     (shut off transfer when target reached, vent ST)

if fill_frac < 0.15
    % slow fill
    U.lambdaE = 0.5;
    U.lambdaV = getVaporizerValveState(P,p1,P.p_ST_slow);
    U.STVentState = getSTVentState(P,p1,P.p_ST_slow);
elseif fill_frac < 0.50
    % fast fill
    U.lambdaE = 1;
    U.lambdaV = getVaporizerValveState(P,p1,P.p_ST_fast);
    U.STVentState = getSTVentState(P,p1,P.p_ST_fast);
elseif fill_frac < 0.80
    % reduced fast fill
    U.lambdaE = 1;
    U.lambdaV = getVaporizerValveState(P,p1,P.p_ST_slow);
    U.STVentState = ET_fill_complete * (p1 > P.p_ST_final);
else
    % topping / post-fill
    U.lambdaE = (1-ET_fill_complete);
    U.lambdaV = (1-ET_fill_complete) * getVaporizerValveState(P,p1,P.p_ST_slow);
    U.STVentState = ET_fill_complete*(1-ST_vent_complete);
end


function state = getSTVentState(P,p1,threshold)
% Determine ST vent valve state (hysteresis band +/- 5%)
if p1 < threshold+0.05*threshold
    % turn off valve
    state = 0;
elseif p1 >threshold-0.05*threshold
    % turn on valve
    state = 1;
else
    % stays at same value
    state = P.STVentState;
end

% NOTE: getETVentState is not called -- ET vent state is controlled by the
% ODE event detection in LH2Simulate.m (VentEvents function). Retained
% here for reference / potential future use.
function state = getETVentState(P,p2)
% Determine ET vent valve state (hysteresis between p_ET_low and p_ET_high)
if p2 < P.p_ET_low
    state = 0;
elseif p2 > P.p_ET_high
    state = 1;
else
    state = P.ETVentState;
end


function state = getVaporizerValveState(P,p1,pSet)
% Goal: provide enough flow to maintain p1 at pSet.
% Opens proportionally when below target, closes above.
if p1<pSet-0.02*pSet
    state = max(0,10*(pSet-p1)/pSet);
    state = min(1,state);
elseif p1>pSet+0.02*pSet
    state = 0;
else
    state = P.VapValveState;
end
