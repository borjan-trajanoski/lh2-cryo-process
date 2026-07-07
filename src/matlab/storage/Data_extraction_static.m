% Data_extraction_static.m
% Post-processing for LH2 static boil-off simulation.
%
% Since the vapor state is now TEMPERATURE (not internal energy),
% post-processing is much simpler -- no REFPROP loops needed for Tv.

fprintf('Extracting derived quantities...\n');

P  = LH2Model;
nL = P.nL;
nV = P.nV;
N  = length(results.t);

%% Liquid density, volume, fill level
results.rho_L = zeros(N,1);
for z = 1:N
    u = results.uL(z,nL) / 1000;
    results.rho_L(z) = -5.12074746E-07*u^3 - 1.56628367E-05*u^2 ...
                       - 1.18436797E-01*u + 7.06218354E+01;
end

results.VL      = results.mL ./ results.rho_L;
results.hL      = results.VL ./ P.A;
results.Vullage = P.VTotal - results.VL;
results.pctFill = results.hL ./ P.H;

%% Vapor density
results.rhov = results.mv ./ results.Vullage;

%% Liquid temperatures (polynomial)
results.TL = zeros(N, nL);
for z = 1:N
    for ii = 1:nL
        u = results.uL(z,ii) / 1000;
        results.TL(z,ii) = 1.44867559E-07*u^3 - 2.53438808E-04*u^2 ...
                          + 1.05449468E-01*u + 2.03423757E+01;
        results.TL(z,ii) = max(13.804, min(32.93, results.TL(z,ii)));
    end
end

%% Vapor pressure via (T,D) flash -- robust
fprintf('Computing vapor pressures (%d points)...\n', N);
results.pv = zeros(N, 1);
for z = 1:N
    Tv_z   = max(results.Tv(z,nV), 14);
    rhov_z = max(results.rhov(z), 0.01);
    q_z = refpropm('q','T',Tv_z,'D',rhov_z,'PARAHYD');
    if q_z >= 0 && q_z < 1
        results.pv(z) = refpropm('P','T',Tv_z,'Q',1,'PARAHYD') * 1e3;
    else
        results.pv(z) = refpropm('P','T',Tv_z,'D',rhov_z,'PARAHYD') * 1e3;
    end

    if mod(z, round(N/10)) == 0
        fprintf('  %d%%\n', round(100*z/N));
    end
end

%% Vent flow (recomputed)
results.Jvvalve_recomp = zeros(N, 1);
for z = 1:N
    results.Jvvalve_recomp(z) = results.VentState(z) ...
        * gasFlow_local(P.S_valve, P.gamma_, results.rhov(z), results.pv(z), P.p_atm);
end

%% BOG: cumulative vented mass and daily rates
results.BOG_cumulative = zeros(N, 1);
results.BOG_rate       = zeros(N, 1);
for z = 2:N
    dt = results.t(z) - results.t(z-1);
    results.BOG_rate(z)       = results.Jvvalve_recomp(z);
    results.BOG_cumulative(z) = results.BOG_cumulative(z-1) + dt * results.BOG_rate(z);
end

%% Total mass and mass loss
results.mTotal   = results.mL + results.mv;
results.massLoss = results.mTotal(1) - results.mTotal;

%% Daily BOG averages
tDays = results.t / 86400;
if tDays(end) > 1
    nDays = floor(tDays(end));
    results.dailyBOG_kg  = zeros(nDays, 1);
    results.dailyBOG_pct = zeros(nDays, 1);
    for d = 1:nDays
        idx = (tDays >= (d-1)) & (tDays < d);
        if any(idx)
            results.dailyBOG_kg(d)  = trapz(results.t(idx), results.BOG_rate(idx));
            results.dailyBOG_pct(d) = 100 * results.dailyBOG_kg(d) / results.mTotal(1);
        end
    end
end

fprintf('Data extraction complete.\n');

%% Local function
function mdot = gasFlow_local(CA, gamma, rho, P1, P2)
    if P1 < P2
        mdot = -gasFlow_local(CA, gamma, rho, P2, P1);
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
