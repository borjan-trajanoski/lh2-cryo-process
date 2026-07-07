function plotStaticBoiloff(data)
% plotStaticBoiloff  Diagnostic and results plots for static boil-off.

P = evalin('base','LH2Model');
psiToPa = 6894.75729;
tDays = data.t / 86400;
nL = P.nL;
nV = P.nV;

% Diagnostic: check for size mismatches or complex data
fprintf('Plot diagnostics: t=%d, TL=%dx%d, Tv=%dx%d, Ts=%d, Tw=%d, pv=%d\n', ...
    length(tDays), size(data.TL,1), size(data.TL,2), ...
    size(data.Tv,1), size(data.Tv,2), length(data.Ts), length(data.Tw), ...
    length(data.pv));

% ---- Figure 1: Temperatures --------------------------------------------
figure('Name','Temperatures','NumberTitle','off');

subplot(2,1,1);
hold on;
plot(tDays, real(data.TL(:,nL)), 'b-',  'DisplayName','T_{liquid}');
plot(tDays, real(data.Ts),       'r--', 'DisplayName','T_{film}');
plot(tDays, real(data.Tv(:,nV)), 'Color',[0.93 0.69 0.13], 'DisplayName','T_{vapor}');
plot(tDays, real(data.Tw),       'm-',  'DisplayName','T_{wall}');
hold off;
ylabel('Temperature (K)');
xlabel('Time (days)');
legend('Location','best');
title('Tank Temperatures');
grid on; xlim([0 tDays(end)]);

subplot(2,1,2);
hold on;
for i = 1:nL
    plot(tDays, real(data.TL(:,i)), 'DisplayName', sprintf('TL layer %d', i));
end
hold off;
ylabel('Temperature (K)'); xlabel('Time (days)');
legend('Location','best'); title('Liquid Boundary Layer Temperatures');
grid on; xlim([0 tDays(end)]);

% ---- Figure 2: Pressure and Vent State ----------------------------------
figure('Name','Pressure & Venting','NumberTitle','off');

subplot(2,1,1);
hold on;
plot(tDays, real(data.pv)/psiToPa, 'b-', 'DisplayName','P_v');
yline(P.p_high/psiToPa, 'r--', 'DisplayName','PRD open');
yline(P.p_low/psiToPa, 'r:', 'DisplayName','PRD close');
hold off;
ylabel('Vapor Pressure (psia)'); xlabel('Time (days)');
legend('Location','best');
title('Vapor Pressure');
grid on; xlim([0 tDays(end)]);

subplot(2,1,2);
plot(tDays, real(data.VentState), 'k-');
ylabel('Vent State (0=closed, 1=open)'); xlabel('Time (days)');
title('Vent Valve State');
grid on; xlim([0 tDays(end)]); ylim([-0.1 1.1]);

% ---- Figure 3: Masses and Fill Level ------------------------------------
figure('Name','Mass & Fill Level','NumberTitle','off');

subplot(2,2,1);
plot(tDays, real(data.mL), 'b-');
ylabel('Liquid Mass (kg)'); xlabel('Time (days)');
title('Liquid Mass'); grid on; xlim([0 tDays(end)]);

subplot(2,2,2);
plot(tDays, real(data.mv), 'r-');
ylabel('Vapor Mass (kg)'); xlabel('Time (days)');
title('Vapor Mass'); grid on; xlim([0 tDays(end)]);

subplot(2,2,3);
plot(tDays, real(data.mTotal), 'k-');
ylabel('Total Mass (kg)'); xlabel('Time (days)');
title('Total Mass'); grid on; xlim([0 tDays(end)]);

subplot(2,2,4);
plot(tDays, 100*real(data.pctFill), 'b-');
ylabel('Fill Level (%)'); xlabel('Time (days)');
title('Tank Fill Level'); grid on; xlim([0 tDays(end)]);

% ---- Figure 4: BOG ------------------------------------------------------
figure('Name','Boil-Off Gas','NumberTitle','off');

subplot(2,1,1);
plot(tDays, real(data.BOG_cumulative), 'b-');
ylabel('Cumulative BOG (kg)'); xlabel('Time (days)');
title('Cumulative Boil-Off Loss');
grid on; xlim([0 tDays(end)]);

subplot(2,1,2);
if isfield(data, 'dailyBOG_kg') && ~isempty(data.dailyBOG_kg)
    nDayPlot = length(data.dailyBOG_kg);
    bar(1:nDayPlot, data.dailyBOG_kg, 'FaceColor', [0.3 0.6 0.9]);
    ylabel('Daily BOG (kg/day)'); xlabel('Day');
    title(sprintf('Daily Boil-Off Rate (avg: %.2f kg/day, %.3f%%/day)', ...
        mean(data.dailyBOG_kg), mean(data.dailyBOG_pct)));
    grid on; xlim([0.5 nDayPlot+0.5]);
else
    plot(tDays, real(data.BOG_rate), 'r-');
    ylabel('Vent Flow (kg/s)'); xlabel('Time (days)');
    title('Instantaneous Vent Flow');
    grid on; xlim([0 tDays(end)]);
end

% ---- Figure 5: Energy Balance -------------------------------------------
figure('Name','Energy Balance','NumberTitle','off');

subplot(2,1,1);
hold on;
plot(tDays, real(data.QdotEW),  'b-',  'DisplayName','Q_{EW}');
plot(tDays, real(data.QdotWV),  'r-',  'DisplayName','Q_{WV}');
plot(tDays, real(data.QdotWL),  'Color',[0 0.5 0], 'DisplayName','Q_{WL}');
plot(tDays, real(-data.QdotVS), 'm--', 'DisplayName','-Q_{VS}');
plot(tDays, real(-data.QdotLS), 'c--', 'DisplayName','-Q_{LS}');
hold off;
ylabel('Heat Flow (W)'); xlabel('Time (days)');
legend('Location','best');
title('Heat Flows'); grid on; xlim([0 tDays(end)]);

subplot(2,1,2);
plot(tDays, real(data.Jcd), 'b-');
ylabel('Condensation (kg/s)'); xlabel('Time (days)');
title('Condensation Rate (positive = condensing)');
grid on; xlim([0 tDays(end)]);

% ---- Figure 6: Densities -----------------------------------------------
figure('Name','Densities','NumberTitle','off');

subplot(2,1,1);
plot(tDays, real(data.rho_L), 'b-');
ylabel('Liquid Density (kg/m^3)'); xlabel('Time (days)');
title('Liquid Density'); grid on; xlim([0 tDays(end)]);

subplot(2,1,2);
plot(tDays, real(data.rhov), 'r-');
ylabel('Vapor Density (kg/m^3)'); xlabel('Time (days)');
title('Vapor Density'); grid on; xlim([0 tDays(end)]);

% ---- Summary to command window ------------------------------------------
fprintf('\n===== BOIL-OFF SUMMARY =====\n');
fprintf('Simulation time:    %.1f days\n', tDays(end));
fprintf('Initial liquid:     %.1f kg (%.1f%% fill)\n', ...
    data.mL(1), 100*data.pctFill(1));
fprintf('Final liquid:       %.1f kg (%.1f%% fill)\n', ...
    data.mL(end), 100*data.pctFill(end));
fprintf('Total mass vented:  %.1f kg\n', data.BOG_cumulative(end));
fprintf('Total mass lost:    %.1f kg (%.2f%%)\n', ...
    data.massLoss(end), 100*data.massLoss(end)/data.mTotal(1));
if isfield(data, 'dailyBOG_kg') && ~isempty(data.dailyBOG_kg)
    fprintf('Avg daily BOG:      %.2f kg/day (%.3f%%/day)\n', ...
        mean(data.dailyBOG_kg), mean(data.dailyBOG_pct));
end
fprintf('============================\n\n');

end
