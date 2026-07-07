function plotLH2Data(data)
% plotLH2Data(data)
%   Plots results from LH2Simulate in thesis-quality format.
%   Style matches PLOT_SETTINGS.py: IEEE, LaTeX rendering, Arial font,
%   1200 DPI, inward ticks on all sides, thick spines.
%
%   Produces 2 publication figures (4 panels total):
%     Figure 1: (a) Pressure vs time   (b) Temperature vs time
%     Figure 2: (a) Trailer mass + BOG  (b) Storage tank depletion

try
    P = evalin('base','LH2Model');
catch ME
    if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
        evalin('base','LH2ModelParams');
        P = evalin('base','LH2Model');
    else
        error(ME.message);
    end
end

%% ===== STYLE SETTINGS (matching PLOT_SETTINGS.py) =====
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

% Colors (from PLOT_SETTINGS.py)
c_red    = [228, 26, 28]  / 255;   % '#e41a1c'
c_green  = [0, 128, 0]    / 255;   % '#008000'
c_blue   = [55, 126, 184]  / 255;  % '#377eb8'
c_orange = [255, 127, 0]  / 255;   % '#ff7f00'
c_purple = [152, 78, 163] / 255;   % '#984ea3'
c_grey   = [0.4, 0.4, 0.4];

% Dimensions
fig_width   = 8;        % inches (2-panel figure = 2 x 4)
fig_height  = 3;        % inches
font_name   = 'Arial';
label_fs    = 14;
tick_fs     = 10;
legend_fs   = 8;
lw          = 1;        % line width
spine_lw    = 1.5;
tick_len    = 4;
dpi         = 1200;

% Time in minutes
t_min = data.t / 60;
t_max = 65;  % [min] show fill + brief post-fill settling
t_ticks = 0:10:60;  % x-axis tick marks

%% ===== Helper: apply axis style =====
    function style_ax(ax)
        set(ax, 'FontName', font_name, 'FontSize', tick_fs);
        set(ax, 'LineWidth', spine_lw);
        set(ax, 'TickDir', 'in');
        set(ax, 'TickLength', [tick_len/max(fig_width,fig_height)/72, 0.01]);
        set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on');
        set(ax, 'Box', 'on');
        set(ax, 'Color', 'w');
    end

%% ================================================================
%  FIGURE 1: Pressure and Temperature
%  (a) Tank pressures vs time   (b) Temperatures vs time
%% ================================================================
fig1 = figure('Units', 'inches', 'Position', [1 1 fig_width fig_height], ...
    'PaperUnits', 'inches', 'PaperSize', [fig_width fig_height], ...
    'PaperPosition', [0 0 fig_width fig_height], 'Color', 'w');

% --- (a) Pressure ---
ax1a = subplot(1,2,1);
hold on;
plot(t_min, data.pv1/1e5, '-', 'Color', c_red,  'LineWidth', lw);
plot(t_min, data.pv2/1e5, '-', 'Color', c_blue, 'LineWidth', lw);
plot([0 t_max], [P.p_ET_high/1e5, P.p_ET_high/1e5], '--', 'Color', c_grey, 'LineWidth', 0.75);
plot([0 t_max], [P.p_ET_low/1e5,  P.p_ET_low/1e5],  ':', 'Color', c_grey, 'LineWidth', 0.75);
hold off;
xlim([0 t_max]);
set(ax1a, 'XTick', t_ticks);
xlabel('Time / min', 'FontSize', label_fs, 'FontName', font_name);
ylabel('Pressure / bar', 'FontSize', label_fs, 'FontName', font_name);
lg1a = legend('$P_\mathrm{v}$ (ST)', '$P_\mathrm{v}$ (ET)', ...
    'ET vent (open)', 'ET vent (close)', ...
    'Location', 'east', 'FontSize', legend_fs);
lg1a.BoxFace.ColorType = 'truecoloralpha';
lg1a.BoxFace.ColorData = uint8([255;255;255;230]);
lg1a.EdgeColor = [0 0 0];
lg1a.LineWidth = 0.75;
style_ax(ax1a);
text(0.03, 0.93, '\textbf{(a)}', 'Units', 'normalized', 'FontSize', label_fs, ...
    'Interpreter', 'latex', 'VerticalAlignment', 'top');

% --- (b) Temperature ---
ax1b = subplot(1,2,2);
hold on;
plot(t_min, data.TL1(:,end), '--', 'Color', c_red,    'LineWidth', lw);
plot(t_min, data.Tv1(:,end), '-',  'Color', c_red,    'LineWidth', lw);
plot(t_min, data.TL2(:,end), '--', 'Color', c_blue,   'LineWidth', lw);
plot(t_min, data.Tv2(:,end), '-',  'Color', c_blue,   'LineWidth', lw);
plot(t_min, data.Tw2,        ':',  'Color', c_purple,  'LineWidth', lw);
hold off;
xlim([0 t_max]);
set(ax1b, 'XTick', t_ticks);
xlabel('Time / min', 'FontSize', label_fs, 'FontName', font_name);
ylabel('Temperature / K', 'FontSize', label_fs, 'FontName', font_name);
lg1b = legend('$T_\mathrm{L}$ (ST)', '$T_\mathrm{V}$ (ST)', ...
    '$T_\mathrm{L}$ (ET)', '$T_\mathrm{V}$ (ET)', ...
    '$T_\mathrm{wall}$ (ET)', ...
    'Location', 'east', 'FontSize', legend_fs);
lg1b.BoxFace.ColorType = 'truecoloralpha';
lg1b.BoxFace.ColorData = uint8([255;255;255;230]);
lg1b.EdgeColor = [0 0 0];
lg1b.LineWidth = 0.75;
style_ax(ax1b);
text(0.03, 0.93, '\textbf{(b)}', 'Units', 'normalized', 'FontSize', label_fs, ...
    'Interpreter', 'latex', 'VerticalAlignment', 'top');

% Save
print(fig1, 'figure_transfer_pressure_temp', '-dpng', ['-r' num2str(dpi)]);
fprintf('Saved: figure_transfer_pressure_temp.png\n');

%% ================================================================
%  FIGURE 2: Trailer Mass Balance (single panel)
%  Left axis: total mass + liquid mass in trailer
%  Right axis: cumulative BOG vented from trailer
%% ================================================================
fig2_width = 4;  % single panel, half width
fig2 = figure('Units', 'inches', 'Position', [1 5 fig2_width fig_height], ...
    'PaperUnits', 'inches', 'PaperSize', [fig2_width fig_height], ...
    'PaperPosition', [0 0 fig2_width fig_height], 'Color', 'w');

ax2 = axes;
yyaxis left;
hold on;
plot(t_min, (data.mL2 + data.mv2), '-', 'Color', c_blue, 'LineWidth', lw);
plot(t_min, data.mL2,              '-', 'Color', c_green, 'LineWidth', lw);
hold off;
ylabel('Mass in trailer / kg', 'FontSize', label_fs, 'FontName', font_name);
ax2.YColor = 'k';

yyaxis right;
plot(t_min, data.Boiloff_ET, '-', 'Color', c_red, 'LineWidth', lw);
ylabel('Cumulative BOG / kg', 'FontSize', label_fs, 'FontName', font_name);
ax2.YColor = c_red;

xlim([0 t_max]);
set(ax2, 'XTick', t_ticks);
xlabel('Time / min', 'FontSize', label_fs, 'FontName', font_name);
lg2 = legend('Total mass', 'Liquid only', 'BOG vented', ...
    'Location', 'west', 'FontSize', legend_fs);
lg2.BoxFace.ColorType = 'truecoloralpha';
lg2.BoxFace.ColorData = uint8([255;255;255;230]);
lg2.EdgeColor = [0 0 0];
lg2.LineWidth = 0.75;
style_ax(ax2);

% Save
print(fig2, 'figure_transfer_mass_balance', '-dpng', ['-r' num2str(dpi)]);
fprintf('Saved: figure_transfer_mass_balance.png\n');

end
