function test_plot_style_manager()

set(groot, ...
    'defaultAxesFontSize', 11, ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultAxesLineWidth', 0.5, ...
    'defaultAxesLabelFontSizeMultiplier', 1, ...
    'defaultAxesTitleFontSizeMultiplier', 1);

sty = plot_style_manager();
t = linspace(0, 40, 1200);

% =====================================================================
% Figure 1: real plotting workflow using add_standard_background
% =====================================================================
figure('Name', 'Integrated Background Test', 'Color', 'w');
hold on;

% synthetic signals
u_sp1  = 0.90 + 0.10 * sin(0.35 * t);
u_res1 = 0.90 + 0.08 * sin(0.35 * t - 0.25) + 0.015 * sin(1.8 * t);

u_sp2  = 0.45 + 0.07 * square(0.22 * t);
u_res2 = 0.45 + 0.05 * square(0.22 * t - 0.18) + 0.012 * cos(1.6 * t);

u_sp3  = -0.05 + 0.14 * sin(0.17 * t + 0.9);
u_res3 = -0.05 + 0.12 * sin(0.17 * t + 0.55) + 0.01 * sin(2.3 * t);

u_sp4  = -0.55 + 0.05 * sawtooth(0.24 * t, 0.5);
u_res4 = -0.55 + 0.04 * sawtooth(0.24 * t - 0.1, 0.5) + 0.01 * cos(1.9 * t);

s = plot_style(sty, 'setpoint1', '--', 1.0);
plot(t, u_sp1, s{:}, 'DisplayName', 'setpoint1');
s = plot_style(sty, 'response1', '-', 1.2);
plot(t, u_res1, s{:}, 'DisplayName', 'response1');
s = plot_style(sty, 'setpoint4', ':', 0.9);
plot(t, u_sp2, s{:}, 'DisplayName', 'setpoint4');
s = plot_style(sty, 'response4', '-.', 1.1);
plot(t, u_res2, s{:}, 'DisplayName', 'response4');
s = plot_style(sty, 'setpoint7', '--', 0.8);
plot(t, u_sp3, s{:}, 'DisplayName', 'setpoint7');
s = plot_style(sty, 'response7', '-', 1.0);
plot(t, u_res3, s{:}, 'DisplayName', 'response7');
s = plot_style(sty, 'setpoint10', ':', 0.8);
plot(t, u_sp4, s{:}, 'DisplayName', 'setpoint10');
s = plot_style(sty, 'response10', '-', 1.0);
plot(t, u_res4, s{:}, 'DisplayName', 'response10');

flight_intervals = [ ...
     0   4   0; ...
     4   8   1; ...
     8  11   2; ...
    11  14  10; ...
    14  18  17; ...
    18  22  22; ...
    22  26  18; ...
    26  30  20; ...
    30  33   5; ...
    33  36   3; ...
    36  38   4; ...
    38  40  14];

flight_labels = { ...
    'flight 0', 'flight 1', 'flight 2', 'flight 10', ...
    'flight 17', 'flight 22', 'flight 18', 'flight 20', ...
    'flight 5', 'flight 3', 'flight 4', 'flight 14'};

vtol_intervals = [ ...
     0  12   1; ...
    12  27   2; ...
    27  34   3; ...
    34  40   2];

vtol_labels = {'vtol 1', 'vtol 2', 'vtol 3', 'vtol 2'};

add_standard_background(flight_intervals, flight_labels, true, vtol_intervals, vtol_labels);

xlim([0, 40]);
ylim([-1.05, 1.10]);
grid on;
xlabel('Time (s)');
ylabel('Signal');
title('Integrated test: lines + add_standard_background');
legend('show', 'Location', 'eastoutside');
set(gca, 'Layer', 'top');

% =====================================================================
% Figure 2: all named foreground colors
% =====================================================================
figure('Name', 'Foreground Color Palette', 'Color', 'w');
tiledlayout(2,1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
for k = 1:10
    y = 11 - k;
    s = plot_style(sty, sprintf('setpoint%d', k), '-', 1.5);
    plot([0, 1], [y, y], s{:}, 'DisplayName', sprintf('setpoint%d', k));
end
xlim([0, 1]); ylim([0.5, 10.5]); grid on;
title('setpoint1 ... setpoint10');
legend('show', 'Location', 'eastoutside');

nexttile; hold on;
for k = 1:10
    y = 11 - k;
    s = plot_style(sty, sprintf('response%d', k), '-', 1.5);
    plot([0, 1], [y, y], s{:}, 'DisplayName', sprintf('response%d', k));
end
xlim([0, 1]); ylim([0.5, 10.5]); grid on;
title('response1 ... response10');
legend('show', 'Location', 'eastoutside');

% =====================================================================
% Figure 3: all background color definitions
% =====================================================================
figure('Name', 'Background Color Palette', 'Color', 'w');
tiledlayout(2,1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
flight_vals = [0, 15, 1, 2, 10, 17, 22, 18, 20, 5, 3, 4, 14];
for k = 1:numel(flight_vals)
    v = flight_vals(k);
    x0 = k - 1;
    patch([x0, x0+1, x0+1, x0], [0, 0, 1, 1], sty.color.(sprintf('flight_mode_%d', v)), ...
        'EdgeColor', 'none', 'FaceAlpha', sty.bg.flight_mode.alpha);
    text(x0 + 0.5, 0.5, sprintf('fm %d', v), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 9, 'Interpreter', 'none');
end
xlim([0, numel(flight_vals)]); ylim([0, 1]);
set(gca, 'YTick', [], 'XTick', []);
title('flight_mode background colors');
box on;

nexttile; hold on;
vtol_vals = [1, 2, 3];
for k = 1:numel(vtol_vals)
    v = vtol_vals(k);
    x0 = k - 1;
    patch([x0, x0+1, x0+1, x0], [0, 0, 1, 1], sty.color.(sprintf('vtol_state_%d', v)), ...
        'EdgeColor', 'none', 'FaceAlpha', sty.bg.vtol_state.alpha);
    text(x0 + 0.5, 0.5, sprintf('vtol %d', v), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 9, 'Interpreter', 'none');
end
xlim([0, numel(vtol_vals)]); ylim([0, 1]);
set(gca, 'YTick', [], 'XTick', []);
title('vtol_state background colors');
box on;

end
