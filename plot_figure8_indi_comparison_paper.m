clear; close all; clc;

%% Figure-eight INDI comparison for the paper
% The script deliberately compares every controller against the original
% trajectory_setpoint. vehicle_local_position_setpoint is a controller-internal
% corrected setpoint and is therefore not a common reference for comparing
% different outer-loop controllers.

%% User settings
root_dir = fileparts(mfilename('fullpath'));
mat_file = fullfile(root_dir, 'data', ...
    'log_103_2026-8-12-18-12-08.mat');
output_dir = fullfile(root_dir, 'results', 'figure8_indi_comparison');

orbit_nav_state = 21;
number_of_periods = 1;
fallback_omega = 0.25;       % rad/s, only used if automatic detection fails
phase_samples = 1201;
assert(number_of_periods == 1, ...
    'Phase-aligned paper plots currently require number_of_periods = 1.');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Load and validate data
S = load(mat_file, 'log');
D = S.log.data;
required_topics = { ...
    'trajectory_setpoint_0', 'vehicle_local_position_0', ...
    'vehicle_attitude_0', 'vehicle_attitude_setpoint_0', ...
    'vehicle_angular_velocity_0', 'vehicle_rates_setpoint_0', ...
    'vehicle_local_position_setpoint_0', 'rate_ctrl_status_0', ...
    'vehicle_status_0'};

for k = 1:numel(required_topics)
    assert(isfield(D, required_topics{k}), ...
        'Required topic is missing: %s', required_topics{k});
end

traj_tbl = D.trajectory_setpoint_0;
local_tbl = D.vehicle_local_position_0;
att_tbl = D.vehicle_attitude_0;
att_sp_tbl = D.vehicle_attitude_setpoint_0;
rate_tbl = D.vehicle_angular_velocity_0;
rate_sp_tbl = D.vehicle_rates_setpoint_0;
local_sp_tbl = D.vehicle_local_position_setpoint_0;
rate_status_tbl = D.rate_ctrl_status_0;
status_tbl = D.vehicle_status_0;

t_traj = topic_time(traj_tbl);
t_local = topic_time(local_tbl);
t_att = topic_time(att_tbl);
t_att_sp = topic_time(att_sp_tbl);
t_rate = topic_time(rate_tbl);
t_rate_sp = topic_time(rate_sp_tbl);
t_local_sp = topic_time(local_sp_tbl);
t_rate_status = topic_time(rate_status_tbl);
t_status = topic_time(status_tbl);

position_ref_all = double([traj_tbl.position_0_, traj_tbl.position_1_, traj_tbl.position_2_]);
velocity_ref_all = double([traj_tbl.velocity_0_, traj_tbl.velocity_1_, traj_tbl.velocity_2_]);
position_all = double([local_tbl.x, local_tbl.y, local_tbl.z]);
velocity_all = double([local_tbl.vx, local_tbl.vy, local_tbl.vz]);
q_all = make_quaternion_continuous(double([att_tbl.q_0_, att_tbl.q_1_, att_tbl.q_2_, att_tbl.q_3_]));
q_sp_all = make_quaternion_continuous(double([att_sp_tbl.q_d_0_, att_sp_tbl.q_d_1_, ...
    att_sp_tbl.q_d_2_, att_sp_tbl.q_d_3_]));
rate_all = double([rate_tbl.xyz_0_, rate_tbl.xyz_1_, rate_tbl.xyz_2_]);
rate_sp_all = double([rate_sp_tbl.roll, rate_sp_tbl.pitch, rate_sp_tbl.yaw]);

%% Find the three formal figure-eight activations
% log103 contains an initial low-speed circle, three formal circles, an
% initial low-speed figure eight, and then the three formal figure-eight
% comparisons. Taking the final three Orbit intervals excludes both initial
% low-speed trials without relying on manually typed timestamps.
is_orbit = double(status_tbl.nav_state) == orbit_nav_state;
start_idx = find(diff([false; is_orbit]) == 1);
end_idx = find(diff([is_orbit; false]) == -1);
orbit_intervals = [t_status(start_idx), t_status(end_idx)];
assert(size(orbit_intervals, 1) >= 3, ...
    'Fewer than three Orbit intervals were found.');
formal_intervals = orbit_intervals(end-2:end, :);

%% Recover the actual post-adjustment period from trajectory_setpoint
all_period_samples = [];

for k = 1:3
    crossings = positive_center_crossings(t_traj, position_ref_all(:, 1), ...
        velocity_ref_all(:, 1), formal_intervals(k, :));

    if numel(crossings) >= 2
        all_period_samples = [all_period_samples; diff(crossings)]; %#ok<AGROW>
    end
end

if isempty(all_period_samples)
    period = 2 * pi / fallback_omega;
    warning('Could not recover the period; using fallback omega %.3f rad/s.', fallback_omega);
else
    period = median(all_period_samples, 'omitnan');
end

omega = 2 * pi / period;
window_duration = number_of_periods * period;
fprintf('Recovered formal figure-eight omega: %.5f rad/s\n', omega);
fprintf('Recovered period: %.5f s\n', period);

%% Select the same phase-aligned middle period from each activation
chronological_cases = repmat(empty_case(), 1, 3);

for k = 1:3
    interval = formal_intervals(k, :);
    crossings = positive_center_crossings(t_traj, position_ref_all(:, 1), ...
        velocity_ref_all(:, 1), interval);
    start_time = choose_centered_complete_window(interval, window_duration);
    time_grid = linspace(start_time, start_time + window_duration, phase_samples).';
    assert(~isempty(crossings), 'No reference phase crossing was found in a formal interval.');
    [~, phase_origin_idx] = min(abs(crossings - mean(time_grid)));
    phase_origin = crossings(phase_origin_idx);
    raw_phase = mod((time_grid - phase_origin) * omega, 2 * pi);
    [phase, plot_order] = sort(raw_phase);

    c = empty_case();
    c.interval = interval;
    c.window = [time_grid(1), time_grid(end)];
    c.time = time_grid;
    c.phase = phase;
    c.plot_time = phase / omega;
    c.plot_order = plot_order;
    c.position_ref = interp_series(t_traj, position_ref_all, time_grid);
    c.velocity_ref = interp_series(t_traj, velocity_ref_all, time_grid);
    c.position = interp_series(t_local, position_all, time_grid);
    c.velocity = interp_series(t_local, velocity_all, time_grid);
    c.q_sp = interp_quaternion(t_att_sp, q_sp_all, time_grid);
    c.q = interp_quaternion(t_att, q_all, time_grid);
    c.rate_sp = interp_series(t_rate_sp, rate_sp_all, time_grid);
    c.rate = interp_series(t_rate, rate_all, time_grid);

    c.acc_indi_fraction = logical_fraction(t_local_sp, ...
        double(local_sp_tbl.acc_indi_active), c.window);
    c.rate_indi_fraction = logical_fraction(t_rate_status, ...
        double(rate_status_tbl.indi_active), c.window);
    c.name = identify_algorithm(c.acc_indi_fraction, c.rate_indi_fraction);

    % Error notation follows main.tex:
    %   e_p = p_r - p, e_v = v_r - v,
    %   e_R = R_c \ominus R, e_omega = omega_r^b - omega^b.
    c.position_error = c.position_ref - c.position;
    c.velocity_error = c.velocity_ref - c.velocity;
    c.position_error_xy = vecnorm(c.position_error(:, 1:2), 2, 2);
    c.velocity_error_xy = vecnorm(c.velocity_error(:, 1:2), 2, 2);

    c.height_error = c.position_error(:, 3);
    c.vertical_velocity_error = c.velocity_error(:, 3);

    q_dot = abs(sum(c.q_sp .* c.q, 2));
    c.attitude_error = 2 * acos(constrain(q_dot, 0, 1));
    c.rate_error = vecnorm(c.rate_sp - c.rate, 2, 2);

    c.metrics = calculate_metrics(c);
    chronological_cases(k) = c;
end

% Always present algorithms in the conceptual order used in the paper.
desired_order = {'Baseline', 'PINDI', 'GINDI'};
cases = repmat(empty_case(), 1, 3);

for k = 1:3
    idx = find(strcmp({chronological_cases.name}, desired_order{k}), 1);
    assert(~isempty(idx), 'Could not identify the %s interval.', desired_order{k});
    cases(k) = chronological_cases(idx);
end

%% Quantitative summary
algorithm = string({cases.name}).';
xy_position_rmse_m = arrayfun(@(c) c.metrics.xy_position_rmse, cases).';
xy_velocity_rmse_mps = arrayfun(@(c) c.metrics.xy_velocity_rmse, cases).';
height_rmse_m = arrayfun(@(c) c.metrics.height_rmse, cases).';
vertical_velocity_rmse_mps = arrayfun(@(c) c.metrics.vertical_velocity_rmse, cases).';
attitude_rmse_deg = arrayfun(@(c) c.metrics.attitude_rmse_deg, cases).';
angular_rate_rmse_radps = arrayfun(@(c) c.metrics.rate_rmse, cases).';
baseline = [xy_position_rmse_m(1), xy_velocity_rmse_mps(1), ...
    height_rmse_m(1), vertical_velocity_rmse_mps(1), ...
    attitude_rmse_deg(1), angular_rate_rmse_radps(1)];
values = [xy_position_rmse_m, xy_velocity_rmse_mps, height_rmse_m, ...
    vertical_velocity_rmse_mps, attitude_rmse_deg, angular_rate_rmse_radps];
improvement_percent = 100 * (1 - values ./ baseline);

metrics_table = table(algorithm, xy_position_rmse_m, xy_velocity_rmse_mps, ...
    height_rmse_m, vertical_velocity_rmse_mps, attitude_rmse_deg, ...
    angular_rate_rmse_radps, ...
    improvement_percent(:, 1), improvement_percent(:, 2), ...
    improvement_percent(:, 3), improvement_percent(:, 4), ...
    improvement_percent(:, 5), improvement_percent(:, 6), ...
    'VariableNames', {'Algorithm', 'XYPositionRMSE_m', 'XYVelocityRMSE_mps', ...
    'HeightRMSE_m', 'VerticalVelocityRMSE_mps', 'AttitudeRMSE_deg', ...
    'AngularRateRMSE_radps', ...
    'XYPositionImprovement_pct', 'XYVelocityImprovement_pct', ...
    'HeightImprovement_pct', 'VerticalVelocityImprovement_pct', ...
    'AttitudeImprovement_pct', 'AngularRateImprovement_pct'});

disp(metrics_table);
writetable(metrics_table, fullfile(output_dir, 'figure8_indi_comparison_metrics.csv'));

for k = 1:3
    fprintf('%-10s window %.3f--%.3f s, acc INDI %.1f%%, rate INDI %.1f%%\n', ...
        cases(k).name, cases(k).window(1), cases(k).window(2), ...
        100 * cases(k).acc_indi_fraction, 100 * cases(k).rate_indi_fraction);
end

%% Paper style
addpath(root_dir);
set(groot, ...
    'defaultAxesFontSize', 6, ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultAxesLineWidth', 0.5, ...
    'defaultAxesLabelFontSizeMultiplier', 1, ...
    'defaultAxesTitleFontSizeMultiplier', 1, ...
    'defaultTextFontName', 'Times New Roman', ...
    'defaultLegendFontName', 'Times New Roman', ...
    'defaultTextInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex');

% Match the comparison styles in plot_case2_for_GINDI.m exactly:
%   algorithm 1: purple dashed,       0.5 pt
%   algorithm 2: orange dash-dot,     0.5 pt
%   algorithm 3: blue solid,          0.7 pt
% The common reference follows that file's setpoint style.
sty = plot_style_manager();
STYLE_REFERENCE = plot_style(sty, 'setpoint1', '--', 0.8);
STYLE_BASELINE = {'Color', [0.494, 0.184, 0.556], ...
    'LineStyle', '--', 'LineWidth', 0.5};
STYLE_RATE_INDI = {'Color', [0.96, 0.62, 0.26], ...
    'LineStyle', '-.', 'LineWidth', 0.5};
STYLE_DUAL_INDI = {'Color', [0.000, 0.447, 0.741], ...
    'LineStyle', '-', 'LineWidth', 0.7};
case_plot_styles = {STYLE_BASELINE, STYLE_RATE_INDI, STYLE_DUAL_INDI};

%% Figure 1: phase-aligned XY trajectories as three LaTeX subfigures
trajectory_reference = cell(1, 3);
trajectory_response = cell(1, 3);
trajectory_file_stems = {'baseline', 'pindi', 'gindi'};

for k = 1:3
    order = cases(k).plot_order;
    trajectory_reference{k} = cases(k).position_ref(order, 1:2);
    trajectory_response{k} = cases(k).position(order, 1:2);
end

for k = 1:3
    fig_xy = figure('Name', sprintf('Figure-eight XY %s', cases(k).name), ...
        'Color', 'w', 'Units', 'centimeters', 'Position', [0, 0, 5.0, 5.0], ...
        'PaperUnits', 'centimeters', 'PaperSize', [5.0, 5.0], ...
        'PaperPosition', [0, 0, 5.0, 5.0], 'PaperPositionMode', 'auto', ...
        'InvertHardcopy', 'off');
    ax = axes(fig_xy); hold(ax, 'on');
    h_ref = plot(ax, trajectory_reference{k}(:, 1), trajectory_reference{k}(:, 2), ...
        STYLE_REFERENCE{:});
    h_res = plot(ax, trajectory_response{k}(:, 1), trajectory_response{k}(:, 2), ...
        case_plot_styles{k}{:});
    plot(ax, trajectory_reference{k}(1, 1), trajectory_reference{k}(1, 2), 'o', ...
        'Color', sty.color.setpoint1, 'MarkerSize', 3, 'HandleVisibility', 'off');
    axis(ax, 'equal'); grid(ax, 'on'); box(ax, 'on');
    xy = [trajectory_reference{k}; trajectory_response{k}];
    xy_min = min(xy, [], 1);
    xy_max = max(xy, [], 1);
    xy_pad = 0.06 * max(xy_max - xy_min);
    xlim(ax, [xy_min(1) - xy_pad, xy_max(1) + xy_pad]);
    ylim(ax, [xy_min(2) - xy_pad, xy_max(2) + xy_pad]);
    xlabel(ax, '$p_x$ (m)');
    ylabel(ax, '$p_y$ (m)');
    lgd = legend(ax, [h_ref, h_res], {'Reference', 'Response'}, ...
        'Location', 'north', 'NumColumns', 1);
    lgd.ItemTokenSize = [20, 7];

    PlotToFile(fig_xy, fullfile(output_dir, sprintf( ...
        'figure8_xy_tracking_%s.pdf', trajectory_file_stems{k})), 5.2, 5.0);
end

%% Figure 2: translational tracking errors versus time
fig_trans = figure('Name', 'Figure-eight translation errors', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [0, 0, 7.0, 9.0], ...
    'PaperUnits', 'centimeters', 'PaperSize', [7.0, 9.0], ...
    'PaperPosition', [0, 0, 7.0, 9.0], 'PaperPositionMode', 'auto', ...
    'InvertHardcopy', 'off');
tl_trans = tiledlayout(fig_trans, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
trans_data = { ...
    % Vector-valued horizontal errors are shown by their Euclidean norms.
    @(c) c.position_error_xy, ...
    @(c) c.velocity_error_xy, ...
    % Scalar vertical errors retain their signs to expose bias and direction.
    @(c) 100 * c.height_error, ...
    @(c) 100 * c.vertical_velocity_error};
trans_labels = {'$||\mathbf{e}_{p,xy}||$ (m)', '$||\mathbf{e}_{v,xy}||$ (m/s)', ...
    '$e_{p,z}$ (cm)', '$e_{v,z}$ (cm/s)'};
trans_error_is_signed = [false, false, true, true];
ax_trans = gobjects(1, 4);

for row = 1:4
    ax_trans(row) = nexttile(tl_trans); hold(ax_trans(row), 'on');
    for k = 1:3
        y = trans_data{row}(cases(k));
        plot(ax_trans(row), cases(k).plot_time, y(cases(k).plot_order), ...
            case_plot_styles{k}{:});
    end
    if trans_error_is_signed(row)
        yline(ax_trans(row), 0, '-', 'Color', [0.75, 0.75, 0.75], ...
            'LineWidth', 0.35, 'HandleVisibility', 'off');
    end
    ylabel(ax_trans(row), trans_labels{row});
    xlim(ax_trans(row), [0, window_duration]);
    xticks(ax_trans(row), 0:5:window_duration);
    grid(ax_trans(row), 'on'); box(ax_trans(row), 'on');
end

lgd = legend(ax_trans(1), {cases.name}, 'Location', 'north', 'NumColumns', 3);
lgd.ItemTokenSize = [20, 7];
xlabel(ax_trans(end), 'Time (s)');
linkaxes(ax_trans, 'x');
PlotToFile(fig_trans, fullfile(output_dir, 'figure8_translation_errors.pdf'), 7.0, 9.0);

%% Figure 3: rotational inner-loop tracking errors
fig_rot = figure('Name', 'Figure-eight rotational errors', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [0, 0, 7.0, 5.5], ...
    'PaperUnits', 'centimeters', 'PaperSize', [7.0, 5.5], ...
    'PaperPosition', [0, 0, 7.0, 5.5], 'PaperPositionMode', 'auto', ...
    'InvertHardcopy', 'off');
tl_rot = tiledlayout(fig_rot, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
% Both rotational errors are vector-valued and are shown by their norms.
rot_data = {@(c) c.attitude_error * 180 / pi, @(c) c.rate_error};
rot_labels = {'$||\mathbf{e}_{R}||$ (deg)', '$||\mathbf{e}_{\omega}||$ (rad/s)'};
ax_rot = gobjects(1, 2);

for row = 1:2
    ax_rot(row) = nexttile(tl_rot); hold(ax_rot(row), 'on');
    for k = 1:3
        y = rot_data{row}(cases(k));
        plot(ax_rot(row), cases(k).plot_time, y(cases(k).plot_order), ...
            case_plot_styles{k}{:});
    end
    ylabel(ax_rot(row), rot_labels{row});
    xlim(ax_rot(row), [0, window_duration]);
    xticks(ax_rot(row), 0:5:window_duration);
    grid(ax_rot(row), 'on'); box(ax_rot(row), 'on');
end

lgd = legend(ax_rot(1), {cases.name}, 'Location', 'north', 'NumColumns', 3);
lgd.ItemTokenSize = [20, 7];
xlabel(ax_rot(end), 'Time (s)');
linkaxes(ax_rot, 'x');
PlotToFile(fig_rot, fullfile(output_dir, 'figure8_rotational_errors.pdf'), 7.0, 5.5);

fprintf('Saved paper outputs to:\n  %s\n', output_dir);

%% Local helpers
function c = empty_case()
    c = struct('name', '', 'interval', [NaN, NaN], 'window', [NaN, NaN], ...
        'time', [], 'phase', [], 'plot_time', [], 'plot_order', [], ...
        'position_ref', [], 'velocity_ref', [], ...
        'position', [], 'velocity', [], 'q_sp', [], 'q', [], ...
        'rate_sp', [], 'rate', [], 'position_error', [], ...
        'velocity_error', [], 'position_error_xy', [], ...
        'velocity_error_xy', [], 'height_error', [], ...
        'vertical_velocity_error', [], 'attitude_error', [], 'rate_error', [], ...
        'acc_indi_fraction', NaN, 'rate_indi_fraction', NaN, ...
        'metrics', struct());
end

function t = topic_time(tbl)
    t = double(tbl.timestamp) * 1e-6;
end

function yq = interp_series(t, y, tq)
    valid = isfinite(t) & all(isfinite(y), 2);
    [tu, ia] = unique(t(valid), 'stable');
    yu = y(valid, :);
    yu = yu(ia, :);
    yq = interp1(tu, yu, tq, 'linear', NaN);
end

function q = normalize_quaternion(q)
    q = q ./ max(vecnorm(q, 2, 2), eps);
end

function q = make_quaternion_continuous(q)
    q = normalize_quaternion(q);
    for k = 2:size(q, 1)
        if dot(q(k, :), q(k-1, :)) < 0
            q(k, :) = -q(k, :);
        end
    end
end

function qq = interp_quaternion(t, q, tq)
    q = make_quaternion_continuous(q);
    qq = interp_series(t, q, tq);
    qq = normalize_quaternion(qq);
end

function crossings = positive_center_crossings(t, x, vx, interval)
    use = t >= interval(1) & t <= interval(2) & ...
        isfinite(x) & isfinite(vx);
    tt = t(use); xx = x(use); vv = vx(use);
    if numel(tt) < 3
        crossings = [];
        return;
    end

    center = 0.5 * (max(xx) + min(xx));
    centered_x = xx - center;
    idx = find(centered_x(1:end-1) <= 0 & centered_x(2:end) > 0 & vv(2:end) > 0);
    crossings = zeros(numel(idx), 1);

    for k = 1:numel(idx)
        i = idx(k);
        ratio = -centered_x(i) / (centered_x(i+1) - centered_x(i));
        crossings(k) = tt(i) + ratio * (tt(i+1) - tt(i));
    end
end

function start_time = choose_centered_complete_window(interval, duration)
    assert(diff(interval) >= duration, ...
        'No complete trajectory period fits inside interval %.3f--%.3f s.', ...
        interval(1), interval(2));
    start_time = mean(interval) - duration / 2;
end

function fraction = logical_fraction(t, x, window)
    use = t >= window(1) & t <= window(2) & isfinite(x);
    assert(any(use), 'No mode-status samples in the selected interval.');
    fraction = mean(x(use) > 0.5);
end

function name = identify_algorithm(acc_fraction, rate_fraction)
    if acc_fraction > 0.9 && rate_fraction > 0.9
        name = 'GINDI';
    elseif acc_fraction < 0.1 && rate_fraction > 0.9
        name = 'PINDI';
    elseif acc_fraction < 0.1 && rate_fraction < 0.1
        name = 'Baseline';
    else
        error('Controller mode is not constant in the selected period: acc=%.2f, rate=%.2f.', ...
            acc_fraction, rate_fraction);
    end
end

function m = calculate_metrics(c)
    m.xy_position_rmse = sqrt(mean(c.position_error_xy.^2, 'omitnan'));
    m.xy_velocity_rmse = sqrt(mean(c.velocity_error_xy.^2, 'omitnan'));
    m.height_rmse = sqrt(mean(c.height_error.^2, 'omitnan'));
    m.vertical_velocity_rmse = sqrt(mean(c.vertical_velocity_error.^2, 'omitnan'));
    m.attitude_rmse_deg = sqrt(mean((c.attitude_error * 180 / pi).^2, 'omitnan'));
    m.rate_rmse = sqrt(mean(c.rate_error.^2, 'omitnan'));
end


function y = constrain(x, lower, upper)
    y = min(max(x, lower), upper);
end
