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
lag_search_limit = 1.0;      % s
lag_search_step = 0.005;     % s
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

    c.position_error = c.position - c.position_ref;
    c.velocity_error = c.velocity - c.velocity_ref;
    c.position_error_xy = vecnorm(c.position_error(:, 1:2), 2, 2);
    c.velocity_error_xy = vecnorm(c.velocity_error(:, 1:2), 2, 2);

    % Positive height error means the aircraft is below the requested height.
    % In NED coordinates h_ref-h = z-z_ref.
    c.height_error = c.position_error(:, 3);
    c.vertical_velocity_error = c.velocity_error(:, 3);

    tangent = c.velocity_ref(:, 1:2);
    tangent = tangent ./ max(vecnorm(tangent, 2, 2), eps);
    normal = [-tangent(:, 2), tangent(:, 1)];
    c.along_error = sum(c.position_error(:, 1:2) .* tangent, 2);
    c.cross_error = sum(c.position_error(:, 1:2) .* normal, 2);

    q_dot = abs(sum(c.q_sp .* c.q, 2));
    c.attitude_error = 2 * acos(constrain(q_dot, 0, 1));
    c.rate_error = vecnorm(c.rate_sp - c.rate, 2, 2);

    [c.equivalent_delay, c.delay_corrected_xy_rmse] = estimate_xy_delay( ...
        t_traj, position_ref_all(:, 1:2), time_grid, c.position(:, 1:2), ...
        lag_search_limit, lag_search_step);

    c.metrics = calculate_metrics(c);
    chronological_cases(k) = c;
end

% Always present algorithms in the conceptual order used in the paper.
desired_order = {'Baseline', 'Rate INDI', 'Dual INDI'};
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
along_track_rmse_m = arrayfun(@(c) c.metrics.along_rmse, cases).';
cross_track_rmse_m = arrayfun(@(c) c.metrics.cross_rmse, cases).';
equivalent_delay_s = arrayfun(@(c) c.equivalent_delay, cases).';
delay_corrected_xy_rmse_m = arrayfun(@(c) c.delay_corrected_xy_rmse, cases).';

baseline = [xy_position_rmse_m(1), xy_velocity_rmse_mps(1), ...
    height_rmse_m(1), vertical_velocity_rmse_mps(1), ...
    attitude_rmse_deg(1), angular_rate_rmse_radps(1)];
values = [xy_position_rmse_m, xy_velocity_rmse_mps, height_rmse_m, ...
    vertical_velocity_rmse_mps, attitude_rmse_deg, angular_rate_rmse_radps];
improvement_percent = 100 * (1 - values ./ baseline);

metrics_table = table(algorithm, xy_position_rmse_m, xy_velocity_rmse_mps, ...
    height_rmse_m, vertical_velocity_rmse_mps, attitude_rmse_deg, ...
    angular_rate_rmse_radps, along_track_rmse_m, cross_track_rmse_m, ...
    equivalent_delay_s, delay_corrected_xy_rmse_m, ...
    improvement_percent(:, 1), improvement_percent(:, 2), ...
    improvement_percent(:, 3), improvement_percent(:, 4), ...
    improvement_percent(:, 5), improvement_percent(:, 6), ...
    'VariableNames', {'Algorithm', 'XYPositionRMSE_m', 'XYVelocityRMSE_mps', ...
    'HeightRMSE_m', 'VerticalVelocityRMSE_mps', 'AttitudeRMSE_deg', ...
    'AngularRateRMSE_radps', 'AlongTrackRMSE_m', 'CrossTrackRMSE_m', ...
    'EquivalentDelay_s', 'DelayCorrectedXYRMSE_m', ...
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
    'defaultAxesFontSize', 7, ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultAxesLineWidth', 0.4, ...
    'defaultAxesLabelFontSizeMultiplier', 1, ...
    'defaultAxesTitleFontSizeMultiplier', 1, ...
    'defaultTextFontName', 'Times New Roman', ...
    'defaultLegendFontName', 'Times New Roman', ...
    'defaultTextInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex');

COLORS.reference = [0.15, 0.15, 0.15];
COLORS.baseline = [0.45, 0.45, 0.45];
COLORS.rate = [0.850, 0.325, 0.098];
COLORS.dual = [0.000, 0.447, 0.741];
case_colors = [COLORS.baseline; COLORS.rate; COLORS.dual];
case_styles = {':', '--', '-'};
line_widths = [0.95, 0.95, 1.15];

%% Figure 1: phase-aligned XY trajectory small multiples
relative_reference = cell(1, 3);
relative_response = cell(1, 3);
all_xy = [];

for k = 1:3
    center = 0.5 * (max(cases(k).position_ref(:, 1:2), [], 1) + ...
        min(cases(k).position_ref(:, 1:2), [], 1));
    order = cases(k).plot_order;
    relative_reference{k} = cases(k).position_ref(order, 1:2) - center;
    relative_response{k} = cases(k).position(order, 1:2) - center;
    all_xy = [all_xy; relative_reference{k}; relative_response{k}]; %#ok<AGROW>
end

xy_min = min(all_xy, [], 1);
xy_max = max(all_xy, [], 1);
xy_pad = 0.06 * max(xy_max - xy_min);

fig_xy = figure('Name', 'Figure-eight XY comparison', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [0, 0, 17.8, 5.8]);
tl_xy = tiledlayout(fig_xy, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:3
    ax = nexttile(tl_xy); hold(ax, 'on');
    h_ref = plot(ax, relative_reference{k}(:, 1), relative_reference{k}(:, 2), ...
        '--', 'Color', COLORS.reference, 'LineWidth', 0.85);
    h_res = plot(ax, relative_response{k}(:, 1), relative_response{k}(:, 2), ...
        case_styles{k}, 'Color', case_colors(k, :), 'LineWidth', line_widths(k));
    plot(ax, relative_reference{k}(1, 1), relative_reference{k}(1, 2), 'o', ...
        'Color', COLORS.reference, 'MarkerSize', 3, 'HandleVisibility', 'off');
    axis(ax, 'equal'); grid(ax, 'on'); box(ax, 'on');
    xlim(ax, [xy_min(1) - xy_pad, xy_max(1) + xy_pad]);
    ylim(ax, [xy_min(2) - xy_pad, xy_max(2) + xy_pad]);
    xlabel(ax, '$p_x-p_{x,c}$ (m)');
    if k == 1
        ylabel(ax, '$p_y-p_{y,c}$ (m)');
    end
    title(ax, cases(k).name);
    text(ax, 0.03, 0.05, sprintf('$e_{p,xy}=%.3f$ m', ...
        cases(k).metrics.xy_position_rmse), 'Units', 'normalized', ...
        'VerticalAlignment', 'bottom');

    if k == 1
        lgd = legend(ax, [h_ref, h_res], {'Reference', 'Response'}, ...
            'Location', 'northoutside', 'NumColumns', 2);
        lgd.ItemTokenSize = [13, 7];
    end
end

PlotToFile(fig_xy, fullfile(output_dir, 'figure8_xy_tracking_comparison.pdf'), 17.8, 5.8);

%% Figure 2: translational tracking errors versus trajectory phase
fig_trans = figure('Name', 'Figure-eight translation errors', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [0, 0, 12.0, 11.0]);
tl_trans = tiledlayout(fig_trans, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
trans_data = { ...
    @(c) c.position_error_xy, ...
    @(c) c.velocity_error_xy, ...
    @(c) 100 * c.height_error, ...
    @(c) 100 * c.vertical_velocity_error};
trans_labels = {'Horizontal position error (m)', 'Horizontal velocity error (m/s)', ...
    '$e_h$ (cm)', '$e_{v_h}$ (cm/s)'};
ax_trans = gobjects(1, 4);

for row = 1:4
    ax_trans(row) = nexttile(tl_trans); hold(ax_trans(row), 'on');
    for k = 1:3
        y = trans_data{row}(cases(k));
        plot(ax_trans(row), cases(k).phase * 180 / pi, y(cases(k).plot_order), ...
            case_styles{k}, 'Color', case_colors(k, :), 'LineWidth', line_widths(k));
    end
    yline(ax_trans(row), 0, '-', 'Color', [0.75, 0.75, 0.75], ...
        'LineWidth', 0.35, 'HandleVisibility', 'off');
    ylabel(ax_trans(row), trans_labels{row});
    xlim(ax_trans(row), [0, 360 * number_of_periods]);
    xticks(ax_trans(row), 0:120:360 * number_of_periods);
    grid(ax_trans(row), 'on'); box(ax_trans(row), 'on');
end

lgd = legend(ax_trans(1), {cases.name}, 'Location', 'northoutside', 'NumColumns', 3);
lgd.ItemTokenSize = [15, 7];
xlabel(ax_trans(end), 'Trajectory phase $\theta$ (deg)');
linkaxes(ax_trans, 'x');
PlotToFile(fig_trans, fullfile(output_dir, 'figure8_translation_errors.pdf'), 12.0, 11.0);

%% Figure 3: along-track and cross-track position errors
fig_path = figure('Name', 'Figure-eight path error decomposition', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [0, 0, 12.0, 6.4]);
tl_path = tiledlayout(fig_path, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
path_data = {@(c) c.along_error, @(c) c.cross_error};
path_labels = {'Along-track error $e_{\parallel}$ (m)', ...
    'Cross-track error $e_{\perp}$ (m)'};
ax_path = gobjects(1, 2);

for row = 1:2
    ax_path(row) = nexttile(tl_path); hold(ax_path(row), 'on');
    for k = 1:3
        y = path_data{row}(cases(k));
        plot(ax_path(row), cases(k).phase * 180 / pi, y(cases(k).plot_order), ...
            case_styles{k}, 'Color', case_colors(k, :), 'LineWidth', line_widths(k));
    end
    yline(ax_path(row), 0, '-', 'Color', [0.75, 0.75, 0.75], ...
        'LineWidth', 0.35, 'HandleVisibility', 'off');
    ylabel(ax_path(row), path_labels{row});
    xlim(ax_path(row), [0, 360 * number_of_periods]);
    xticks(ax_path(row), 0:120:360 * number_of_periods);
    grid(ax_path(row), 'on'); box(ax_path(row), 'on');
end

lgd = legend(ax_path(1), {cases.name}, 'Location', 'northoutside', 'NumColumns', 3);
lgd.ItemTokenSize = [15, 7];
xlabel(ax_path(end), 'Trajectory phase $\theta$ (deg)');
linkaxes(ax_path, 'x');
PlotToFile(fig_path, fullfile(output_dir, 'figure8_path_error_decomposition.pdf'), 12.0, 6.4);

%% Figure 4: rotational inner-loop tracking errors
fig_rot = figure('Name', 'Figure-eight rotational errors', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [0, 0, 12.0, 6.4]);
tl_rot = tiledlayout(fig_rot, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
rot_data = {@(c) c.attitude_error * 180 / pi, @(c) c.rate_error};
rot_labels = {'Attitude geodesic error $e_R$ (deg)', ...
    'Angular-rate error (rad/s)'};
ax_rot = gobjects(1, 2);

for row = 1:2
    ax_rot(row) = nexttile(tl_rot); hold(ax_rot(row), 'on');
    for k = 1:3
        y = rot_data{row}(cases(k));
        plot(ax_rot(row), cases(k).phase * 180 / pi, y(cases(k).plot_order), ...
            case_styles{k}, 'Color', case_colors(k, :), 'LineWidth', line_widths(k));
    end
    ylabel(ax_rot(row), rot_labels{row});
    xlim(ax_rot(row), [0, 360 * number_of_periods]);
    xticks(ax_rot(row), 0:120:360 * number_of_periods);
    grid(ax_rot(row), 'on'); box(ax_rot(row), 'on');
end

lgd = legend(ax_rot(1), {cases.name}, 'Location', 'northoutside', 'NumColumns', 3);
lgd.ItemTokenSize = [15, 7];
xlabel(ax_rot(end), 'Trajectory phase $\theta$ (deg)');
linkaxes(ax_rot, 'x');
PlotToFile(fig_rot, fullfile(output_dir, 'figure8_rotational_errors.pdf'), 12.0, 6.4);

fprintf('Saved paper outputs to:\n  %s\n', output_dir);

%% Local helpers
function c = empty_case()
    c = struct('name', '', 'interval', [NaN, NaN], 'window', [NaN, NaN], ...
        'time', [], 'phase', [], 'plot_order', [], 'position_ref', [], 'velocity_ref', [], ...
        'position', [], 'velocity', [], 'q_sp', [], 'q', [], ...
        'rate_sp', [], 'rate', [], 'position_error', [], ...
        'velocity_error', [], 'position_error_xy', [], ...
        'velocity_error_xy', [], 'height_error', [], ...
        'vertical_velocity_error', [], 'along_error', [], ...
        'cross_error', [], 'attitude_error', [], 'rate_error', [], ...
        'acc_indi_fraction', NaN, 'rate_indi_fraction', NaN, ...
        'equivalent_delay', NaN, 'delay_corrected_xy_rmse', NaN, ...
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
        name = 'Dual INDI';
    elseif acc_fraction < 0.1 && rate_fraction > 0.9
        name = 'Rate INDI';
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
    m.along_rmse = sqrt(mean(c.along_error.^2, 'omitnan'));
    m.cross_rmse = sqrt(mean(c.cross_error.^2, 'omitnan'));
end

function [best_delay, best_rmse] = estimate_xy_delay(t_ref, p_ref, t, p, limit, step)
    delays = (-limit:step:limit).';
    costs = NaN(size(delays));

    for k = 1:numel(delays)
        % Positive delay means the response matches an earlier reference.
        delayed_reference = interp_series(t_ref, p_ref, t - delays(k));
        error = p - delayed_reference;
        costs(k) = sqrt(mean(sum(error.^2, 2), 'omitnan'));
    end

    [best_rmse, idx] = min(costs);
    best_delay = delays(idx);
end


function y = constrain(x, lower, upper)
    y = min(max(x, lower), upper);
end
