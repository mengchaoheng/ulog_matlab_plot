clear; close all; clc;

%% User settings
root_dir = fileparts(mfilename('fullpath'));
mat_file = fullfile(root_dir, 'data', ...
    'log_37_2026-8-10-19-19-22.mat');
output_dir = fullfile(root_dir, 'results');
output_trajectory_pdf = fullfile(output_dir, 'figure8_tracking_realflight.pdf');
output_velocity_pdf = fullfile(output_dir, 'figure8_velocity_tracking_realflight.pdf');
output_attitude_pdf = fullfile(output_dir, 'figure8_attitude_tracking_realflight.pdf');
output_rate_pdf = fullfile(output_dir, 'figure8_rate_tracking_realflight.pdf');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Translational flight window: two complete figure-eight periods.
t_start  = 57.60;
period   = 31.40;
n_cycles = 2;
t_end    = t_start + n_cycles * period;

% Rotational-state window: one complete period after the switch transient.
state_t_start  = 60.00;
state_n_cycles = 1;
state_t_end    = state_t_start + state_n_cycles * period;

% Plot-only downsampling, following /Users/mch/Proj/PINDI/test5_2_1_1.m.
% All RMSE values below are calculated from the full-rate data.
attitude_response_step = 10;  % about 245 Hz -> 24.5 Hz
attitude_setpoint_step = 5;   % about 122 Hz -> 24.4 Hz
rate_response_step     = 20;  % about 242 Hz -> 12.1 Hz
rate_setpoint_step     = 10;  % about 245 Hz -> 24.5 Hz

%% Load data
S = load(mat_file, 'log');
D = S.log.data;

pos_sp_tbl = D.vehicle_local_position_setpoint_0;
pos_pv_tbl = D.vehicle_local_position_0;
t_pos_sp = double(pos_sp_tbl.timestamp) * 1e-6;
t_pos_pv = double(pos_pv_tbl.timestamp) * 1e-6;
pos_sp = double([pos_sp_tbl.x, pos_sp_tbl.y, -pos_sp_tbl.z]);
pos_pv = double([pos_pv_tbl.x, pos_pv_tbl.y, -pos_pv_tbl.z]);
vel_sp = double([pos_sp_tbl.vx, pos_sp_tbl.vy, -pos_sp_tbl.vz]);
vel_pv = double([pos_pv_tbl.vx, pos_pv_tbl.vy, -pos_pv_tbl.vz]);

att_sp_tbl = D.vehicle_attitude_setpoint_0;
att_pv_tbl = D.vehicle_attitude_0;
t_att_sp = double(att_sp_tbl.timestamp) * 1e-6;
t_att_pv = double(att_pv_tbl.timestamp) * 1e-6;
q_sp = double([att_sp_tbl.q_d_0_, att_sp_tbl.q_d_1_, ...
    att_sp_tbl.q_d_2_, att_sp_tbl.q_d_3_]);
q_pv = double([att_pv_tbl.q_0_, att_pv_tbl.q_1_, ...
    att_pv_tbl.q_2_, att_pv_tbl.q_3_]);
eul_sp = quat2eul(q_sp, 'ZYX');
eul_pv = quat2eul(q_pv, 'ZYX');
att_sp = [eul_sp(:, 3), eul_sp(:, 2), unwrap(eul_sp(:, 1))];
att_pv = [eul_pv(:, 3), eul_pv(:, 2), unwrap(eul_pv(:, 1))];

rate_sp_tbl = D.vehicle_rates_setpoint_0;
rate_pv_tbl = D.vehicle_angular_velocity_0;
t_rate_sp = double(rate_sp_tbl.timestamp) * 1e-6;
t_rate_pv = double(rate_pv_tbl.timestamp) * 1e-6;
rate_sp = double([rate_sp_tbl.roll, rate_sp_tbl.pitch, rate_sp_tbl.yaw]);
rate_pv = double([rate_pv_tbl.xyz_0_, rate_pv_tbl.xyz_1_, ...
    rate_pv_tbl.xyz_2_]);

%% Window masks and full-rate tracking metrics
use_pos_sp = t_pos_sp >= t_start & t_pos_sp <= t_end & ...
    all(isfinite(pos_sp), 2);
use_pos_pv = t_pos_pv >= t_start & t_pos_pv <= t_end & ...
    all(isfinite(pos_pv), 2);
use_vel_sp = t_pos_sp >= t_start & t_pos_sp <= t_end & ...
    all(isfinite(vel_sp), 2);
use_vel_pv = t_pos_pv >= t_start & t_pos_pv <= t_end & ...
    all(isfinite(vel_pv), 2);
use_att_sp = t_att_sp >= state_t_start & t_att_sp <= state_t_end & ...
    all(isfinite(att_sp), 2);
use_att_pv = t_att_pv >= state_t_start & t_att_pv <= state_t_end & ...
    all(isfinite(att_pv), 2);
use_rate_sp = t_rate_sp >= state_t_start & t_rate_sp <= state_t_end & ...
    all(isfinite(rate_sp), 2);
use_rate_pv = t_rate_pv >= state_t_start & t_rate_pv <= state_t_end & ...
    all(isfinite(rate_pv), 2);

assert(nnz(use_pos_sp) > 2 && nnz(use_pos_pv) > 2, ...
    'No valid position data in the selected window.');

pos_rmse = tracking_rmse(t_pos_sp(use_pos_sp), pos_sp(use_pos_sp, :), ...
    t_pos_pv(use_pos_pv), pos_pv(use_pos_pv, :), false);
vel_rmse = tracking_rmse(t_pos_sp(use_vel_sp), vel_sp(use_vel_sp, :), ...
    t_pos_pv(use_vel_pv), vel_pv(use_vel_pv, :), false);
att_rmse = tracking_rmse(t_att_sp(use_att_sp), att_sp(use_att_sp, :), ...
    t_att_pv(use_att_pv), att_pv(use_att_pv, :), true) * 180 / pi;
rate_rmse = tracking_rmse(t_rate_sp(use_rate_sp), rate_sp(use_rate_sp, :), ...
    t_rate_pv(use_rate_pv), rate_pv(use_rate_pv, :), false) * 180 / pi;

%% Paper style
set(groot, ...
    'defaultAxesFontSize', 6, ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultAxesLineWidth', 0.35, ...
    'defaultAxesLabelFontSizeMultiplier', 1, ...
    'defaultTextFontName', 'Times New Roman', ...
    'defaultLegendFontName', 'Times New Roman', ...
    'defaultTextInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex');

% RGB color, line style, and line width (pt).
reference_color = [0.850, 0.325, 0.098];
response_color = [0.000, 0.447, 0.741];
reference_style = {'Color', reference_color, 'LineStyle', '--', 'LineWidth', 0.8};
response_style = {'Color', response_color, 'LineStyle', '-', 'LineWidth', 0.5};

%% Figure 1: 3-D position tracking
fig_traj = figure('Name', '3-D trajectory tracking', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [0, 0, 7, 5]);
ax_traj = axes(fig_traj); hold(ax_traj, 'on');
% A dense 3-D line can make MATLAB's native '--' renderer appear solid.
% Draw the reference as explicit arc-length dash segments, as in
% plot_paper_figures.m, and use a NaN proxy for its legend entry.
plot3DashedHidden(ax_traj, pos_sp(use_pos_sp, :).', 1.00, 0.50, ...
    'Color', reference_color, 'LineWidth', 0.8);
h_traj_sp = plot3(ax_traj, NaN, NaN, NaN, reference_style{:});
h_traj_res = plot3(ax_traj, pos_pv(use_pos_pv, 1), ...
    pos_pv(use_pos_pv, 2), pos_pv(use_pos_pv, 3), response_style{:});
xlabel(ax_traj, '$p_x$ (m)');
ylabel(ax_traj, '$p_y$ (m)');
zlabel(ax_traj, '$-p_z$ (m)');
grid(ax_traj, 'on'); box(ax_traj, 'on');
view(ax_traj, 42, 25);
zlim(ax_traj, [3, 6]);
zticks(ax_traj, 0:2:6);
daspect(ax_traj, [1, 1, 1]);
lgd = legend(ax_traj, [h_traj_sp, h_traj_res], ...
    {'Reference', 'Response'}, 'Location', 'best', 'NumColumns', 2);
lgd.ItemTokenSize = [20, 8];
lgd.Position = [0.288192769567176 0.757763726455958 ...
    0.450828360552763 0.0675616197183098];
export_paper_figure(fig_traj, output_trajectory_pdf, 7, 5);

%% Figure 2: velocity tracking
fig_vel = figure('Name', 'Velocity tracking', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [0, 0, 7, 6]);
tl_vel = tiledlayout(fig_vel, 3, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
vel_labels = {'$v_x$ (m/s)', '$v_y$ (m/s)', '$-v_z$ (m/s)'};
ax_vel = gobjects(1, 3);
for k = 1:3
    ax_vel(k) = nexttile(tl_vel); hold(ax_vel(k), 'on');
    h_vel_res = plot(ax_vel(k), t_pos_pv(use_vel_pv) - t_start, ...
        vel_pv(use_vel_pv, k), response_style{:});
    h_vel_sp = plot(ax_vel(k), t_pos_sp(use_vel_sp) - t_start, ...
        vel_sp(use_vel_sp, k), reference_style{:});
    ylabel(ax_vel(k), vel_labels{k});
    xlim(ax_vel(k), [0, n_cycles * period]);
    grid(ax_vel(k), 'on'); box(ax_vel(k), 'on');
end
ylim(ax_vel(1), [-2.5, 2.5]);
ylim(ax_vel(2), [-2.5, 2.5]);
ylim(ax_vel(3), [-0.2, 0.2]);
lgd_vel = legend(ax_vel(3), [h_vel_sp, h_vel_res], ...
    {'Reference', 'Response'}, 'Location', 'best', 'NumColumns', 2);
lgd_vel.ItemTokenSize = [14, 8];
xlabel(ax_vel(3), 'Time (s)');
linkaxes(ax_vel, 'x');
export_paper_figure(fig_vel, output_velocity_pdf, 7, 6);

%% Figure 3: attitude tracking
idx_att_sp = decimated_indices(use_att_sp, attitude_setpoint_step);
idx_att_pv = decimated_indices(use_att_pv, attitude_response_step);
idx_rate_sp = decimated_indices(use_rate_sp, rate_setpoint_step);
idx_rate_pv = decimated_indices(use_rate_pv, rate_response_step);

fig_att = figure('Name', 'Attitude tracking', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [0, 0, 7, 6]);
tl_att = tiledlayout(fig_att, 3, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
att_labels = {'$\varphi$ (deg)', '$\theta$ (deg)', '$\psi$ (deg)'};
ax_att = gobjects(1, 3);

for k = 1:3
    ax_att(k) = nexttile(tl_att); hold(ax_att(k), 'on');
    h_att_res = plot(ax_att(k), t_att_pv(idx_att_pv) - state_t_start, ...
        att_pv(idx_att_pv, k) * 180 / pi, response_style{:});
    h_att_sp = plot(ax_att(k), t_att_sp(idx_att_sp) - state_t_start, ...
        att_sp(idx_att_sp, k) * 180 / pi, reference_style{:});
    ylabel(ax_att(k), att_labels{k});
    xlim(ax_att(k), [0, state_n_cycles * period]);
    grid(ax_att(k), 'on'); box(ax_att(k), 'on');
end
ylim(ax_att(1), [-5, 5]);
ylim(ax_att(2), [-10, 10]);
ylim(ax_att(3), [24, 27]);
lgd_att = legend(ax_att(3), [h_att_sp, h_att_res], ...
    {'Reference', 'Response'}, 'Location', 'best', 'NumColumns', 2);
lgd_att.ItemTokenSize = [14, 8];
xlabel(ax_att(3), 'Time (s)');
linkaxes(ax_att, 'x');
export_paper_figure(fig_att, output_attitude_pdf, 7, 6);

%% Figure 4: angular-rate tracking
fig_rate = figure('Name', 'Angular-rate tracking', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [0, 0, 7, 6]);
tl_rate = tiledlayout(fig_rate, 3, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
rate_labels = {'$\omega_x^b$ (deg/s)', '$\omega_y^b$ (deg/s)', ...
    '$\omega_z^b$ (deg/s)'};
ax_rate = gobjects(1, 3);

for k = 1:3
    ax_rate(k) = nexttile(tl_rate); hold(ax_rate(k), 'on');
    h_rate_res = plot(ax_rate(k), t_rate_pv(idx_rate_pv) - state_t_start, ...
        rate_pv(idx_rate_pv, k) * 180 / pi, response_style{:});
    h_rate_sp = plot(ax_rate(k), t_rate_sp(idx_rate_sp) - state_t_start, ...
        rate_sp(idx_rate_sp, k) * 180 / pi, reference_style{:});
    ylabel(ax_rate(k), rate_labels{k});
    xlim(ax_rate(k), [0, state_n_cycles * period]);
    grid(ax_rate(k), 'on'); box(ax_rate(k), 'on');
end
ylim(ax_rate(1), [-10, 10]);
ylim(ax_rate(2), [-15, 15]);
ylim(ax_rate(3), [-5, 5]);
lgd_rate = legend(ax_rate(3), [h_rate_sp, h_rate_res], ...
    {'Reference', 'Response'}, 'Location', 'best', 'NumColumns', 2);
lgd_rate.ItemTokenSize = [14, 8];
xlabel(ax_rate(3), 'Time (s)');
linkaxes(ax_rate, 'x');
export_paper_figure(fig_rate, output_rate_pdf, 7, 6);

%% Console summary
fprintf('Trajectory/velocity window: %.2f--%.2f s (%d cycles)\n', ...
    t_start, t_end, n_cycles);
fprintf('Attitude/rate window: %.2f--%.2f s (%d cycle)\n', ...
    state_t_start, state_t_end, state_n_cycles);
fprintf('Position RMSE [px py -pz] = [%.3f %.3f %.3f] m\n', pos_rmse);
fprintf('Velocity RMSE [vx vy -vz] = [%.3f %.3f %.3f] m/s\n', vel_rmse);
fprintf('Attitude RMSE [roll pitch yaw] = [%.2f %.2f %.2f] deg\n', ...
    att_rmse);
fprintf(['Rate RMSE [omega_x^b omega_y^b omega_z^b] = ' ...
    '[%.2f %.2f %.2f] deg/s\n'], rate_rmse);
fprintf('Saved:\n  %s\n  %s\n  %s\n  %s\n', output_trajectory_pdf, ...
    output_velocity_pdf, output_attitude_pdf, output_rate_pdf);

%% Local helpers
function idx = decimated_indices(mask, step)
    idx = find(mask);
    idx = idx(1:step:end);
end

function rmse = tracking_rmse(t_sp, y_sp, t_pv, y_pv, wrap_yaw)
    [t_sp, ia] = unique(t_sp, 'stable');
    y_sp = y_sp(ia, :);
    y_sp_i = interp1(t_sp, y_sp, t_pv, 'linear', NaN);
    err = y_pv - y_sp_i;
    if wrap_yaw
        err(:, 3) = atan2(sin(err(:, 3)), cos(err(:, 3)));
    end
    rmse = sqrt(mean(err.^2, 1, 'omitnan'));
end

function h = plot3DashedHidden(ax, p, dashPeriod, duty, varargin)
%PLOT3DASHEDHIDDEN Draw a 3-D polyline with reliable arc-length dashes.
% Native dashed rendering can disappear on dense vector trajectories, so
% each visible dash is emitted as its own short solid line segment.

    p = double(p);
    valid = all(isfinite(p), 1);
    p = p(:, valid);

    if size(p, 2) < 2
        h = gobjects(0);
        return;
    end

    ds = sqrt(sum(diff(p, 1, 2).^2, 1));
    s = [0, cumsum(ds)];
    [s, ia] = unique(s, 'stable');
    p = p(:, ia);

    dashLength = dashPeriod * duty;
    totalLength = s(end);
    holdState = ishold(ax);
    hold(ax, 'on');
    h = gobjects(0);

    for s0 = 0:dashPeriod:totalLength
        s1 = min(s0 + dashLength, totalLength);
        q0 = interpPointByArcLength(p, s, s0);
        q1 = interpPointByArcLength(p, s, s1);
        inside = s > s0 & s < s1;
        q = [q0, p(:, inside), q1];
        h(end+1) = plot3(ax, q(1,:), q(2,:), q(3,:), ...
            'LineStyle', '-', 'HandleVisibility', 'off', varargin{:});
    end

    if ~holdState
        hold(ax, 'off');
    end
end

function q = interpPointByArcLength(p, s, sq)
    q = zeros(3, 1);
    for k = 1:3
        q(k) = interp1(s, p(k,:), sq, 'linear');
    end
end

function export_paper_figure(fig, filename, width_cm, height_cm)
    set(fig, 'Units', 'centimeters', 'Position', [0, 0, width_cm, height_cm]);
    exportgraphics(fig, filename, 'ContentType', 'vector', ...
        'BackgroundColor', 'white', 'Resolution', 1500, ...
        'Width', width_cm, 'Height', height_cm, ...
        'Padding', 'tight', 'Units', 'centimeters');
end
