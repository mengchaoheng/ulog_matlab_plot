%% plot_case2_position_rmse_boxplot_color_only.m
% Position tracking RMSE boxplot for case2 repeated experiments.
% Each MAT file corresponds to one algorithm and contains 50 repeated trials.
%
% Trial timing rule:
%   first trial starts at 50 s;
%   each trial lasts 20 s;
%   each pause lasts 8 s;
%   the rule repeats 50 times.

clc; clear; close all;
addpath(genpath(pwd));

%% ========================================================================
%  User settings
% =========================================================================
base_dir = './data/case2';

algorithm_files = { ...
    '11_48_19.mat', ...   % PX4
    '11_52_03.mat', ...   % PINDI
    '11_59_33.mat'};      % GINDI

custom_labels = {'PX4', 'PINDI', 'GINDI'};

experiment_start_time = 50;   % [s]
experiment_duration   = 20;   % [s]
pause_duration        = 8;    % [s]
num_experiments       = 50;
min_samples_per_trial = 5;

save_plots = true;
result_dir = './results/case2_position_rmse';

%% ========================================================================
%  Figure style, kept consistent with test5_2_2RMS.m
% =========================================================================
set(groot, ...
    'defaultAxesFontSize', 8, ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultAxesLineWidth', 0.5, ...
    'defaultAxesLabelFontSizeMultiplier', 1, ...
    'defaultAxesTitleFontSizeMultiplier', 1);

alg_colors = [ ...
    0.15 0.40 0.75;   % PX4
    0.70 0.30 0.05;   % PINDI
    0.20 0.55 0.25];  % GINDI

%% ========================================================================
%  Compute RMSE samples
% =========================================================================
trial_start_times = experiment_start_time + ...
    (0:num_experiments-1) * (experiment_duration + pause_duration);
trial_end_times = trial_start_times + experiment_duration;

all_rmse_position = [];
all_rmse_xyz = [];
group_ids = [];
rmse_table_all = table();

for alg_id = 1:numel(algorithm_files)
    file_path = fullfile(base_dir, algorithm_files{alg_id});

    if ~isfile(file_path)
        error('MAT file not found: %s', file_path);
    end

    data = load(file_path);
    [t_pos, xyz, t_sp, xyz_sp] = extract_position_data(data, file_path);

    [rmse_position_each, rmse_xyz_each, n_samples_each] = compute_trial_position_rmse( ...
        t_pos, xyz, t_sp, xyz_sp, trial_start_times, trial_end_times, min_samples_per_trial);

    valid_idx = isfinite(rmse_position_each);
    rmse_position_valid = rmse_position_each(valid_idx);
    rmse_xyz_valid = rmse_xyz_each(valid_idx, :);
    valid_trial_id = find(valid_idx);

    all_rmse_position = [all_rmse_position; rmse_position_valid(:)];
    all_rmse_xyz = [all_rmse_xyz; rmse_xyz_valid];
    group_ids = [group_ids; repmat(alg_id, numel(rmse_position_valid), 1)];

    rmse_table_alg = table( ...
        repmat(string(custom_labels{alg_id}), numel(valid_trial_id), 1), ...
        valid_trial_id(:), ...
        trial_start_times(valid_trial_id).', ...
        trial_end_times(valid_trial_id).', ...
        n_samples_each(valid_idx), ...
        rmse_xyz_valid(:,1), ...
        rmse_xyz_valid(:,2), ...
        rmse_xyz_valid(:,3), ...
        rmse_position_valid(:), ...
        'VariableNames', {'Algorithm', 'Trial', 'StartTime_s', 'EndTime_s', ...
        'NumSamples', 'RMSE_x_m', 'RMSE_y_m', 'RMSE_z_m', 'RMSE_position_m'});

    if isempty(rmse_table_all)
        rmse_table_all = rmse_table_alg;
    else
        rmse_table_all = [rmse_table_all; rmse_table_alg];
    end

    fprintf('%s: %d valid trials out of %d.\n', ...
        custom_labels{alg_id}, numel(rmse_position_valid), num_experiments);
end

if isempty(all_rmse_position)
    error('No valid position RMSE samples were computed. Check timestamps and experiment timing settings.');
end

%% ========================================================================
%  Boxplot
% =========================================================================
fig1 = figure(1);
set(fig1, 'Color', 'w');
hold on;

h_box = gobjects(numel(custom_labels), 1);
for alg_id = 1:numel(custom_labels)
    idx = (group_ids == alg_id);
    h_box(alg_id) = boxchart(ones(sum(idx),1) * alg_id, all_rmse_position(idx), ...
        'BoxFaceColor', alg_colors(alg_id, :), ...
        'BoxEdgeColor', 'k', ...
        'BoxWidth', 0.35, ...
        'LineWidth', 0.8, ...
        'MarkerStyle', 'x', ...
        'MarkerColor', alg_colors(alg_id, :), ...
        'MarkerSize', 4);
end

set(gca, ...
    'XTick', 1:numel(custom_labels), ...
    'XTickLabel', custom_labels);
xlim([0.5, numel(custom_labels)+0.5]);
grid on;
ylabel('Position RMSE (m)');

% Optional title can be enabled if needed.
% title('Position tracking RMSE across repeated trials');

%% ========================================================================
%  Save figure and data
% =========================================================================
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

save(fullfile(result_dir, 'position_rmse_boxplot_data.mat'), ...
    'all_rmse_position', 'all_rmse_xyz', 'group_ids', 'custom_labels', ...
    'algorithm_files', 'trial_start_times', 'trial_end_times', ...
    'experiment_start_time', 'experiment_duration', 'pause_duration', ...
    'num_experiments', 'rmse_table_all');

writetable(rmse_table_all, fullfile(result_dir, 'position_rmse_samples.csv'));

if save_plots
    pdf_path = fullfile(result_dir, 'position_rmse_boxplot_case2.pdf');
    png_path = fullfile(result_dir, 'position_rmse_boxplot_case2.png');

    if exist('PlotToFileColorPDF', 'file') == 2
        PlotToFileColorPDF(fig1, pdf_path, 8, 5);
    else
        set(fig1, 'Units', 'centimeters', 'Position', [2 2 8 5]);
        exportgraphics(fig1, pdf_path, 'ContentType', 'vector');
    end

    exportgraphics(fig1, png_path, 'Resolution', 300);
end

fprintf('\nPosition RMSE samples saved to:\n');
disp(fullfile(result_dir, 'position_rmse_boxplot_data.mat'));
disp(fullfile(result_dir, 'position_rmse_samples.csv'));

%% ========================================================================
%  Local functions
% =========================================================================
function [t_pos, xyz, t_sp, xyz_sp] = extract_position_data(data, file_path)
    if isfield(data, 'log') && isfield(data.log, 'data')
        log_data = data.log.data;

        if ~isfield(log_data, 'vehicle_local_position_0')
            error('vehicle_local_position_0 not found in %s.', file_path);
        end
        if ~isfield(log_data, 'vehicle_local_position_setpoint_0')
            error('vehicle_local_position_setpoint_0 not found in %s.', file_path);
        end

        pos_tbl = log_data.vehicle_local_position_0;
        sp_tbl  = log_data.vehicle_local_position_setpoint_0;

        t_pos = normalize_time(pos_tbl.timestamp);
        xyz = [pos_tbl.x, pos_tbl.y, pos_tbl.z];

        t_sp = normalize_time(sp_tbl.timestamp);
        xyz_sp = [sp_tbl.x, sp_tbl.y, sp_tbl.z];

    elseif all(isfield(data, {'vehicle_local_position_t', 'XYZ', ...
            'vehicle_local_position_setpoint_t', 'XYZ_setpoint'}))

        t_pos = normalize_time(data.vehicle_local_position_t);
        xyz = data.XYZ;
        t_sp = normalize_time(data.vehicle_local_position_setpoint_t);
        xyz_sp = data.XYZ_setpoint;

    else
        error(['Position variables not found in %s. Expected either raw log.data ' ...
            'or variables vehicle_local_position_t, XYZ, vehicle_local_position_setpoint_t, XYZ_setpoint.'], file_path);
    end

    [t_pos, pos_order] = sort(t_pos(:));
    xyz = double(xyz(pos_order, :));

    [t_sp, sp_order] = sort(t_sp(:));
    xyz_sp = double(xyz_sp(sp_order, :));

    valid_pos = isfinite(t_pos) & all(isfinite(xyz), 2);
    valid_sp = isfinite(t_sp) & all(isfinite(xyz_sp), 2);

    t_pos = t_pos(valid_pos);
    xyz = xyz(valid_pos, :);
    t_sp = t_sp(valid_sp);
    xyz_sp = xyz_sp(valid_sp, :);

    [t_sp, ia] = unique(t_sp, 'stable');
    xyz_sp = xyz_sp(ia, :);
end

function t = normalize_time(t_raw)
    t = double(t_raw(:));

    % PX4 ULog timestamps are usually in microseconds. Data extracted by
    % load_data_main.m may already be in seconds.
    finite_idx = isfinite(t);
    if any(finite_idx) && max(abs(t(finite_idx))) > 1e5
        t = t * 1e-6;
    end
end

function [rmse_position, rmse_xyz, n_samples] = compute_trial_position_rmse( ...
        t_pos, xyz, t_sp, xyz_sp, trial_start_times, trial_end_times, min_samples_per_trial)

    num_trials = numel(trial_start_times);
    rmse_position = nan(num_trials, 1);
    rmse_xyz = nan(num_trials, 3);
    n_samples = zeros(num_trials, 1);

    t_sp_min = min(t_sp);
    t_sp_max = max(t_sp);

    for trial_id = 1:num_trials
        t0 = trial_start_times(trial_id);
        t1 = trial_end_times(trial_id);

        valid_t0 = max(t0, t_sp_min);
        valid_t1 = min(t1, t_sp_max);

        idx = (t_pos >= valid_t0) & (t_pos < valid_t1);

        if nnz(idx) < min_samples_per_trial
            warning('Trial %d has %d samples in [%.3f, %.3f] s and was skipped.', ...
                trial_id, nnz(idx), t0, t1);
            continue;
        end

        xyz_sp_interp = interp1(t_sp, xyz_sp, t_pos(idx), 'linear');
        err_xyz = xyz(idx, :) - xyz_sp_interp;
        valid_err = all(isfinite(err_xyz), 2);

        if nnz(valid_err) < min_samples_per_trial
            warning('Trial %d has %d valid interpolated samples and was skipped.', ...
                trial_id, nnz(valid_err));
            continue;
        end

        err_xyz = err_xyz(valid_err, :);
        n_samples(trial_id) = size(err_xyz, 1);

        rmse_xyz(trial_id, :) = sqrt(mean(err_xyz.^2, 1));
        rmse_position(trial_id) = sqrt(mean(sum(err_xyz.^2, 2)));
    end
end
