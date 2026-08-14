clear;
close all;
clc;

%% Load data
root_dir = fileparts(mfilename('fullpath'));
mat_file = fullfile(root_dir, 'data', '12_09_34.mat');
output_file = fullfile(root_dir, 'results', 'tracking_error_xyz.pdf');

assert(isfile(mat_file), 'MAT file not found: %s', mat_file);
S = load(mat_file, 'log');
assert(isfield(S, 'log') && isfield(S.log, 'data'), ...
    'Expected variable log.data is missing from %s.', mat_file);
D = S.log.data;

required_topics = {'vehicle_local_position_0', ...
    'vehicle_local_position_setpoint_0'};
for k = 1:numel(required_topics)
    assert(isfield(D, required_topics{k}), ...
        'Required topic is missing: %s', required_topics{k});
end

pos_tbl = D.vehicle_local_position_0;
pos_sp_tbl = D.vehicle_local_position_setpoint_0;
vehicle_local_position_t = double(pos_tbl.timestamp) * 1e-6;
vehicle_local_position_setpoint_t = double(pos_sp_tbl.timestamp) * 1e-6;
XYZ = double([pos_tbl.x, pos_tbl.y, pos_tbl.z]);
XYZ_setpoint = double([pos_sp_tbl.x, pos_sp_tbl.y, pos_sp_tbl.z]);

%% Paper style
set(groot, ...
    'defaultAxesFontSize', 6, ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultAxesLineWidth', 0.5, ...
    'defaultAxesLabelFontSizeMultiplier', 1, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex');
 
%% Three-axis tracking errors in the three controller windows
fig = figure('Name', 'Three-axis Tracking Error', ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [0, 0, 7, 7], ...
    'PaperUnits', 'centimeters', ...
    'PaperSize', [7, 7], ...
    'PaperPosition', [0, 0, 7, 7], ...
    'PaperPositionMode', 'auto', ...
    'InvertHardcopy', 'off');

% Order: Baseline, PINDI, GINDI. Each style explicitly specifies RGB,
% line style, and line width in points.
segment_start_times = [30, 65, 100];
segment_duration = 20;
segment_names = {'Baseline', 'PINDI', 'GINDI'};
segment_styles = { ...
    {'Color', [0.494, 0.184, 0.556], 'LineStyle', '--', 'LineWidth', 0.5}, ...
    {'Color', [0.960, 0.620, 0.260], 'LineStyle', '-.', 'LineWidth', 0.5}, ...
    {'Color', [0.000, 0.447, 0.741], 'LineStyle', '-',  'LineWidth', 0.7}};

t_res = vehicle_local_position_t(:);
t_sp  = vehicle_local_position_setpoint_t(:);

[t_sp_unique, ia] = unique(t_sp, 'stable');
XYZ_sp_unique = XYZ_setpoint(ia, :);

XYZ_sp_interp = interp1(t_sp_unique, XYZ_sp_unique, t_res, 'linear', 'extrap');

% Tracking error: response minus setpoint
e_xyz = XYZ - XYZ_sp_interp;

ylabels = {'$e_x$ (m)', '$e_y$ (m)', '$e_z$ (m)'};
ax = gobjects(1, 3);

for i = 1:3
    ax(i) = subplot(3, 1, i); hold on;

    for k = 1:numel(segment_start_times)
        t0 = segment_start_times(k);
        t1 = t0 + segment_duration;
        idx = (t_res >= t0) & (t_res <= t1);

        if nnz(idx) >= 2
            t_seg = t_res(idx) - t0;
            s = segment_styles{k};
            plot(t_seg, e_xyz(idx, i), s{:}, 'DisplayName', segment_names{k});
        end
    end

    grid on;
    ylabel(ylabels{i});
    xlim([0, segment_duration]);

    if i == 1
        lgd = legend('Location', 'best', 'NumColumns', 3);
        lgd.ItemTokenSize = [20, 18];
        drawnow;
        lgd.Units = 'normalized';
        pos = lgd.Position;
        pos(1) = pos(1) + 0.03;
        pos(2) = pos(2) + 0.025;
        lgd.Position = pos;
    end

    xlabel('Time (s)');
end

linkaxes(ax, 'x');
drawnow;
if ~exist(fileparts(output_file), 'dir')
    mkdir(fileparts(output_file));
end
export_paper_figure(fig, output_file, 7, 7);

function export_paper_figure(fig, filename, width_cm, height_cm)
    set(fig, 'Units', 'centimeters', 'Position', [0, 0, width_cm, height_cm]);
    exportgraphics(fig, filename, 'ContentType', 'vector', ...
        'BackgroundColor', 'white', 'Resolution', 1500, ...
        'Width', width_cm, 'Height', height_cm, ...
        'Padding', 'tight', 'Units', 'centimeters');
end
