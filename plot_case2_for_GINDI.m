clear all;
close all;
clc;

run('load_data_main.m');   %12_09_34.ulg
addpath(genpath(pwd));

%% =========================================================================
%  Global figure / font style
% =========================================================================
set(groot, ...
    'defaultAxesFontSize', 6, ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultAxesLineWidth', 0.5, ...
    'defaultAxesLabelFontSizeMultiplier', 1, ...
    'defaultAxesTitleFontSizeMultiplier', 1);

set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');

%% =========================================================================
%  Unified named-color manager
% =========================================================================
sty = plot_style_manager();

% ---- Local line-width definition ----
LW_XS = 0.5;
LW_S  = 0.5;
LW_M  = 0.8;
LW_L  = 1.1;
LW_XL = 1.4;

LW_SP = LW_M;
LW_RES = LW_S;
LW_MAIN = LW_M;
LW_MAIN_BOLD = LW_L;
LW_MULTI = LW_M;
LW_MULTI_BOLD = LW_XL;

STYLE_SP  = plot_style(sty, 'setpoint1', '--', LW_SP);
STYLE_RES = plot_style(sty, 'response1', '-',  LW_RES);
STYLE_AXIS1 = plot_style(sty, 'setpoint2', '-',  LW_MAIN);
STYLE_AXIS2 = plot_style(sty, 'response2', '--', LW_MAIN);
STYLE_AXIS3 = plot_style(sty, 'setpoint3', '-.', LW_MAIN);
STYLE_AXIS4 = plot_style(sty, 'response3', ':',  LW_MAIN);
STYLE_AXIS5 = plot_style(sty, 'setpoint4', '-',  LW_MAIN);
STYLE_AXIS1_BOLD = plot_style(sty, 'setpoint2', '-',  LW_MAIN_BOLD);
STYLE_AXIS2_BOLD = plot_style(sty, 'response2', '--', LW_MAIN_BOLD);
STYLE_AXIS3_BOLD = plot_style(sty, 'setpoint3', '-.', LW_MAIN_BOLD);
STYLE_AXIS4_BOLD = plot_style(sty, 'response3', ':',  LW_MAIN_BOLD);
STYLE_AXIS5_BOLD = plot_style(sty, 'setpoint4', '-',  LW_MAIN_BOLD);
STYLE_AXIS1_THIN = plot_style(sty, 'setpoint2', '-',  LW_XS);
STYLE_AXIS2_THIN = plot_style(sty, 'response2', '--', LW_XS);
STYLE_AXIS3_THIN = plot_style(sty, 'setpoint3', '-.', LW_XS);
STYLE_AXIS1_DASH = plot_style(sty, 'setpoint2', '--', LW_XS);
STYLE_AXIS2_DASH = plot_style(sty, 'response2', '-.', LW_XS);
STYLE_AXIS4_DOT  = plot_style(sty, 'response3', ':',  LW_MAIN);
STYLE_THRUST     = plot_style(sty, 'neutral1', '-',  LW_MAIN_BOLD);
STYLE_THRUST_ALT = plot_style(sty, 'neutral1', '--', LW_MAIN);
STYLE_THRUST_BOLD = plot_style(sty, 'neutral1', '-', LW_XL);

C_AXIS1 = sty.color.setpoint2;
C_AXIS2 = sty.color.response2;
C_AXIS3 = sty.color.setpoint3;
C_AXIS4 = sty.color.response3;
C_AXIS5 = sty.color.setpoint4;
C_AXIS6 = sty.color.response4;
C_AXIS7 = sty.color.setpoint5;
COLOR_WARN_FILL = sty.color.warn_fill;
 
%% =========================================================================
%  Three-axis tracking error trajectories in three algorithm windows
% =========================================================================
if(exist('XYZ', 'var') && exist('XYZ_setpoint', 'var') && ...
        exist('vehicle_local_position_t', 'var') && ...
        exist('vehicle_local_position_setpoint_t', 'var'))

    fig = figure('Name', 'Three-axis Tracking Error', ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [0, 0, 7, 7], ...
    'PaperUnits', 'centimeters', ...
    'PaperSize', [7, 7], ...
    'PaperPosition', [0, 0, 7, 7], ...
    'PaperPositionMode', 'auto', ...
    'InvertHardcopy', 'off');

    % Start times are arranged as [PX4, PINDI, GINDI]
    start_time=30;
    segment_start_times = [start_time, start_time+35, start_time+35*2];    % s
    segment_duration    = 20;                % s
    segment_names       = {'Geometric Control', 'PINDI', 'GINDI'};
    STYLE_1 = {'Color', [0.494, 0.184, 0.556], ...
                  'LineStyle', '--', ...
                  'LineWidth', 0.5};
    STYLE_2 = {'Color', [0.96, 0.62, 0.26], ...
                  'LineStyle', '-.', ...
                  'LineWidth', 0.5};
    STYLE_3 = {'Color', [0.000, 0.447, 0.741], ...
                  'LineStyle', '-', ...
                  'LineWidth', 0.7};

    segment_styles      = {STYLE_1, STYLE_2,  STYLE_3};

    t_res = vehicle_local_position_t(:);
    t_sp  = vehicle_local_position_setpoint_t(:);

    [t_sp_unique, ia] = unique(t_sp, 'stable');
    XYZ_sp_unique = XYZ_setpoint(ia, :);

    XYZ_sp_interp = interp1(t_sp_unique, XYZ_sp_unique, t_res, 'linear', 'extrap');

    % Tracking error: response minus setpoint
    e_xyz = XYZ - XYZ_sp_interp;

    ylabels = {'$e_x$ (m)', '$e_y$ (m)', '$e_z$ (m)'};
    ax = [];

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
            % lgd.Box = 'off';
        
            drawnow;
        
            lgd.Units = 'normalized';
            pos = lgd.Position;
            pos(1) = pos(1) + 0.03;    % 往左，按需调
            pos(2) = pos(2) + 0.025;   % 往上，保留你觉得合适的高度
            lgd.Position = pos;
        end

        % if i == 3
            xlabel('Time (s)');
        % end
    end

    linkaxes(ax, 'x');
    drawnow;
    PlotToFile(fig, 'results/tracking_error_xyz.pdf', 7, 7);
end