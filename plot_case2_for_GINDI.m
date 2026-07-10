clear all;
close all;
clc;

run('load_data_main.m'); % load('flight_data.mat'); 20_00_18/06_44_06/09_29_53////10_44_21/13_31_44
addpath(genpath(pwd));

%% =========================================================================
%  Global figure / font style
% =========================================================================
set(groot, ...
    'defaultAxesFontSize', 5, ...
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
LW_M  = 0.7;
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
%  Figure 1 2, 3, 4, 5: (vehicle_angular_velocity, Attitude, Vel, Pos, Control)
% =========================================================================
% --- Figure 1: vehicle_angular_velocity ---
if(0)%(exist('vehicle_angular_velocity', 'var') && exist('vehicle_rates_setpoint', 'var'))
    figure('Name', 'Rates', 'Color', 'w');
    % titles = {'Roll Rate', 'Pitch Rate', 'Yaw Rate'};
    ylabels = {'p (deg/s)', 'q (deg/s)', 'r (deg/s)'};
    ax = [];
    for i = 1:3
        ax(i) = subplot(3,1,i); hold on;step = 30;
        plot(vehicle_rates_setpoint_t(1:step:end), vehicle_rates_setpoint(1:step:end,i)*r2d, STYLE_SP{:});step = 30;
        plot(vehicle_angular_velocity_t(1:step:end), vehicle_angular_velocity(1:step:end,i)*r2d, STYLE_RES{:});
        grid on; ylabel(ylabels{i}); % title(titles{i});
        if i==1, legend('Setpoint','Response','NumColumns',2, 'Location', 'best'); end
        add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
        
    end
    linkaxes(ax, 'x'); xlabel('Time (s)');
    % PlotToFile(gcf, 'results/rates.pdf', 12, 5);
end

% --- Figure 2: Attitude ---
if(0)%(exist('Roll', 'var') && exist('Roll_setpoint', 'var'))
    figure('Name', 'Attitude', 'Color', 'w');step = 10;
    d_sp = {Roll_setpoint(1:step:end), Pitch_setpoint(1:step:end), Yaw_setpoint(1:step:end)};step = 1;
    d_res = {Roll(1:step:end), Pitch(1:step:end), Yaw(1:step:end)};
    % titles = {'Roll', 'Pitch', 'Yaw'};

    % \varphi corresponds to Roll, \theta corresponds to Pitch, \phi corresponds to Yaw
    ylabels = {'$\varphi$ (deg)', '$\theta$ (deg)', '$\phi$ (deg)'};

    ax = [];
    for i = 1:3
        ax(i) = subplot(3,1,i); hold on;step = 10;
        plot(vehicle_attitude_setpoint_t(1:step:end), d_sp{i}*r2d, STYLE_SP{:});step = 1;
        plot(vehicle_attitude_t(1:step:end), d_res{i}*r2d, STYLE_RES{:});
        grid on; ylabel(ylabels{i});  %title(titles{i}); 
        if i==1, legend('Setpoint','Response','NumColumns',2, 'Location', 'best'); end
        add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
    end
    linkaxes(ax, 'x'); xlabel('Time (s)');
    % PlotToFile(gcf, 'results/att.pdf', 12, 5);
end

% --- Figure 3a: Velocity & TECS (Dynamic subplot count: 3 or 4) ---
if (0)%exist('V_XYZ', 'var')
    figure('Name', 'Velocity', 'Color', 'w');  
    n_rows = 3;
    ax = [];

    % --- Subplot 1: Velocity X ---
    ax(1) = subplot(n_rows, 1, 1); hold on;
    if exist('V_XYZ_setpoint', 'var'), plot(vehicle_local_position_setpoint_t, V_XYZ_setpoint(:,1), STYLE_SP{:}); end
    plot(vehicle_local_position_t, V_XYZ(:,1), STYLE_RES{:});
    grid on; ylabel('$v_x$ (m/s)');   %title('Velocity X');
    add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

    % --- Subplot 2: Velocity Y ---
    ax(2) = subplot(n_rows, 1, 2); hold on;
    if exist('V_XYZ_setpoint', 'var'), plot(vehicle_local_position_setpoint_t, V_XYZ_setpoint(:,2), STYLE_SP{:}); end
    plot(vehicle_local_position_t, V_XYZ(:,2), STYLE_RES{:});
    grid on; ylabel('$v_y$ (m/s)');   %title('Velocity Y');
    add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

    % --- Subplot 3: Velocity Z ---
    ax(3) = subplot(n_rows, 1, 3); hold on;
    if exist('V_XYZ_setpoint', 'var'), plot(vehicle_local_position_setpoint_t, V_XYZ_setpoint(:,3), STYLE_SP{:}); end
    plot(vehicle_local_position_t, V_XYZ(:,3), STYLE_RES{:});
    if exist('V_XYZ_setpoint', 'var'), legend('Setpoint','Response','NumColumns',2, 'Location', 'best'); end
    grid on; ylabel('$v_z$ (m/s)'); xlabel('Time (s)'); %title('Velocity Z');
    add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);



    linkaxes(ax, 'x');   
    % PlotToFile(gcf, 'results/vel.pdf', 12, 5);
end

% --- Subplot 3b: TECS Height Rate (Only plotted if exists) ---
if (0)%exist('tecs_t', 'var')
    figure('Name', 'tecs', 'Color', 'w');  
    n_rows = 6;
    ax = [];
    
    ax(1) = subplot(n_rows, 1, 1); hold on;
    plot(tecs_t, altitude_sp, STYLE_SP{:});hold on;
    s = plot_style(sty, 'setpoint7', ':', 1.4);
    plot(tecs_t, altitude_reference, s{:});
    plot(vehicle_local_position_t, -XYZ(:,3)+ref_alt, STYLE_RES{:});hold on;
    grid on; ylabel('$h$ (m)');xlabel('Time (s)');title('TECS altitude');
    legend('altitude sp', 'altitude reference','altitude', 'Location', 'best');
    add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

    ax(2) = subplot(n_rows, 1, 2); hold on;
     
    s = plot_style(sty, 'setpoint4', '--', 1);
    plot(tecs_t, height_rate_reference, s{:});  
    s = plot_style(sty, 'response4', '-', 1.5);
    plot(tecs_t, height_rate_direct, s{:}); 
    plot(tecs_t, height_rate_setpoint, STYLE_SP{:}); 
    plot(tecs_t, height_rate, STYLE_RES{:});
    grid on; ylabel('$\dot{h}$ (m/s)'); xlabel('Time (s)'); title('TECS height rate');
    legend('height rate reference', 'height rate direct','height rate setpoint','height rate', 'Location', 'best');
    add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

    ax(3) = subplot(n_rows, 1, 3); hold on;
    % s = plot_style(sty, 'setpoint7', '--', 0.8);
    % plot(tecs_t, equivalent_airspeed_sp, s{:});
    plot(tecs_t, true_airspeed_sp, STYLE_SP{:});
    plot(tecs_t, true_airspeed_filtered, STYLE_RES{:});
    grid on; ylabel('$V_a$ (m/s)'); xlabel('Time (s)'); title('TECS airspeed');
    legend('true airspeed sp', 'true airspeed filtered', 'Location', 'best');
    add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

    ax(4) = subplot(n_rows, 1, 4); hold on;
    plot(tecs_t, true_airspeed_derivative_sp, STYLE_SP{:});
    plot(tecs_t, true_airspeed_derivative, STYLE_RES{:});
    s = plot_style(sty, 'response4', '-.', 1.1);
    % plot(tecs_t, true_airspeed_derivative_raw, s{:});
    grid on; ylabel('$\dot{V}_a$ (m/s)'); title('TECS airspeed derivative');
    legend('true airspeed derivative sp', 'true airspeed derivative', 'Location', 'best');
    add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

    ax(5) = subplot(n_rows, 1, 5); hold on;
    plot(tecs_t, total_energy_rate_sp, STYLE_SP{:});
    plot(tecs_t, total_energy_rate, STYLE_RES{:});
    grid on; ylabel('$\dot{E}$'); xlabel('Time (s)');title('TECS total energy rate');
    legend('total energy rate sp','total energy rate', 'Location', 'best');
    add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

    ax(6) = subplot(n_rows, 1, 6); hold on;
    plot(tecs_t, total_energy_balance_rate_sp, STYLE_SP{:});
    plot(tecs_t, total_energy_balance_rate, STYLE_RES{:});
    grid on; ylabel('$\dot{B}$'); xlabel('Time (s)'); title('TECS total energy balance rate');
    legend('total energy balance rate sp','total energy balance rate', 'Location', 'best');
    add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

    % PlotToFile(gcf, 'results/tecs.pdf', 20, 24);
end

% --- Figure 4: Position ---
if(0)%(exist('XYZ', 'var'))
    figure('Name', 'Position', 'Color', 'w'); 
    ylabels = {'x (m)', 'y (m)', 'z (m)'};
    has_sp = exist('XYZ_setpoint', 'var');
    ax = [];
    for i = 1:3
        ax(i) = subplot(3,1,i); hold on;
        if has_sp, plot(vehicle_local_position_setpoint_t, XYZ_setpoint(:,i), STYLE_SP{:}); end
        plot(vehicle_local_position_t, XYZ(:,i), STYLE_RES{:});
        grid on; ylabel(ylabels{i}); 
        if i==3
            % title('Position'); 
            legend('Setpoint','Response','NumColumns',2, 'Location', 'northeast'); 
        end
        add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
    end
    linkaxes(ax, 'x');   xlabel('Time (s)');
    % PlotToFile(gcf, 'results/pos.pdf', 12, 5);
end

%% =========================================================================
%  Trajectory  
% =========================================================================
%% =========================================================================
%  Tracking error in three algorithm windows
% =========================================================================
if(0)%(exist('XYZ', 'var') && exist('XYZ_setpoint', 'var') && ...
     %   exist('vehicle_local_position_t', 'var') && ...
     %   exist('vehicle_local_position_setpoint_t', 'var'))

    figure('Name', 'Tracking Error', 'Color', 'w'); hold on;

    % Start times are arranged as [PX4, PINDI, GINDI]
    start_time=30;
    segment_start_times = [start_time, start_time+35, start_time+35*2];      % s
    segment_duration    = 20;                % s
    segment_names       = {'PX4', 'PINDI', 'GINDI'};
    segment_styles      = {STYLE_AXIS1_BOLD, STYLE_AXIS2_BOLD, STYLE_AXIS3_BOLD};

    t_res = vehicle_local_position_t(:);
    t_sp  = vehicle_local_position_setpoint_t(:);

    [t_sp_unique, ia] = unique(t_sp, 'stable');
    XYZ_sp_unique = XYZ_setpoint(ia, :);

    XYZ_sp_interp = interp1(t_sp_unique, XYZ_sp_unique, t_res, 'linear', 'extrap');

    e_xyz  = XYZ - XYZ_sp_interp;
    e_norm = sqrt(sum(e_xyz.^2, 2));

    for k = 1:numel(segment_start_times)
        t0 = segment_start_times(k);
        t1 = t0 + segment_duration;
        idx = (t_res >= t0) & (t_res <= t1);

        if nnz(idx) >= 2
            t_seg = t_res(idx) - t0;
            s = segment_styles{k};
            plot(t_seg, e_norm(idx), s{:}, 'DisplayName', segment_names{k});
        end
    end

    grid on;
    xlabel('Time from segment start (s)');
    ylabel('Tracking error (m)');
    xlim([0, segment_duration]);
    legend('Location', 'best');
    % PlotToFile(gcf, 'results/tracking_error.pdf', 8, 3.5);
end


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
    segment_names       = {'PX4', 'PINDI', 'GINDI'};
    STYLE_1 = plot_style(sty, 'setpoint1', '--',  LW_MAIN);
    STYLE_2 = plot_style(sty, 'setpoint8', '-.', LW_MAIN);
    STYLE_3 = plot_style(sty, 'response4', '-', LW_MAIN);
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
            % lgd.ItemTokenSize = [8, 4];
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