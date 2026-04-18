clear all;
close all;
clc;

run('load_data_main.m'); % load('flight_data.mat'); 20_00_18
addpath(genpath(pwd));

%% =========================================================================
%  Global figure / font style
% =========================================================================
set(groot, ...
    'defaultAxesFontSize', 7, ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultAxesLineWidth', 0.5, ...
    'defaultAxesLabelFontSizeMultiplier', 1, ...
    'defaultAxesTitleFontSizeMultiplier', 1);

set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');

%% =========================================================================
%  Unified style manager
% =========================================================================
sty = make_plot_style();

% 为兼容你原先 Figure 1-5 的写法，保留这两个名字
STYLE_SP  = sty.setpoint;
STYLE_RES = sty.response;

%% =========================================================================
%  Parameter settings
% =========================================================================
MAX_MOTORS = 12; 
MAX_SERVOS = 8;  
plot_together = 0;   % 1: Servos and motors combined display, 0: independent display
verbose = 0;         % 1: Display more diagnostics
control_fig = 0;     % 1: Display control quantities
n_raw_plot = 8;      % For versions before 1.13, plot first 8 channels of pwm

%% =========================================================================
%  Figure 1 2, 3, 4, 5: (vehicle_angular_velocity, Attitude, Vel, Pos, Control)
% =========================================================================
% --- Figure 1: vehicle_angular_velocity ---
if(exist('vehicle_angular_velocity', 'var') && exist('vehicle_rates_setpoint', 'var'))
    figure('Name', 'Rates', 'Color', 'w');
    % titles = {'Roll Rate', 'Pitch Rate', 'Yaw Rate'};
    ylabels = {'p (deg/s)', 'q (deg/s)', 'r (deg/s)'};
    ax = [];
    for i = 1:3
        ax(i) = subplot(3,1,i); hold on;step = 30;
        plot(vehicle_rates_setpoint_t(1:step:end), vehicle_rates_setpoint(1:step:end,i)*r2d, STYLE_SP{:});step = 30;
        plot(vehicle_angular_velocity_t(1:step:end), vehicle_angular_velocity(1:step:end,i)*r2d, STYLE_RES{:});
        grid on; ylabel(ylabels{i}); % title(titles{i});
        if i==1, legend('Setpoint', 'Response', 'Location', 'best'); end
        add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
        % lines = findobj(gca, 'Type', 'line');
        % uistack(lines, 'top');
        
    end
    linkaxes(ax, 'x'); xlabel('Time (s)');
    % PlotToFile(gcf, 'results/rates.pdf', 12, 6.8);
end

% --- Figure 2: Attitude ---
if(exist('Roll', 'var') && exist('Roll_setpoint', 'var'))
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
        if i==1, legend('Setpoint','Response', 'Location', 'best'); end
        add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
    end
    linkaxes(ax, 'x'); xlabel('Time (s)');
    % PlotToFile(gcf, 'results/att.pdf', 12, 6.8);
end

% --- Figure 3: Velocity ---
% if(exist('V_XYZ', 'var'))
%     figure('Name', 'Velocity', 'Color', 'w'); 
%     ylabels = {'V_X (m/s)', 'V_Y (m/s)', 'V_Z (m/s)'};
%     has_sp = exist('V_XYZ_setpoint', 'var');
%     ax = [];
%     for i = 1:3
%         ax(i) = subplot(3,1,i); hold on;
%         if has_sp, plot(vehicle_local_position_setpoint_t, V_XYZ_setpoint(:,i), STYLE_SP{:}); end
%         plot(vehicle_local_position_t, V_XYZ(:,i), STYLE_RES{:});
%         grid on; ylabel(ylabels{i}); if i==1, title('Velocity'); legend('Setpoint','Response'); end
%         add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
%     end
%     linkaxes(ax, 'x'); xlabel('Time (s)');
% end

% --- Figure 3a: Velocity & TECS (Dynamic subplot count: 3 or 4) ---
if exist('V_XYZ', 'var')
    figure('Name', 'Velocity', 'Color', 'w');  
    n_rows = 3;
    ax = [];

    % --- Subplot 1: Velocity X ---
    ax(1) = subplot(n_rows, 1, 1); hold on;
    if exist('V_XYZ_setpoint', 'var'), plot(vehicle_local_position_setpoint_t, V_XYZ_setpoint(:,1), STYLE_SP{:}); end
    plot(vehicle_local_position_t, V_XYZ(:,1), STYLE_RES{:});
    grid on; ylabel('$v_x$ (m/s)');  %title('Velocity X');
    add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

    % --- Subplot 2: Velocity Y ---
    ax(2) = subplot(n_rows, 1, 2); hold on;
    if exist('V_XYZ_setpoint', 'var'), plot(vehicle_local_position_setpoint_t, V_XYZ_setpoint(:,2), STYLE_SP{:}); end
    plot(vehicle_local_position_t, V_XYZ(:,2), STYLE_RES{:});
    grid on; ylabel('$v_y$ (m/s)');  %title('Velocity Y');
    add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

    % --- Subplot 3: Velocity Z ---
    ax(3) = subplot(n_rows, 1, 3); hold on;
    if exist('V_XYZ_setpoint', 'var'), plot(vehicle_local_position_setpoint_t, V_XYZ_setpoint(:,3), STYLE_SP{:}); end
    plot(vehicle_local_position_t, V_XYZ(:,3), STYLE_RES{:});
    if exist('V_XYZ_setpoint', 'var'), legend('Setpoint', 'Response', 'Location', 'southwest','NumColumns',2); end
    grid on; ylabel('$v_z$ (m/s)');  %title('Velocity Z');
    add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);



    linkaxes(ax, 'x');  xlabel('Time (s)');
    % PlotToFile(gcf, 'results/vel.pdf', 12, 6.8);
end

% --- Subplot 3b: TECS Height Rate (Only plotted if exists) ---
if exist('tecs_t', 'var')
    figure('Name', 'tecs', 'Color', 'w');  
    n_rows = 5;
    ax = [];
    ax(1) = subplot(n_rows, 1, 1); hold on;
    % plot(tecs_t, altitude_sp);
    plot(tecs_t, height_rate_reference, sty.axis2_dash{:});
    plot(tecs_t, height_rate_direct);
    plot(tecs_t, height_rate_setpoint, STYLE_SP{:});
    plot(tecs_t, height_rate, STYLE_RES{:});
    grid on; ylabel('$\dot{h}$ (m/s)');  %title('TECS Height Rate');
    % legend('$\dot{h}_r$', '$\dot{h}$', 'Location', 'best');
    legend('height rate reference', 'height rate direct','height rate setpoint','height rate', 'Location', 'best');
    add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

    ax(2) = subplot(n_rows, 1, 2); hold on;
    plot(tecs_t, equivalent_airspeed_sp);
    plot(tecs_t, true_airspeed_sp, STYLE_SP{:});
    plot(tecs_t, true_airspeed_filtered, STYLE_RES{:});
    grid on; ylabel('$V_a$ (m/s)');  %title('TECS airspeed');
    % legend('$V_{a,d}$', '$V_a$', 'Location', 'best');
    legend('equivalent airspeed sp', 'true airspeed sp', 'true airspeed filtered', 'Location', 'best');
    add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

    ax(3) = subplot(n_rows, 1, 3); hold on;
    plot(tecs_t, true_airspeed_derivative_sp);
    plot(tecs_t, true_airspeed_derivative, STYLE_RES{:});
    plot(tecs_t, true_airspeed_derivative_raw);
    grid on; ylabel('$\dot{V}_a$ (m/s)'); xlabel('Time (s)'); %title('TECS airspeed_derivative');
    % legend('$V_{a,d}$', '$V_a$', 'Location', 'best');
    % legend('Setpoint', 'Response', 'Location', 'best','NumColumns',2);
    legend('true airspeed derivative sp', 'true airspeed derivative', 'true airspeed derivative raw', 'Location', 'best');
    add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

    ax(4) = subplot(n_rows, 1, 4); hold on;
    plot(tecs_t, total_energy_rate_sp, STYLE_RES{:});
    plot(tecs_t, total_energy_rate);
    grid on; ylabel('$\dot{E}$'); xlabel('Time (s)'); 
    % legend('Setpoint', 'Response', 'Location', 'best','NumColumns',2);
    legend('total_energy_rate_sp','total_energy_rate', 'Location', 'best');
    add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

    ax(5) = subplot(n_rows, 1, 5); hold on;
    plot(tecs_t, total_energy_balance_rate_sp, STYLE_RES{:});
    plot(tecs_t, total_energy_balance_rate);
    grid on; ylabel('$\dot{B}$'); xlabel('Time (s)'); 
    % legend('Setpoint', 'Response', 'Location', 'best','NumColumns',2);
    legend('total energy balance rate sp','total energy balance rate', 'Location', 'best');
    add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

    % PlotToFile(gcf, 'results/tecs.pdf', 12, 6.8);
end

% --- Figure 4: Position ---
if(exist('XYZ', 'var'))
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
            legend('Setpoint','Response', 'Location', 'best','NumColumns',2); 
        end
        add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
    end
    linkaxes(ax, 'x'); xlabel('Time (s)');
    % PlotToFile(gcf, 'results/pos.pdf', 12, 6.8);
end

% %% =========================================================================
% %  Trajectory  
% % =========================================================================
% if(exist('XYZ', 'var') && exist('XYZ_setpoint', 'var'))
%     figure('Name', 'Trajectory', 'Color', 'w');
%     step = 10;
%     plot3(XYZ_setpoint(1:step:end,2), XYZ_setpoint(1:step:end,1), -XYZ_setpoint(1:step:end,3), STYLE_SP{:}); hold on;
%     plot3(XYZ(:,2), XYZ(:,1), -XYZ(:,3), STYLE_RES{:});
%     % title('Trajectory'); 
%     xlabel('y (m)'); ylabel('x (m)'); zlabel('-z (m)'); grid on; view(45, 30);legend('Setpoint', 'Response', 'Location', [0.397851640472838 0.245923662405286 0.241278108465608 0.125550660792952]);
%     % PlotToFile(gcf, 'results/traj.pdf', 8, 3.5);
% end

%% =========================================================================
%  Trajectory colored by VTOL state
%  FW gust region: [100, 170] s uses darker color
% =========================================================================
if exist('XYZ', 'var') && exist('vehicle_local_position_t', 'var') && ...
   exist('vis_is_vtol', 'var') && vis_is_vtol && ...
   exist('vis_vtol_intervals', 'var') && ~isempty(vis_vtol_intervals)
   
    figure('Name', 'Trajectory by VTOL State', 'Color', 'w'); hold on;
    
    % --- setpoint trajectory ---
    if exist('XYZ_setpoint', 'var')
        step_sp = 10;
        plot3(XYZ_setpoint(1:step_sp:end,2), ...
              XYZ_setpoint(1:step_sp:end,1), ...
             -XYZ_setpoint(1:step_sp:end,3), ...
              STYLE_SP{:}, 'DisplayName', 'Setpoint');
    end
    
    % --- gust interval in FW ---
    gust_t0 = 100;
    gust_t1 = 170;
    
    % --- colors aligned with draw_status_band(...,'vtol_state') and sty ---
    c_hover    = [0.65, 0.75, 0.90]; % Hover 浅蓝
    c_trans    = [1.00, 0.80, 0.50]; % Trans 浅橙
    
    c_fw_pre   = [0.466, 0.674, 0.188]; % FW (横风前) - sty.c_axis3 绿色
    c_fw_gust  = [0.70,  0.22,  0.40];  % FW (横风段) - sty.c_response 紫红
    c_fw_post  = [0.494, 0.184, 0.556]; % FW (横风后) - sty.c_axis4 紫色
    
    t_pos = vehicle_local_position_t(:);
    
    % legend flags
    shown_hover    = false;
    shown_trans    = false;
    shown_fw_pre   = false;
    shown_fw_gust  = false;
    shown_fw_post  = false;
    
    for k = 1:size(vis_vtol_intervals, 1)
        t_s = vis_vtol_intervals(k, 1);
        t_e = vis_vtol_intervals(k, 2);
        s   = vis_vtol_intervals(k, 3);
        
        % ---------------- Hover ----------------
        if s == 1
            idx = find(t_pos >= t_s & t_pos <= t_e);
            if numel(idx) >= 2
                if ~shown_hover
                    plot3(XYZ(idx,2), XYZ(idx,1), -XYZ(idx,3), ...
                        'Color', c_hover, 'LineWidth', 0.5, 'LineStyle', '-', 'DisplayName', 'Hover');
                    shown_hover = true;
                else
                    plot3(XYZ(idx,2), XYZ(idx,1), -XYZ(idx,3), ...
                        'Color', c_hover, 'LineWidth', 0.5, 'LineStyle', '-', 'HandleVisibility', 'off');
                end
            end
            
        % ---------------- FW ----------------
        elseif s == 2
            % part A: FW before gust
            t1a = t_s;
            t1b = min(t_e, gust_t0);
            if t1b > t1a
                idx = find(t_pos >= t1a & t_pos <= t1b);
                if numel(idx) >= 2
                    if ~shown_fw_pre
                        plot3(XYZ(idx,2), XYZ(idx,1), -XYZ(idx,3), ...
                            'Color', c_fw_pre, 'LineWidth', 0.5, 'LineStyle', '--', 'DisplayName', 'FW (Pre-gust)');
                        shown_fw_pre = true;
                    else
                        plot3(XYZ(idx,2), XYZ(idx,1), -XYZ(idx,3), ...
                            'Color', c_fw_pre, 'LineWidth', 0.5, 'LineStyle', '--', 'HandleVisibility', 'off');
                    end
                end
            end
            
            % part B: FW with gust
            t2a = max(t_s, gust_t0);
            t2b = min(t_e, gust_t1);
            if t2b > t2a
                idx = find(t_pos >= t2a & t_pos <= t2b);
                if numel(idx) >= 2
                    if ~shown_fw_gust
                        plot3(XYZ(idx,2), XYZ(idx,1), -XYZ(idx,3), ...
                            'Color', c_fw_gust, 'LineWidth', 0.8, 'LineStyle', '-', 'DisplayName', 'FW (Gust)');
                        shown_fw_gust = true;
                    else
                        plot3(XYZ(idx,2), XYZ(idx,1), -XYZ(idx,3), ...
                            'Color', c_fw_gust, 'LineWidth', 0.8, 'LineStyle', '-', 'HandleVisibility', 'off');
                    end
                end
            end
            
            % part C: FW after gust
            t3a = max(t_s, gust_t1);
            t3b = t_e;
            if t3b > t3a
                idx = find(t_pos >= t3a & t_pos <= t3b);
                if numel(idx) >= 2
                    if ~shown_fw_post
                        plot3(XYZ(idx,2), XYZ(idx,1), -XYZ(idx,3), ...
                            'Color', c_fw_post, 'LineWidth', 0.5, 'LineStyle', '--', 'DisplayName', 'FW (Post-gust)');
                        shown_fw_post = true;
                    else
                        plot3(XYZ(idx,2), XYZ(idx,1), -XYZ(idx,3), ...
                            'Color', c_fw_post, 'LineWidth', 0.5, 'LineStyle', '--', 'HandleVisibility', 'off');
                    end
                end
            end
            
        % ---------------- Transition ----------------
        elseif s == 3
            idx = find(t_pos >= t_s & t_pos <= t_e);
            if numel(idx) >= 2
                if ~shown_trans
                    plot3(XYZ(idx,2), XYZ(idx,1), -XYZ(idx,3), ...
                        'Color', c_trans, 'LineWidth', 1, 'LineStyle', '-', 'DisplayName', 'Trans');
                    shown_trans = true;
                else
                    plot3(XYZ(idx,2), XYZ(idx,1), -XYZ(idx,3), ...
                        'Color', c_trans, 'LineWidth', 1, 'LineStyle', '-', 'HandleVisibility', 'off');
                end
            end
        end
    end
    
    % % optional markers
    % plot3(XYZ(1,2),   XYZ(1,1),   -XYZ(1,3),   'ko', 'MarkerSize', 5, ...
    %     'MarkerFaceColor', 'k', 'DisplayName', 'Start');
    % 
    % % 终点使用灰色方块，明确区分起点且视觉上不显得突兀
    % plot3(XYZ(end,2), XYZ(end,1), -XYZ(end,3), 's', 'MarkerSize', 6, ...
    %     'MarkerEdgeColor', [0.4, 0.4, 0.4], ...
    %     'MarkerFaceColor', [0.7, 0.7, 0.7], ...
    %     'DisplayName', 'End');
    % start / end
    plot3(XYZ(1,2), XYZ(1,1), -XYZ(1,3), 'ko', 'MarkerSize', 5, ...
        'MarkerFaceColor', 'k', 'DisplayName', 'Start');

    plot3(XYZ(end,2), XYZ(end,1), -XYZ(end,3), 's', 'MarkerSize', 6, ...
        'MarkerEdgeColor', [0.4, 0.4, 0.4], ...
        'MarkerFaceColor', [0.7, 0.7, 0.7], ...
        'DisplayName', 'End');

    % gust entry / exit
    [~, idx_gust_in]  = min(abs(t_pos - gust_t0));
    [~, idx_gust_out] = min(abs(t_pos - gust_t1));

    plot3(XYZ(idx_gust_in,2), XYZ(idx_gust_in,1), -XYZ(idx_gust_in,3), ...
        'd', 'MarkerSize', 7, ...
        'MarkerEdgeColor', [0.25, 0.25, 0.25], ...
        'MarkerFaceColor', c_fw_gust, ...
        'LineWidth', 0.8, ...
        'DisplayName', 'Gust entry');

    plot3(XYZ(idx_gust_out,2), XYZ(idx_gust_out,1), -XYZ(idx_gust_out,3), ...
        'p', 'MarkerSize', 8, ...
        'MarkerEdgeColor', [0.25, 0.25, 0.25], ...
        'MarkerFaceColor', c_fw_post, ...
        'LineWidth', 0.8, ...
        'DisplayName', 'Gust exit');

    xlabel('y (m)');
    ylabel('x (m)');
    zlabel('-z (m)');
    grid on;
    view([31.6 26.3]);
    legend('Location', 'eastoutside');
    % PlotToFile(gcf, 'results/traj.png', 12, 5);
end



%% =========================================================================
%  Control figures
% =========================================================================
if control_fig
    % ---------------------------------------------------------------------
    % Actuator Controls
    % ---------------------------------------------------------------------
    if ~isempty(actuator_controls_0.time)
        figure('Name', 'Actuator Controls', 'Color', 'w');

        has_group1 = ~isempty(actuator_controls_1.time);
        if has_group1
            n_cols = 2;
            layout_title = 'Actuator Controls (Left: Group 0, Right: Group 1)';
        else
            n_cols = 1;
            layout_title = 'Actuator Controls (Group 0)';
        end

        line_cols = {sty.axis1, sty.axis2, sty.axis3};
        titles = {'All', 'Roll', 'Pitch', 'Yaw', 'Thrust'};
        ax_all = [];

        for row = 1:5
            idx_left = (row - 1) * n_cols + 1;
            ax = subplot(5, n_cols, idx_left); hold on;
            ax_all = [ax_all, ax];

            if row == 1
                plot(actuator_controls_0.time, actuator_controls_0.roll,  line_cols{1}{:}, 'DisplayName', 'Roll');
                plot(actuator_controls_0.time, actuator_controls_0.pitch, line_cols{2}{:}, 'DisplayName', 'Pitch');
                plot(actuator_controls_0.time, actuator_controls_0.yaw,   line_cols{3}{:}, 'DisplayName', 'Yaw');

                if ~isempty(actuator_controls_0.thrust_z_neg)
                    plot(actuator_controls_0.time_thrust, actuator_controls_0.thrust_z_neg, ...
                        sty.thrust{:}, 'DisplayName', 'Thrust (up)');
                end
                if ~isempty(actuator_controls_0.thrust_x) && any(actuator_controls_0.thrust_x ~= 0)
                    plot(actuator_controls_0.time_thrust, actuator_controls_0.thrust_x, ...
                        sty.thrust_alt{:}, 'DisplayName', 'Thrust (fwd)');
                end
                title('Group 0 (Main)');
                if n_cols == 1
                    legend('Location', 'best', 'NumColumns', 5);
                end

            elseif row >= 2 && row <= 4
                data_map = {actuator_controls_0.roll, actuator_controls_0.pitch, actuator_controls_0.yaw};
                plot(actuator_controls_0.time, data_map{row-1}, line_cols{row-1}{:});
                ylabel(titles{row});

            elseif row == 5
                if ~isempty(actuator_controls_0.thrust_z_neg)
                    plot(actuator_controls_0.time_thrust, actuator_controls_0.thrust_z_neg, ...
                        sty.thrust{:}, 'DisplayName', 'Up');
                end
                if ~isempty(actuator_controls_0.thrust_x)
                    plot(actuator_controls_0.time_thrust, actuator_controls_0.thrust_x, ...
                        sty.thrust_alt{:}, 'DisplayName', 'Fwd');
                end
                ylabel('Thrust');
                legend('show', 'Location', 'best');
            end
            grid on;
            add_standard_background(vis_flight_intervals, vis_flight_names, ...
                vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

            if has_group1
                idx_right = (row - 1) * n_cols + 2;
                ax = subplot(5, n_cols, idx_right); hold on;
                ax_all = [ax_all, ax];

                if row == 1
                    plot(actuator_controls_1.time, actuator_controls_1.roll,  line_cols{1}{:}, 'DisplayName', 'Roll');
                    plot(actuator_controls_1.time, actuator_controls_1.pitch, line_cols{2}{:}, 'DisplayName', 'Pitch');
                    plot(actuator_controls_1.time, actuator_controls_1.yaw,   line_cols{3}{:}, 'DisplayName', 'Yaw');
                    if ~isempty(actuator_controls_1.thrust_x)
                        plot(actuator_controls_1.time, actuator_controls_1.thrust_x, ...
                            sty.thrust_alt{:}, 'DisplayName', 'Thrust (fwd)');
                    end
                    title('Group 1 (Aux/FW)');

                elseif row >= 2 && row <= 4
                    data_map = {actuator_controls_1.roll, actuator_controls_1.pitch, actuator_controls_1.yaw};
                    plot(actuator_controls_1.time, data_map{row-1}, line_cols{row-1}{:});

                elseif row == 5
                    if ~isempty(actuator_controls_1.thrust_x)
                        plot(actuator_controls_1.time, actuator_controls_1.thrust_x, ...
                            sty.thrust_alt{:}, 'DisplayName', 'Fwd');
                    end
                end
                grid on;
                add_standard_background(vis_flight_intervals, vis_flight_names, ...
                    vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
            end
        end

        linkaxes(ax_all, 'x');
        xlabel(ax_all(end), 'Time (s)');
        sgtitle(layout_title);
        % PlotToFile(gcf, 'results/control.png', 20, 20);
    end

    % ---------------------------------------------------------------------
    % Figure 6/7/8: Actuators & PWM
    % ---------------------------------------------------------------------
    if dynamic_control_alloc

        if plot_together
            if isfield(log.data, 'actuator_motors_0') || isfield(log.data, 'actuator_servos_0')

                n_motors = 0;
                if isfield(log, 'params') && isfield(log.params, 'CA_ROTOR_COUNT')
                    n_motors = double(log.params.CA_ROTOR_COUNT);
                else
                    if exist('motors', 'var')
                        n_motors = size(motors, 2);
                    end
                end

                n_servos = 0;
                if isfield(log, 'params') && isfield(log.params, 'CA_SV_CS_COUNT')
                    n_servos = double(log.params.CA_SV_CS_COUNT);
                else
                    if exist('servos', 'var')
                        n_servos = size(servos, 2);
                    end
                end

                n_motors = min(n_motors, MAX_MOTORS);
                n_servos = min(n_servos, MAX_SERVOS);

                has_motor_data = (n_motors > 0) && isfield(log.data, 'actuator_motors_0') && exist('motors', 'var');
                has_servo_data = (n_servos > 0) && isfield(log.data, 'actuator_servos_0') && exist('servos', 'var');

                total_subplots = double(has_motor_data) + double(has_servo_data);

                if total_subplots > 0
                    figure('Name', 'Actuator Outputs (Merged)', 'Color', 'w');
                    current_plot_idx = 1;

                    if has_motor_data
                        ax_m = subplot(total_subplots, 1, current_plot_idx); hold on;
                        colors_motor = make_channel_colors(n_motors, sty);
                        t_m = log.data.actuator_motors_0.timestamp * 1e-6;
                        for i = 1:n_motors
                            if i <= size(motors, 2)
                                plot(t_m, motors(:, i), 'Color', colors_motor(i,:), ...
                                    'LineWidth', sty.lw_multi, ...
                                    'DisplayName', sprintf('Motor %d', i));
                            end
                        end
                        grid on;
                        ylabel('Motors Output');
                        title(sprintf('Actuator: Motors (Total %d)', n_motors));
                        legend('show');
                        add_standard_background(vis_flight_intervals, vis_flight_names, ...
                            vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

                        if has_servo_data
                            set(gca, 'XTickLabel', []);
                        else
                            xlabel('Time (s)');
                        end
                        current_plot_idx = current_plot_idx + 1;
                    end

                    if has_servo_data
                        ax_s = subplot(total_subplots, 1, current_plot_idx); hold on;
                        colors_servo = make_channel_colors(n_servos, sty);
                        t_s = log.data.actuator_servos_0.timestamp * 1e-6;
                        for i = 1:n_servos
                            if i <= size(servos, 2)
                                plot(t_s, servos(:, i), 'Color', colors_servo(i,:), ...
                                    'LineWidth', sty.lw_multi, ...
                                    'DisplayName', sprintf('Servo %d', i));
                            end
                        end
                        grid on;
                        ylabel('Servos Output');
                        title(sprintf('Actuator: Servos (Total %d)', n_servos));
                        legend('show');
                        xlabel('Time (s)');
                        add_standard_background(vis_flight_intervals, vis_flight_names, ...
                            vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
                    end

                    if has_motor_data && has_servo_data
                        linkaxes([ax_m, ax_s], 'x');
                    end
                    % PlotToFile(gcf, 'results/motors_servos.png', 20, 20);
                end
            end

        else
            % -----------------------------------------------------------------
            % Figure 6: Motors
            % -----------------------------------------------------------------
            if isfield(log.data, 'actuator_motors_0') && exist('motors', 'var')

                if isfield(log, 'params') && isfield(log.params, 'CA_ROTOR_COUNT')
                    n_motors = double(log.params.CA_ROTOR_COUNT);
                else
                    n_motors = size(motors, 2);
                end
                n_motors = min(n_motors, MAX_MOTORS);

                if n_motors > 0
                    figure('Name', 'Actuator Motors', 'Color', 'w');
                    motor_colors = make_channel_colors(n_motors, sty);

                    for i = 1:n_motors
                        subplot(n_motors, 1, i); hold on;

                        if i <= size(motors, 2)
                            plot(log.data.actuator_motors_0.timestamp*1e-6, motors(:, i), ...
                                'Color', motor_colors(i,:), 'LineWidth', sty.lw_multi);
                        end

                        grid on;
                        ylabel(sprintf('Motor %d', i));
                        add_standard_background(vis_flight_intervals, vis_flight_names, ...
                            vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

                        if i == 1
                            title(sprintf('Actuator Motors (Total: %d)', n_motors));
                        end
                        if i == n_motors
                            xlabel('Time (s)');
                        else
                            set(gca, 'XTickLabel', []);
                        end
                    end
                    % PlotToFile(gcf, 'results/motors.png', 20, 20);
                end
            end

            % -----------------------------------------------------------------
            % Figure 7: Servos
            % -----------------------------------------------------------------
            if isfield(log.data, 'actuator_servos_0') && exist('servos', 'var')

                if isfield(log, 'params') && isfield(log.params, 'CA_SV_CS_COUNT')
                    n_servos = double(log.params.CA_SV_CS_COUNT);
                else
                    n_servos = size(servos, 2);
                end
                n_servos = min(n_servos, MAX_SERVOS);

                if n_servos > 0
                    figure('Name', 'Actuator Servos', 'Color', 'w');
                    servo_colors = make_channel_colors(n_servos, sty);

                    for i = 1:n_servos
                        subplot(n_servos, 1, i); hold on;

                        if i <= size(servos, 2)
                            plot(log.data.actuator_servos_0.timestamp*1e-6, servos(:, i), ...
                                'Color', servo_colors(i,:), 'LineWidth', sty.lw_multi);
                        end

                        grid on;
                        ylabel(sprintf('Servo %d', i));
                        add_standard_background(vis_flight_intervals, vis_flight_names, ...
                            vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

                        if i == 1
                            title(sprintf('Actuator Servos (Total: %d)', n_servos));
                        end
                        if i == n_servos
                            xlabel('Time (s)');
                        else
                            set(gca, 'XTickLabel', []);
                        end
                    end
                    % PlotToFile(gcf, 'results/servos.png', 20, 20);
                end
            end
        end

        % -----------------------------------------------------------------
        % Figure 8: PWM Outputs
        % -----------------------------------------------------------------
        if exist('active_channels', 'var') && ~isempty(active_channels)
            figure('Name', 'Actuator Outputs (PWM)', 'Color', 'w');

            is_motor = strcmp({active_channels.type}, 'Motor');
            idx_m = find(is_motor);
            idx_s = find(~is_motor);

            n_plots = double(~isempty(idx_m)) + double(~isempty(idx_s));
            cur = 1;
            ax_list = [];

            if ~isempty(idx_m)
                ax = subplot(n_plots, 1, cur); hold on;
                cols = make_channel_colors(length(idx_m), sty);
                for k = 1:length(idx_m)
                    info = active_channels(idx_m(k));
                    data = log.data.actuator_outputs_0.(info.col_name);
                    plot(outputs_t, data, 'Color', cols(k,:), 'LineWidth', sty.lw_multi, ...
                        'DisplayName', sprintf('%s (Ch%d)', info.name, info.idx));
                end
                grid on; ylabel('PWM (us)'); title('Motors PWM'); legend('show');
                add_standard_background(vis_flight_intervals, vis_flight_names, ...
                    vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
                ax_list = [ax_list, ax];
                cur = cur + 1;
            end

            if ~isempty(idx_s)
                ax = subplot(n_plots, 1, cur); hold on;
                cols = make_channel_colors(length(idx_s), sty);
                for k = 1:length(idx_s)
                    info = active_channels(idx_s(k));
                    data = log.data.actuator_outputs_0.(info.col_name);
                    plot(outputs_t, data, 'Color', cols(k,:), 'LineWidth', sty.lw_multi, ...
                        'DisplayName', sprintf('%s (Ch%d)', info.name, info.idx));
                end
                grid on; ylabel('PWM (us)'); title('Servos/Other PWM'); legend('show');
                add_standard_background(vis_flight_intervals, vis_flight_names, ...
                    vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
                ax_list = [ax_list, ax];
            end

            if ~isempty(ax_list)
                linkaxes(ax_list, 'x');
            end
            xlabel('Time (s)');
            % PlotToFile(gcf, 'results/PWM.png', 20, 20);
        else
            fprintf('Figure 8 Skipped: No active PWM channels identified.\n');
        end

    else
        % -----------------------------------------------------------------
        % Original plotting fallback
        % -----------------------------------------------------------------
        figure('Name', 'Actuator Outputs (Raw)', 'Color', 'w');

        ax_raw = gobjects(1, n_raw_plot);
        raw_colors = make_channel_colors(n_raw_plot, sty);

        for i = 1:n_raw_plot
            ax_raw(i) = subplot(n_raw_plot, 1, i); hold on;

            y_data = outputs(:, i);

            if ~all(isnan(y_data)) && any(y_data ~= 0)
                plot(outputs_t, y_data, 'Color', raw_colors(i,:), 'LineWidth', sty.lw_multi);
            end

            ylabel(sprintf('Out %d', i-1));
            grid on;

            if i == 1
                title(sprintf('Actuator Outputs (First %d Raw)', n_raw_plot));
            end
            if i < n_raw_plot
                set(gca, 'XTickLabel', []);
            else
                xlabel('Time (s)');
            end

            axis tight;
            add_standard_background(vis_flight_intervals, vis_flight_names, ...
                vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
        end

        linkaxes(ax_raw, 'x');
        % PlotToFile(gcf, 'results/pwm.png', 20, 20);
    end
end

%% =========================================================================
%  Verbose figures
% =========================================================================
if verbose

    % ---------------------------------------------------------------------
    % Angular Accel
    % ---------------------------------------------------------------------
    if exist('vehicle_angular_acceleration', 'var') || ...
       (isfield(log.data, 'vehicle_angular_velocity_0') && ...
        ismember('xyz_derivative_0_', log.data.vehicle_angular_velocity_0.Properties.VariableNames))

        figure('Name', 'Ang Acc vs rates', 'Color', 'w');
        titles = {'Roll Ang Acc', 'Pitch Ang Acc', 'Yaw Ang Acc'};
        ax = gobjects(1,3);

        for i = 1:3
            ax(i) = subplot(3,1,i); hold on;
            plot(vehicle_angular_acceleration_t, vehicle_angular_acceleration(:,i), sty.aux1{:});
            plot(vehicle_angular_velocity_t, vehicle_angular_velocity(:,i), sty.response_fade{:});
            grid on; ylabel('rad/s$^2$'); title(titles{i});
            add_standard_background(vis_flight_intervals, vis_flight_names, ...
                vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
        end
        linkaxes(ax, 'x'); xlabel('Time (s)');
        % PlotToFile(gcf, 'results/Angular_Accel.png', 20, 10);
    end

    % ---------------------------------------------------------------------
    % Manual Control Inputs
    % ---------------------------------------------------------------------
    if isfield(log.data, 'manual_control_setpoint_0')
        figure('Name', 'Manual Control Inputs', 'Color', 'w'); hold on;
        plot(rc_t, rc_roll,     sty.axis1_bold{:}, 'DisplayName', 'Roll');
        plot(rc_t, rc_pitch,    sty.axis2_bold{:}, 'DisplayName', 'Pitch');
        plot(rc_t, rc_yaw,      sty.axis3_bold{:}, 'DisplayName', 'Yaw');
        plot(rc_t, rc_throttle, sty.thrust_bold{:}, 'DisplayName', 'Throttle');

        grid on;
        legend('show', 'Location', 'best', 'NumColumns', 4);
        ylabel('Norm Input [-1, 1]');
        title('Manual Control Inputs (Sticks)');
        xlabel('Time (s)');
        ylim([-1.1, 1.1]);
        add_standard_background(vis_flight_intervals, vis_flight_names, ...
            vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
        % PlotToFile(gcf, 'results/Manual_input.png', 20, 10);
    end

    % ---------------------------------------------------------------------
    % FFT
    % ---------------------------------------------------------------------
    if exist('actuator_controls_0', 'var')
        figure('Color', 'w', 'Name', 'Actuator Controls FFT');
        ctrl_data = [actuator_controls_0.roll, actuator_controls_0.pitch, actuator_controls_0.yaw];
        markers = {
            'MC_DTERM_CUTOFF',  'D-Term Cutoff';
            'IMU_DGYRO_CUTOFF', 'D-Gyro Cutoff';
            'IMU_GYRO_CUTOFF',  'Gyro Cutoff'
        };
        draw_fft_analysis(actuator_controls_0.time, ctrl_data, {'Roll', 'Pitch', 'Yaw'}, ...
            'Actuator Controls FFT (Torque Setpoint)', log.params, markers);
        % PlotToFile(gcf, 'results/control_fft.png', 20, 10);
    end

    if exist('vehicle_angular_velocity', 'var')
        figure('Color', 'w', 'Name', 'Angular Velocity FFT');
        markers = {
            'IMU_GYRO_CUTOFF',  'Gyro Cutoff';
            'IMU_GYRO_NF_FREQ', 'Notch Freq'
        };
        draw_fft_analysis(vehicle_angular_velocity_t, vehicle_angular_velocity, ...
            {'Rollrate', 'Pitchrate', 'Yawrate'}, 'Angular Velocity FFT', log.params, markers);
        % PlotToFile(gcf, 'results/rate_fft.png', 20, 10);
    end

    if exist('vehicle_angular_acceleration', 'var')
        figure('Color', 'w', 'Name', 'Angular Acceleration FFT');
        markers = {
            'IMU_DGYRO_CUTOFF', 'D-Gyro Cutoff';
            'IMU_GYRO_NF_FREQ', 'Notch Freq'
        };
        draw_fft_analysis(vehicle_angular_acceleration_t, vehicle_angular_acceleration, ...
            {'Roll Acc', 'Pitch Acc', 'Yaw Acc'}, 'Angular Acceleration FFT', log.params, markers);
        % PlotToFile(gcf, 'results/acc_fft.png', 20, 10);
    end

    % ---------------------------------------------------------------------
    % Raw Acceleration
    % ---------------------------------------------------------------------
    if exist('raw_acc', 'var')
        figure('Name', 'Raw Acceleration', 'Color', 'w'); hold on;
        plot(raw_acc_t, raw_acc(:,1), sty.axis1_thin{:}, 'DisplayName', 'Acc X');
        plot(raw_acc_t, raw_acc(:,2), sty.axis2_thin{:}, 'DisplayName', 'Acc Y');
        plot(raw_acc_t, raw_acc(:,3), sty.axis3_thin{:}, 'DisplayName', 'Acc Z');

        grid on;
        legend('show', 'Location', 'best');
        ylabel('Acceleration [m/s$^2$]');
        title('Raw Acceleration');
        xlabel('Time (s)');
        add_standard_background(vis_flight_intervals, vis_flight_names, ...
            vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
        % PlotToFile(gcf, 'results/raw_acc.png', 20, 10);
    end

    % ---------------------------------------------------------------------
    % Vibration Metrics
    % ---------------------------------------------------------------------
    if exist('vib_data', 'var') && ~isempty(vib_data)
        figure('Name', 'Vibration Metrics', 'Color', 'w'); hold on;
        line_colors = make_channel_colors(length(vib_data), sty);
        for k = 1:length(vib_data)
            id = vib_data(k).id;
            plot(vib_data(k).t, vib_data(k).val, ...
                'Color', line_colors(k,:), ...
                'LineWidth', sty.lw_multi_bold, ...
                'DisplayName', sprintf('Accel %d Vib [m/s$^2$]', id));
        end
        grid on;
        ylabel('Vibration Level [m/s$^2$]');
        title('Vibration Metrics');
        xlabel('Time (s)');
        legend('show', 'Location', 'best');
        add_standard_background(vis_flight_intervals, vis_flight_names, ...
            vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
        set(gca, 'Layer', 'top');
        % PlotToFile(gcf, 'results/Vibration_Metrics.png', 20, 10);
    end

    % ---------------------------------------------------------------------
    % Spectrogram
    % ---------------------------------------------------------------------
    if exist('raw_acc', 'var') && ~isempty(raw_acc)
        figure('Name', 'Accel PSD', 'Color', 'w');
        draw_spec_analysis(raw_acc_t, raw_acc, 'Acceleration Power Spectral Density (Sum X+Y+Z)');
        % PlotToFile(gcf, 'results/Acceleration_Spectrogram.png', 20, 10);
    end

    if exist('vehicle_angular_velocity', 'var')
        figure('Name', 'Gyro PSD', 'Color', 'w');
        draw_spec_analysis(vehicle_angular_velocity_t, vehicle_angular_velocity, ...
            'Angular Velocity PSD (Sum X+Y+Z)');
        % PlotToFile(gcf, 'results/Gyro_Spectrogram.png', 20, 10);
    end

    if exist('vehicle_angular_acceleration', 'var')
        figure('Name', 'AngAcc PSD', 'Color', 'w');
        draw_spec_analysis(vehicle_angular_acceleration_t, vehicle_angular_acceleration, ...
            'Angular Acceleration PSD (Sum X+Y+Z)');
        % PlotToFile(gcf, 'results/AngAcc_Spectrogram.png', 20, 10);
    end

    % ---------------------------------------------------------------------
    % Raw Angular Speed
    % ---------------------------------------------------------------------
    if exist('raw_gyro', 'var')
        figure('Name', 'Raw Angular Speed', 'Color', 'w'); hold on;
        plot(raw_gyro_t, raw_gyro(:,1), sty.axis1_thin{:}, 'DisplayName', 'X');
        plot(raw_gyro_t, raw_gyro(:,2), sty.axis2_thin{:}, 'DisplayName', 'Y');
        plot(raw_gyro_t, raw_gyro(:,3), sty.axis3_thin{:}, 'DisplayName', 'Z');

        grid on; legend('show');
        ylabel('Angular Speed [deg/s]');
        title('Raw Angular Speed (Gyroscope)');
        xlabel('Time (s)');
        add_standard_background(vis_flight_intervals, vis_flight_names, ...
            vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
        % PlotToFile(gcf, 'results/raw_gyro.png', 20, 10);
    end

    % ToDO:
    % %% =========================================================================
    % %  FIFO Acceleration Analysis
    % %  Corresponds to Python: FIFO accel (Raw, PSD, Sampling Regularity)
    % % =========================================================================
    % if exist('fifo_acc', 'var')
    %     for k = 1:length(fifo_acc)
    %         id = fifo_acc(k).id;
    %         t = fifo_acc(k).t;
    %         d = fifo_acc(k).d;
    % 
    %         % 1. Raw Acceleration (FIFO)
    %         f_num = 20 + k*3 - 2; % Dynamic figure numbering: 21, 24, 27...
    %         figure(f_num); set(gcf, 'Color', 'w', 'Name', sprintf('FIFO Accel %d Raw', id));
    %         hold on;
    %         plot(t, d(:,1), 'r-', 'LineWidth', 0.1);
    %         plot(t, d(:,2), 'k-', 'LineWidth', 0.1);
    %         plot(t, d(:,3), 'b-', 'LineWidth', 0.1);
    %         title(sprintf('Raw Acceleration (FIFO, IMU%d)', id));
    %         ylabel('[m/s^2]'); grid on; legend('X','Y','Z');
    %         add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
    % 
    %         % 2. PSD (Spectrogram)
    %         % Call previous draw_spec_analysis
    %         % Note: pass microsecond timestamps or let the function handle seconds
    %         figure(f_num+1); set(gcf, 'Color', 'w', 'Name', sprintf('FIFO Accel %d PSD', id));
    %         draw_spec_analysis(t, d, sprintf('Acceleration PSD (FIFO, IMU%d)', id));
    % 
    %         % 3. Sampling Regularity
    %         % Python logic: diff(raw_message_timestamps)
    %         % Check packet arrival intervals rather than sample intervals
    %         figure(f_num+2); set(gcf, 'Color', 'w', 'Name', sprintf('FIFO Accel %d Sampling', id));
    %         if length(fifo_acc(k).raw_t) > 1
    %             dt_diff = diff(fifo_acc(k).raw_t) * 1e6; % Convert to microseconds
    %             plot(fifo_acc(k).raw_t(2:end), dt_diff, 'b.-');
    %             ylabel('Delta t [us]'); title(sprintf('Sampling Regularity (FIFO, IMU%d)', id));
    %             xlabel('Time (s)'); grid on;
    % 
    %             % Mark dropouts (Python logic: plot_dropouts)
    %             % Draw a reference line, e.g., 2000us (500Hz) or 1000us (1kHz)
    %             avg_dt = median(dt_diff);
    %             yline(avg_dt, 'r--', sprintf('Median: %.0f us', avg_dt));
    %         end
    %         add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
    %     end
    % end
    % 
    % %% =========================================================================
    % %  FIFO Gyro Analysis
    % %  Corresponds to Python: FIFO gyro (Raw, PSD)
    % % =========================================================================
    % if exist('fifo_gyro', 'var')
    %     % To avoid figure number conflicts above, start from 40
    %     start_fig = 40;
    %     for k = 1:length(fifo_gyro)
    %         id = fifo_gyro(k).id;
    %         t = fifo_gyro(k).t;
    %         d = fifo_gyro(k).d; % Already in deg/s
    % 
    %         % 1. Raw Gyro (FIFO)
    %         f_num = start_fig + (k-1)*2;
    %         figure(f_num); set(gcf, 'Color', 'w', 'Name', sprintf('FIFO Gyro %d Raw', id));
    %         hold on;
    %         plot(t, d(:,1), 'r-', 'LineWidth', 0.1);
    %         plot(t, d(:,2), 'k-', 'LineWidth', 0.1);
    %         plot(t, d(:,3), 'b-', 'LineWidth', 0.1);
    %         title(sprintf('Raw Gyro (FIFO, IMU%d)', id));
    %         ylabel('[deg/s]'); grid on; legend('X','Y','Z');
    %         add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
    % 
    %         % 2. PSD
    %         figure(f_num+1); set(gcf, 'Color', 'w', 'Name', sprintf('FIFO Gyro %d PSD', id));
    %         % PSD: use original units (rad/s) or deg/s? Python computes PSD without unit conversion,
    %         % and vibration metrics usually care about frequency distribution, so consistent units are fine.
    %         draw_spec_analysis(t, d, sprintf('Gyro PSD (FIFO, IMU%d)', id));
    %     end
    % end
    
    
    
    %% =========================================================================
    %  Raw Magnetic Field Strength
    %  Corresponds to Python: magnetometer_ga_topic
    % =========================================================================
    if exist('mag_data', 'var')
        figure('Name', 'Raw Magnetic Field', 'Color', 'w'); hold on;
        plot(mag_t, mag_data(:,1), sty.axis1{:}, 'DisplayName', 'X');
        plot(mag_t, mag_data(:,2), sty.axis2{:}, 'DisplayName', 'Y');
        plot(mag_t, mag_data(:,3), sty.axis3{:}, 'DisplayName', 'Z');

        grid on; legend('show');
        ylabel('Magnetic Field [Gauss]');
        title('Raw Magnetic Field Strength');
        xlabel('Time (s)');
        add_standard_background(vis_flight_intervals, vis_flight_names, ...
            vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
        % PlotToFile(gcf, 'results/Raw_Magnetic.png', 20, 10);
    end

    % ---------------------------------------------------------------------
    % Distance Sensor
    % ---------------------------------------------------------------------
    if exist('dist_val', 'var') || exist('dist_bottom', 'var')
        figure('Name', 'Distance Sensor', 'Color', 'w'); hold on;

        if exist('dist_val', 'var')
            plot(dist_sensor_t, dist_val, sty.response{:}, 'DisplayName', 'Distance');
        end
        if exist('dist_bottom', 'var')
            plot(dist_bottom_t, dist_bottom, sty.setpoint{:}, 'DisplayName', 'Est. Dist Bottom');
        end

        grid on; legend('show');
        ylabel('Distance [m]');
        title('Distance Sensor');
        xlabel('Time (s)');
        add_standard_background(vis_flight_intervals, vis_flight_names, ...
            vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
        % PlotToFile(gcf, 'results/distance_sensor.png', 20, 10);
    end

    % ---------------------------------------------------------------------
    % GPS Uncertainty
    % ---------------------------------------------------------------------
    if exist('gps_info', 'var')
        figure('Name', 'GPS Uncertainty', 'Color', 'w'); hold on;

        plot(gps_t, gps_info.eph, sty.axis1{:}, 'DisplayName', 'H Pos Accuracy (EPH) [m]');
        plot(gps_t, gps_info.epv, sty.axis2{:}, 'DisplayName', 'V Pos Accuracy (EPV) [m]');

        if ismember('hdop', gps_info.Properties.VariableNames)
            plot(gps_t, gps_info.hdop, sty.axis1_dash{:}, 'DisplayName', 'HDOP [m]');
            plot(gps_t, gps_info.vdop, sty.axis2_dash{:}, 'DisplayName', 'VDOP [m]');
        end

        plot(gps_t, gps_info.s_variance, sty.axis4_dot{:}, 'DisplayName', 'Speed Accuracy [m/s]');
        plot(gps_t, gps_info.satellites, sty.axis5_bold{:}, 'DisplayName', 'Satellites Used');
        plot(gps_t, gps_info.fix_type, sty.axis3_bold{:}, 'DisplayName', 'Fix Type');

        grid on;
        legend('show', 'Location', 'best', 'NumColumns', 2);
        ylabel('Value');
        title('GPS Uncertainty & Status');
        xlabel('Time (s)');
        ylim([0, 40]);
        add_standard_background(vis_flight_intervals, vis_flight_names, ...
            vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
        % PlotToFile(gcf, 'results/GPS_Uncertainty.png', 20, 10);
    end

    % ---------------------------------------------------------------------
    % GPS Noise & Jamming
    % ---------------------------------------------------------------------
    if exist('gps_info', 'var')
        figure('Name', 'GPS Noise & Jamming', 'Color', 'w'); hold on;
        plot(gps_t, gps_info.noise,   sty.axis2{:}, 'DisplayName', 'Noise per ms');
        plot(gps_t, gps_info.jamming, sty.axis1{:}, 'DisplayName', 'Jamming Indicator');

        grid on; legend('show');
        ylabel('Value');
        title('GPS Noise & Jamming');
        xlabel('Time (s)');
        add_standard_background(vis_flight_intervals, vis_flight_names, ...
            vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
        % PlotToFile(gcf, 'results/GPS_Noise_Jamming.png', 20, 10);
    end

    % ---------------------------------------------------------------------
    % Thrust and Magnetic Field
    % ---------------------------------------------------------------------
    figure('Name', 'Thrust and Magnetic Field', 'Color', 'w');
    mag_mag = sqrt(mag_data(:,1).^2 + mag_data(:,2).^2 + mag_data(:,3).^2);

    hold on;
    plot(mag_t, mag_mag, sty.axis1{:}, 'DisplayName', 'Magnetic Field Norm');

    if ~isempty(actuator_controls_0.thrust)
        plot(actuator_controls_0.time_thrust, actuator_controls_0.thrust, ...
            sty.thrust{:}, 'DisplayName', 'Thrust');
    end

    if log.data.vehicle_status_0.is_vtol(1) && ~dynamic_control_alloc && ~isempty(thrust_sp_1)
        plot(actuator_controls_1.time, actuator_controls_1.thrust_x, ...
            sty.thrust_alt{:}, 'DisplayName', 'Thrust (Fixed-wing)');
    end

    add_standard_background(vis_flight_intervals, vis_flight_names, ...
        vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

    grid on;
    xlabel('Time (s)');
    ylabel('Magnitude');
    title('Thrust and Magnetic Field');
    legend('show', 'Location', 'best');
    % PlotToFile(gcf, 'results/Thrust_Magnetic_Field.png', 20, 10);

    % ---------------------------------------------------------------------
    % Power (Battery & System)
    % ---------------------------------------------------------------------
    if exist('bat_v', 'var')
        figure('Name', 'Power', 'Color', 'w');

        ax1 = subplot(2,1,1); hold on;
        yyaxis left
        plot(bat_t, bat_v, sty.axis2_bold{:}, 'DisplayName', 'Voltage [V]');
        ylabel('Voltage [V]');

        yyaxis right
        plot(bat_t, bat_i, sty.axis1{:}, 'DisplayName', 'Current [A]');
        ylabel('Current [A]');

        title('Battery Voltage & Current');
        grid on;
        legend('show', 'Location', 'best');
        add_standard_background(vis_flight_intervals, vis_flight_names, ...
            vis_is_vtol, vis_vtol_intervals, vis_vtol_names);

        ax2 = subplot(2,1,2); hold on;
        yyaxis left
        plot(bat_t, bat_discharged/100, sty.axis4{:}, 'DisplayName', 'Discharged [mAh/100]');
        ylabel('Discharged [mAh/100]');

        yyaxis right
        plot(bat_t, bat_remaining*100, sty.axis3_bold{:}, 'DisplayName', 'Remaining [\%]');
        ylabel('Remaining [\%]');
        ylim([0, 105]);

        xlabel('Time (s)');
        grid on;
        legend('show', 'Location', 'best');
        linkaxes([ax1, ax2], 'x');
        add_standard_background(vis_flight_intervals, vis_flight_names, ...
            vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
        % PlotToFile(gcf, 'results/Battery_System.png', 20, 20);
    end

    % ---------------------------------------------------------------------
    % Temperature
    % ---------------------------------------------------------------------
    if exist('temp_data', 'var') && ~isempty(temp_data)
        figure('Name', 'Temperature', 'Color', 'w'); hold on;
        colors = make_channel_colors(length(temp_data), sty);

        for i = 1:length(temp_data)
            plot(temp_data(i).t, temp_data(i).val, ...
                'LineWidth', sty.lw_multi_bold, ...
                'Color', colors(i,:), ...
                'DisplayName', temp_data(i).name);
        end

        grid on; legend('show', 'Location', 'best');
        ylabel('Temperature [$^\circ$C]');
        title('System Temperatures');
        xlabel('Time (s)');
        add_standard_background(vis_flight_intervals, vis_flight_names, ...
            vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
        % PlotToFile(gcf, 'results/Temperature.png', 20, 10);
    end

    % ---------------------------------------------------------------------
    % Estimator Flags
    % ---------------------------------------------------------------------
    if isfield(log.data, 'estimator_status_0')
        figure('Name', 'Estimator Flags (Dynamic)', 'Color', 'w'); hold on;

        plot_count = 0;
        max_plots = 8;
        colors = make_channel_colors(max_plots, sty);
        active_legends = {};

        for i = 1:size(candidates, 1)
            lbl = candidates{i, 1};
            data = candidates{i, 2};

            if max(data) > 0.1
                plot_count = plot_count + 1;
                plot(est_t, data, 'Color', colors(plot_count,:), 'LineWidth', sty.lw_multi_bold);
                active_legends{end+1} = lbl;

                if plot_count >= max_plots
                    break;
                end
            end
        end

        if plot_count == 0
            plot(est_t, candidates{1,2}, sty.axis4{:});
            active_legends{end+1} = candidates{1,1};
            title('Estimator Flags (All Good)');
        else
            title(sprintf('Estimator Flags (Top %d Active)', plot_count));
        end

        grid on;
        legend(active_legends, 'Location', 'best');
        ylabel('Flag Value');
        xlabel('Time (s)');
        add_standard_background(vis_flight_intervals, vis_flight_names, ...
            vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
        % PlotToFile(gcf, 'results/Estimator_Flags.png', 20, 10);
    end

    % ---------------------------------------------------------------------
    % Failsafe Flags
    % ---------------------------------------------------------------------
    if isfield(log.data, 'failsafe_flags_0')
        figure('Name', 'Failsafe Flags', 'Color', 'w'); hold on;

        if exist('vs_failsafe', 'var')
            area(vs_t, double(vs_failsafe), 'FaceColor', sty.fill_warn, ...
                'EdgeColor', 'none', 'DisplayName', 'In Failsafe Mode');
        end

        plot_idx = 0;
        colors = make_channel_colors(10, sty);

        for i = 1:length(fs_cols)
            col_name = fs_cols{i};
            if strcmp(col_name, 'timestamp') || startsWith(col_name, 'mode_req')
                continue;
            end

            data = fs_table.(col_name);
            if max(data) > 0
                plot_idx = plot_idx + 1;
                plot(fs_t, double(data) * 0.8 + plot_idx, ...
                    'LineWidth', sty.lw_multi_bold, ...
                    'Color', colors(mod(plot_idx-1,10)+1,:), ...
                    'DisplayName', strrep(col_name, '_', ' '));
            end
        end

        if plot_idx == 0
            text(mean(fs_t), 0.5, 'No Failsafe Flags Triggered', 'HorizontalAlignment', 'center');
        end

        grid on; legend('show');
        title('Failsafe Flags Triggered');
        ylabel('Flags (Stacked)');
        xlabel('Time (s)');
        ylim([0, plot_idx + 2]);
        add_standard_background(vis_flight_intervals, vis_flight_names, ...
            vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
        % PlotToFile(gcf, 'results/Failsafe_Flags.png', 20, 10);
    end

    % ---------------------------------------------------------------------
    % CPU & RAM
    % ---------------------------------------------------------------------
    if isfield(log.data, 'cpuload_0')
        figure('Name', 'CPU & RAM', 'Color', 'w'); hold on;
        plot(cpu_t, ram_usage, sty.axis2_bold{:}, 'DisplayName', 'RAM Usage');
        plot(cpu_t, cpu_load,  sty.axis1_bold{:}, 'DisplayName', 'CPU Load');

        grid on;
        legend('show', 'Location', 'best');
        ylabel('Load / Usage [0-1]');
        title('CPU & RAM');
        xlabel('Time (s)');
        ylim([0, 1]);
        add_standard_background(vis_flight_intervals, vis_flight_names, ...
            vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
        % PlotToFile(gcf, 'results/CPU.png', 20, 10);
    end

    % ---------------------------------------------------------------------
    % Sampling Regularity
    % ---------------------------------------------------------------------
    if isfield(log.data, 'sensor_combined_0')
        figure('Name', 'Sampling Regularity', 'Color', 'w'); hold on;

        sc = log.data.sensor_combined_0;
        t_raw = double(sc.timestamp);
        dt_seq = diff(t_raw);
        t_plot_sc = t_raw(2:end) * 1e-6;

        plot(t_plot_sc, dt_seq, sty.axis3_thin{:}, ...
            'DisplayName', 'delta t (between 2 logged samples)');

        if isfield(log.data, 'estimator_status_0')
            es = log.data.estimator_status_0;
            t_es = double(es.timestamp) * 1e-6;
            slip_us = double(es.time_slip) * 1e6;

            plot(t_es, slip_us, sty.axis4_bold{:}, ...
                'DisplayName', 'Estimator time slip (cumulative)');
        end

        grid on;
        legend('show', 'Location', 'best');
        ylabel('[us]');
        title('Sampling Regularity of Sensor Data');
        xlabel('Time (s)');
        
        % Python limits Y axis range to [0, 25000] us (i.e. 0-25ms)
        % This is very effective for observing 250Hz (4000us) or 1kHz (1000us) data
        ylim([0, 25000]); 
        
        % Add background
        add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names);
        % PlotToFile(gcf, 'results/Sampling_Regularity.png', 20, 10);
    end
end

%% =========================================================================
%  Local functions
% =========================================================================
function sty = make_plot_style()
    % -------- Core pair used in Fig.1-5 --------
    sty.c_setpoint = [0.12, 0.12, 0.12];
    sty.c_response = [0.70, 0.22, 0.40];

    % -------- Multi-signal palette --------
    sty.c_axis1 = [0.00, 0.447, 0.741];
    sty.c_axis2 = [0.850, 0.325, 0.098];
    sty.c_axis3 = [0.466, 0.674, 0.188];
    sty.c_axis4 = [0.494, 0.184, 0.556];
    sty.c_axis5 = [0.929, 0.694, 0.125];
    sty.c_axis6 = [0.301, 0.745, 0.933];
    sty.c_axis7 = [0.635, 0.078, 0.184];

    % -------- Line widths --------
    sty.lw_sp = 0.8;
    sty.lw_res = 0.5;
    sty.lw_main = 0.8;
    sty.lw_main_bold = 1;
    sty.lw_thin = 0.5;
    sty.lw_multi = 0.9;
    sty.lw_multi_bold = 1.2;

    % -------- Frequently used styles --------
    sty.setpoint = {'Color', sty.c_setpoint, 'LineStyle', '--', 'LineWidth', sty.lw_sp};
    sty.response = {'Color', sty.c_response, 'LineStyle', '-',  'LineWidth', sty.lw_res};
    sty.response_fade = {'Color', [0.70, 0.22, 0.40, 0.45], 'LineStyle', '-.', 'LineWidth', sty.lw_main};

    sty.axis1 = {'Color', sty.c_axis1, 'LineStyle', '-', 'LineWidth', sty.lw_main};
    sty.axis2 = {'Color', sty.c_axis2, 'LineStyle', '--', 'LineWidth', sty.lw_main};
    sty.axis3 = {'Color', sty.c_axis3, 'LineStyle', '-.', 'LineWidth', sty.lw_main};
    sty.axis4 = {'Color', sty.c_axis4, 'LineStyle', ':', 'LineWidth', sty.lw_main};
    sty.axis5 = {'Color', sty.c_axis5, 'LineStyle', '-', 'LineWidth', sty.lw_main};

    sty.axis1_bold = {'Color', sty.c_axis1, 'LineStyle', '-', 'LineWidth', sty.lw_main_bold};
    sty.axis2_bold = {'Color', sty.c_axis2, 'LineStyle', '--', 'LineWidth', sty.lw_main_bold};
    sty.axis3_bold = {'Color', sty.c_axis3, 'LineStyle', '-.', 'LineWidth', sty.lw_main_bold};
    sty.axis4_bold = {'Color', sty.c_axis4, 'LineStyle', ':', 'LineWidth', sty.lw_main_bold};
    sty.axis5_bold = {'Color', sty.c_axis5, 'LineStyle', '-', 'LineWidth', sty.lw_main_bold};

    sty.axis1_thin = {'Color', sty.c_axis1, 'LineStyle', '-', 'LineWidth', sty.lw_thin};
    sty.axis2_thin = {'Color', sty.c_axis2, 'LineStyle', '--', 'LineWidth', sty.lw_thin};
    sty.axis3_thin = {'Color', sty.c_axis3, 'LineStyle', '-.', 'LineWidth', sty.lw_thin};

    sty.axis1_dash = {'Color', sty.c_axis1, 'LineStyle', '--', 'LineWidth', sty.lw_thin};
    sty.axis2_dash = {'Color', sty.c_axis2, 'LineStyle', '-.', 'LineWidth', sty.lw_thin};
    sty.axis4_dot  = {'Color', sty.c_axis4, 'LineStyle', ':',  'LineWidth', sty.lw_main};

    sty.aux1 = {'Color', sty.c_axis4, 'LineStyle', '--', 'LineWidth', sty.lw_thin};

    sty.thrust = {'Color', [0.10, 0.10, 0.10], 'LineStyle', '-',  'LineWidth', sty.lw_main_bold};
    sty.thrust_alt = {'Color', [0.10, 0.10, 0.10], 'LineStyle', '--', 'LineWidth', sty.lw_main};
    sty.thrust_bold = {'Color', [0.10, 0.10, 0.10], 'LineStyle', '-', 'LineWidth', 1.4};

    sty.fill_warn = [1.00, 0.85, 0.85];
end

function cols = make_channel_colors(n, sty)
    base = [sty.c_axis1;
        sty.c_axis2;
        sty.c_axis3;
        sty.c_axis4;
        sty.c_axis5;
        sty.c_axis6;
        sty.c_axis7];

    if n <= size(base,1)
        cols = base(1:n, :);
    else
        cols = lines(n);
    end
end


function idx_out = subsample_idx(idx_in, step)
    if isempty(idx_in)
        idx_out = idx_in;
        return;
    end
    idx_out = idx_in(1:step:end);
    if idx_out(end) ~= idx_in(end)
        idx_out = [idx_out; idx_in(end)];
    end
end