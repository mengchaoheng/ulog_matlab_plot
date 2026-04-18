close all;
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
    PlotToFile(gcf, 'results/traj.png', 12, 4);
end