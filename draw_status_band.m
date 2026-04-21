%% =========================================================================
%  Helper Function 3: Generic Status Band Drawing (Support specifying Y axis range)
% =========================================================================
function draw_status_band(intervals, labels, y_range, color_mode)
    % y_range: [y_bottom, y_top] specify drawing area
    % color_mode: 'flight_mode' or 'vtol_state' to select color scheme

    if isempty(intervals), return; end

    sty = plot_style_manager();
    default_c = sty.color.bg_default;

    switch color_mode
        case 'flight_mode'
            cfg = sty.bg.flight_mode;
        case 'vtol_state'
            cfg = sty.bg.vtol_state;
        otherwise
            cfg = struct();
            cfg.colors = containers.Map('KeyType', 'double', 'ValueType', 'any');
            cfg.alpha = 0.2;
    end

    colors = cfg.colors;
    alpha_val = cfg.alpha;

    hold on;
    y_b = y_range(1);
    y_t = y_range(2);

    if strcmp(color_mode, 'flight_mode')
        text_y = y_t - (y_t - y_b) * 0.05;
        text_valign = 'top';
    else
        text_y = y_b + (y_t - y_b) / 2;
        text_valign = 'middle';
    end

    for i = 1:size(intervals, 1)
        t_s = intervals(i, 1);
        t_e = intervals(i, 2);
        val = intervals(i, 3);

        if isKey(colors, val)
            c = colors(val);
        else
            c = default_c;
        end

        p = patch([t_s t_e t_e t_s], [y_b y_b y_t y_t], c);
        set(p, 'EdgeColor', 'none', 'FaceAlpha', alpha_val, 'HandleVisibility', 'off');

        if (t_e - t_s) > 0.5
            text(t_s + (t_e - t_s) / 2, text_y, labels{i}, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', text_valign, ...
                'FontSize', 5, 'FontName', 'Times New Roman', ...
                'Color', [0.2 0.2 0.2], ...
                'Interpreter', 'none', 'Clipping', 'on');
        end
    end

    set(gca, 'Layer', 'top');
end
