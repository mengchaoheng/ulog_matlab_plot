function style_cell = plot_style(sty, color_name, line_style, line_width)
%PLOT_STYLE Build a line style cell from named color, line style, and line width.

    if ~isfield(sty, 'color') || ~isfield(sty.color, color_name)
        error('plot_style:UnknownColor', 'Unknown color name: %s', color_name);
    end

    style_cell = {'Color', sty.color.(color_name), ...
                  'LineStyle', line_style, ...
                  'LineWidth', line_width};
end
