function add_standard_background(vis_flight_intervals, vis_flight_names, vis_is_vtol, vis_vtol_intervals, vis_vtol_names)
    ax = gca;

    % 1. 当前坐标范围
    yl = ylim(ax);
    y_h = yl(2) - yl(1);

    % 2. 画背景带
    draw_status_band(vis_flight_intervals, vis_flight_names, yl, 'flight_mode');

    if vis_is_vtol
        y_band_top = yl(1) + 0.1 * y_h;
        draw_status_band(vis_vtol_intervals, vis_vtol_names, [yl(1), y_band_top], 'vtol_state');
    end

    % 3. 将背景对象压到底层：按"所有 children 的合法排列"重排
    ch = ax.Children;
    if isempty(ch)
        set(ax, 'Layer', 'top');
        return;
    end

    tags = get(ch, 'Tag');
    if ischar(tags)
        tags = {tags};
    end

    is_bg = strcmp(tags, 'status_bg');
    ax.Children = [ch(~is_bg); ch(is_bg)];

    % 4. 坐标轴刻度和网格保持在上层
    set(ax, 'Layer', 'top');
end