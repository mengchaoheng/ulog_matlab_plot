function sty = plot_style_manager()
%PLOT_STYLE_MANAGER Centralized named color definition.
%   This file only manages named colors. Line style and line width are
%   provided at the call site through plot_style(sty, color_name, line_style, line_width).
%   flight_mode and vtol_state colors are kept unchanged.

    sty.color = struct();

    % ==================== neutral / utility colors ====================
    sty.color.neutral1  = [0.10, 0.10, 0.10];
    sty.color.neutral2  = [0.35, 0.35, 0.35];
    sty.color.neutral3  = [0.70, 0.70, 0.70];
    sty.color.warn_fill = [1.00, 0.85, 0.85];
    sty.color.bg_default = [0.95, 0.95, 0.95];

    % ==================== setpoint family ====================
    % Rich cool-color family. Organized as hue pairs with dark/light levels.
    % blue
    sty.color.setpoint1  = [0.08, 0.22, 0.48];
    sty.color.setpoint2  = [0.36, 0.58, 0.88];
    % cyan
    sty.color.setpoint3  = [0.00, 0.42, 0.52];
    sty.color.setpoint4  = [0.36, 0.76, 0.86];
    % teal
    sty.color.setpoint5  = [0.00, 0.46, 0.40];
    sty.color.setpoint6  = [0.36, 0.76, 0.66];
    % green
    sty.color.setpoint7  = [0.18, 0.50, 0.22];
    sty.color.setpoint8  = [0.56, 0.78, 0.42];
    % indigo / violet
    sty.color.setpoint9  = [0.32, 0.30, 0.62];
    sty.color.setpoint10 = [0.68, 0.62, 0.86];

    % ==================== response family ====================
    % Rich warm-color family. Organized as hue pairs with dark/light levels.
    % red
    sty.color.response1  = [0.62, 0.12, 0.16];
    sty.color.response2  = [0.90, 0.42, 0.40];
    % orange
    sty.color.response3  = [0.78, 0.34, 0.08];
    sty.color.response4  = [0.96, 0.62, 0.26];
    % amber / gold
    sty.color.response5  = [0.66, 0.46, 0.06];
    sty.color.response6  = [0.90, 0.74, 0.22];
    % magenta / rose
    sty.color.response7  = [0.62, 0.18, 0.42];
    sty.color.response8  = [0.88, 0.54, 0.72];
    % brown / sand
    sty.color.response9  = [0.52, 0.30, 0.12];
    sty.color.response10 = [0.82, 0.66, 0.46];

% ==================== background colors ====================
    sty.color.flight_mode_0  = [0.88, 0.92, 0.98];
    sty.color.flight_mode_15 = [0.88, 0.92, 0.98];
    sty.color.flight_mode_1  = [0.78, 0.98, 0.98];
    sty.color.flight_mode_2  = [0.78, 0.90, 1.00];
    sty.color.flight_mode_10 = [0.90, 0.90, 0.80];
    sty.color.flight_mode_17 = [1.00, 0.96, 0.70];
    sty.color.flight_mode_22 = [1.00, 0.96, 0.70];
    sty.color.flight_mode_18 = [0.78, 0.98, 0.82];
    sty.color.flight_mode_20 = [0.78, 0.98, 0.82];
    sty.color.flight_mode_5  = [1.00, 0.88, 0.70];
    sty.color.flight_mode_3  = [0.88, 0.85, 1.00];
    sty.color.flight_mode_4  = [0.92, 0.85, 0.92];
    sty.color.flight_mode_14 = [0.98, 0.82, 0.98];

    sty.bg.flight_mode.colors = containers.Map('KeyType', 'double', 'ValueType', 'any');
    sty.bg.flight_mode.colors(0)  = sty.color.flight_mode_0;
    sty.bg.flight_mode.colors(15) = sty.color.flight_mode_15;
    sty.bg.flight_mode.colors(1)  = sty.color.flight_mode_1;
    sty.bg.flight_mode.colors(2)  = sty.color.flight_mode_2;
    sty.bg.flight_mode.colors(10) = sty.color.flight_mode_10;
    sty.bg.flight_mode.colors(17) = sty.color.flight_mode_17;
    sty.bg.flight_mode.colors(22) = sty.color.flight_mode_22;
    sty.bg.flight_mode.colors(18) = sty.color.flight_mode_18;
    sty.bg.flight_mode.colors(20) = sty.color.flight_mode_20;
    sty.bg.flight_mode.colors(5)  = sty.color.flight_mode_5;
    sty.bg.flight_mode.colors(3)  = sty.color.flight_mode_3;
    sty.bg.flight_mode.colors(4)  = sty.color.flight_mode_4;
    sty.bg.flight_mode.colors(14) = sty.color.flight_mode_14;
    sty.bg.flight_mode.alpha = 0.2;

    sty.color.vtol_state_1 = [0.65, 0.75, 0.90];
    sty.color.vtol_state_2 = [0.65, 0.85, 0.65];
    sty.color.vtol_state_3 = [1.00, 0.80, 0.50];

    sty.bg.vtol_state.colors = containers.Map('KeyType', 'double', 'ValueType', 'any');
    sty.bg.vtol_state.colors(1) = sty.color.vtol_state_1;
    sty.bg.vtol_state.colors(2) = sty.color.vtol_state_2;
    sty.bg.vtol_state.colors(3) = sty.color.vtol_state_3;
    sty.bg.vtol_state.alpha = 0.80;
end
