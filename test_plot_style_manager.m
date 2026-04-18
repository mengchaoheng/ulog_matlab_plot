function test_plot_style_manager()
    sty = local_make_plot_style();
    t = linspace(0, 10, 400);

    figure('Name', 'Style Manager Test', 'Color', 'w');

    subplot(2,1,1); hold on;
    plot(t, sin(t),             sty.setpoint{:},  'DisplayName', 'sty.setpoint');
    plot(t, sin(t + 0.15),      sty.response{:},  'DisplayName', 'sty.response');
    plot(t, 0.8*cos(t),         sty.reference{:}, 'DisplayName', 'sty.reference');
    plot(t, 0.8*cos(t + 0.2),   sty.estimate{:},  'DisplayName', 'sty.estimate');
    plot(t, 0.5*sin(2*t),       sty.raw{:},       'DisplayName', 'sty.raw');
    plot(t, 0.5*cos(2*t),       sty.command{:},   'DisplayName', 'sty.command');
    plot(t, 0.3*sin(3*t),       sty.state{:},     'DisplayName', 'sty.state');
    plot(t, 0.3*cos(3*t),       sty.auxiliary{:}, 'DisplayName', 'sty.auxiliary');
    plot(t, 0.2*sin(4*t),       sty.event{:},     'DisplayName', 'sty.event');
    grid on; xlabel('Time (s)'); ylabel('Amplitude');
    title('Semantic foreground styles');
    legend('show', 'Location', 'eastoutside');

    subplot(2,1,2); hold on;
    plot(t, 1.00 + 0*t, sty.axis1{:},      'DisplayName', 'sty.axis1');
    plot(t, 0.85 + 0*t, sty.axis2{:},      'DisplayName', 'sty.axis2');
    plot(t, 0.70 + 0*t, sty.axis3{:},      'DisplayName', 'sty.axis3');
    plot(t, 0.55 + 0*t, sty.axis4{:},      'DisplayName', 'sty.axis4');
    plot(t, 0.40 + 0*t, sty.axis5{:},      'DisplayName', 'sty.axis5');
    plot(t, 0.25 + 0*t, sty.axis1_bold{:}, 'DisplayName', 'sty.axis1_bold');
    plot(t, 0.10 + 0*t, sty.axis3_thin{:}, 'DisplayName', 'sty.axis3_thin');
    plot(t,-0.05 + 0*t, sty.thrust{:},     'DisplayName', 'sty.thrust');
    plot(t,-0.20 + 0*t, sty.thrust_alt{:}, 'DisplayName', 'sty.thrust_alt');
    plot(t,-0.35 + 0*t, sty.thrust_bold{:},'DisplayName', 'sty.thrust_bold');
    ylim([-0.5, 1.1]);
    grid on; xlabel('Time (s)'); ylabel('Style preview');
    title('Legacy compatibility styles');
    legend('show', 'Location', 'eastoutside');
end

function sty = local_make_plot_style()
    sty.c_setpoint = [0.12, 0.12, 0.12];
    sty.c_response = [0.70, 0.22, 0.40];
    sty.c_ref      = [0.00, 0.447, 0.741];
    sty.c_est      = [0.850, 0.325, 0.098];
    sty.c_raw      = [0.466, 0.674, 0.188];
    sty.c_cmd      = [0.494, 0.184, 0.556];
    sty.c_state    = [0.929, 0.694, 0.125];
    sty.c_aux      = [0.301, 0.745, 0.933];
    sty.c_event    = [0.635, 0.078, 0.184];
    sty.c_black    = [0.10, 0.10, 0.10];

    sty.lw_xs = 0.5;
    sty.lw_s  = 0.7;
    sty.lw_m  = 0.9;
    sty.lw_l  = 1.1;
    sty.lw_xl = 1.4;

    sty.lw_sp         = sty.lw_m;
    sty.lw_res        = sty.lw_xs;
    sty.lw_main       = sty.lw_m;
    sty.lw_main_bold  = sty.lw_l;
    sty.lw_thin       = sty.lw_xs;
    sty.lw_multi      = sty.lw_m;
    sty.lw_multi_bold = sty.lw_xl;

    sty.setpoint  = {'Color', sty.c_setpoint, 'LineStyle', '--', 'LineWidth', sty.lw_sp};
    sty.response  = {'Color', sty.c_response, 'LineStyle', '-',  'LineWidth', sty.lw_res};
    sty.reference = {'Color', sty.c_ref,      'LineStyle', '--', 'LineWidth', sty.lw_s};
    sty.estimate  = {'Color', sty.c_est,      'LineStyle', '-',  'LineWidth', sty.lw_m};
    sty.raw       = {'Color', sty.c_raw,      'LineStyle', ':',  'LineWidth', sty.lw_xs};
    sty.command   = {'Color', sty.c_cmd,      'LineStyle', '-.', 'LineWidth', sty.lw_m};
    sty.state     = {'Color', sty.c_state,    'LineStyle', '-',  'LineWidth', sty.lw_l};
    sty.auxiliary = {'Color', sty.c_aux,      'LineStyle', '--', 'LineWidth', sty.lw_xs};
    sty.event     = {'Color', sty.c_event,    'LineStyle', ':',  'LineWidth', sty.lw_l};

    sty.c_axis1 = sty.c_ref;
    sty.c_axis2 = sty.c_est;
    sty.c_axis3 = sty.c_raw;
    sty.c_axis4 = sty.c_cmd;
    sty.c_axis5 = sty.c_state;

    sty.axis1 = {'Color', sty.c_axis1, 'LineStyle', '-',  'LineWidth', sty.lw_main};
    sty.axis2 = {'Color', sty.c_axis2, 'LineStyle', '--', 'LineWidth', sty.lw_main};
    sty.axis3 = {'Color', sty.c_axis3, 'LineStyle', '-.', 'LineWidth', sty.lw_main};
    sty.axis4 = {'Color', sty.c_axis4, 'LineStyle', ':',  'LineWidth', sty.lw_main};
    sty.axis5 = {'Color', sty.c_axis5, 'LineStyle', '-',  'LineWidth', sty.lw_main};
    sty.axis1_bold = {'Color', sty.c_axis1, 'LineStyle', '-',  'LineWidth', sty.lw_main_bold};
    sty.axis3_thin = {'Color', sty.c_axis3, 'LineStyle', '-.', 'LineWidth', sty.lw_thin};
    sty.thrust      = {'Color', sty.c_black, 'LineStyle', '-',  'LineWidth', sty.lw_main_bold};
    sty.thrust_alt  = {'Color', sty.c_black, 'LineStyle', '--', 'LineWidth', sty.lw_main};
    sty.thrust_bold = {'Color', sty.c_black, 'LineStyle', '-',  'LineWidth', sty.lw_xl};
end
