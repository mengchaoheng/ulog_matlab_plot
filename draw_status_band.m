%% =========================================================================
%  Helper Function 3: Generic Status Band Drawing (Support specifying Y axis range)
% =========================================================================
function draw_status_band(intervals, labels, y_range, color_mode)
    % y_range: [y_bottom, y_top] specify drawing area
    % color_mode: 'flight_mode' or 'vtol_state' to select color scheme
    
    if isempty(intervals), return; end
    
    % --- Color Scheme ---
    colors = containers.Map('KeyType', 'double', 'ValueType', 'any');
    default_c = [0.95, 0.95, 0.95]; 
    
    if strcmp(color_mode, 'flight_mode')
       % Flight mode colors (soft background colors)
        colors(0)  = [0.88, 0.92, 0.98]; % Manual (= Stabilized)
        colors(15) = [0.88, 0.92, 0.98]; % Stabilized
        colors(1)  = [0.78, 0.98, 0.98]; % Altitude
        colors(2)  = [0.78, 0.90, 1.00]; % Position
        colors(10) = [0.90, 0.90, 0.80]; % Acro 
        
        colors(17) = [1.00, 0.96, 0.70]; % Takeoff 
        colors(22) = [1.00, 0.96, 0.70]; % VTOL Takeoff
        colors(18) = [0.78, 0.98, 0.82]; % Land 
        colors(20) = [0.78, 0.98, 0.82]; % Precland
        colors(5)  = [1.00, 0.88, 0.70]; % RTL
        
        colors(3)  = [0.88, 0.85, 1.00]; % Mission
        colors(4)  = [0.92, 0.85, 0.92]; % Loiter
        colors(14) = [0.98, 0.82, 0.98]; % Offboard
        
        alpha_val = 0.2; 
        
    elseif strcmp(color_mode, 'vtol_state')
       % VTOL state colors (more vivid, for bottom band)
        colors(1) = [0.65, 0.75, 0.90]; % MC: 
        colors(2) = [0.65, 0.85, 0.65]; % FW: 
        colors(3) = [1.00, 0.80, 0.50]; % Transition: 
        
        alpha_val = 0.80; 
    end
    
    hold on;
    y_b = y_range(1); % bottom
    y_t = y_range(2); % top
    
    % === [Modification Start] Text Position Calculation ===
    % Determine text vertical position and alignment based on mode
    if strcmp(color_mode, 'flight_mode')
        % Flight mode: place at top (offset down 2% of height from top to avoid line overlap)
        text_y = y_t - (y_t - y_b) * 0.05;
        text_valign = 'top';   % Vertical alignment: top
    else
        % VTOL or other modes: vertically centered (since VTOL band is already at bottom)
        text_y = y_b + (y_t - y_b) / 2;
        text_valign = 'middle'; % Vertical alignment: middle
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
        
        % Draw rectangle
        p = patch([t_s t_e t_e t_s], [y_b y_b y_t y_t], c);
        set(p, 'EdgeColor', 'none', 'FaceAlpha', alpha_val, 'HandleVisibility', 'off');
        
        % Text label 
        % Only display text when duration is long enough (>2s)
        if (t_e - t_s) > 2.0
            text(t_s + (t_e-t_s)/2, text_y, labels{i}, ... 
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', text_valign, ...      
                'FontSize', 7, 'FontName', 'Times New Roman', ... %  
                'Color', [0.3 0.3 0.3], ... %  
                'Interpreter', 'none', 'Clipping', 'on');
        end
    end
    % Ensure grid and curves are on top layer
    set(gca, 'Layer', 'top');
end