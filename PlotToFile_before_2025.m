function PlotToFile(fig, filename, width_cm, height_cm)

    if nargin < 4
        error('PlotToFile requires 4 inputs: figure_handle, filename, width_cm, height_cm.');
    end

    if ~ishandle(fig) || ~strcmp(get(fig, 'Type'), 'figure')
        error('The first input is not a valid figure handle.');
    end

    out_dir = fileparts(filename);
    if ~isempty(out_dir) && ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    % --- figure size first ---
    set(fig, 'Units', 'centimeters');
    set(fig, 'Position', [2, 2, width_cm, height_cm]);
    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperPosition', [0, 0, width_cm, height_cm]);
    set(fig, 'PaperSize', [width_cm, height_cm]);
    set(fig, 'Color', 'w');

    % --- flush pending graphics update ---
    drawnow;
    pause(0.05);

    % --- only process existing legends ---
    lgds = findall(fig, 'Type', 'legend');
    for k = 1:numel(lgds)
        lgd = lgds(k);

        % only re-optimize legends that are auto/best managed
        try
            if strcmpi(lgd.Location, 'best')
                lgd.Location = 'best';
                drawnow;
                lgd.Units = 'normalized';
                pos = lgd.Position;
                lgd.Location = 'none';
                lgd.Position = pos;
            end
        catch
        end
    end

    drawnow;
    pause(0.05);

    [~, ~, ext] = fileparts(filename);
    ext = lower(ext);

    % R2025a
    if ~verLessThan('matlab', '25.1')
        exportgraphics(fig, filename, ...
            'ContentType', 'vector', ...
            'BackgroundColor', 'white', ...
            'Resolution', 1200, ...
            'Width',  width_cm, ...
            'Height', height_cm, ...
            'Padding', 'tight', ...
            'Units', 'centimeters');
        return;
    end

    switch ext
        case {'.png', '.jpg', '.jpeg', '.tif', '.tiff'}
            exportgraphics(fig, filename, ...
                'BackgroundColor', 'white', ...
                'Resolution', 1200);

        case '.pdf'
            exportgraphics(fig, filename, ...
                'ContentType', 'vector', ...
                'BackgroundColor', 'white');

        otherwise
            error('Unsupported file extension: %s', ext);
    end
end