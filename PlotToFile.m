function PlotToFile(gcf, filename, width_cm, height_cm)

    if nargin < 4
        error('PlotToFile requires 4 inputs: figure_handle, filename, width_cm, height_cm.');
    end

    if ~ishandle(gcf) || ~strcmp(get(gcf, 'Type'), 'figure')
        error('The first input is not a valid figure handle.');
    end

    out_dir = fileparts(filename);
    if ~isempty(out_dir) && ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    % R2025a 
    if ~verLessThan('matlab', '25.1')
        set(gcf, 'Units', 'centimeters');
        set(gcf, 'Position', [0, 0, width_cm, height_cm]);

        exportgraphics(gcf, filename, ...
            'ContentType', 'vector', ...
            'BackgroundColor', 'white', ...
            'Resolution', 1200, ...
            'Width',  width_cm, ...
            'Height', height_cm, ...
            'Padding', "tight", ...
            'Units', "centimeters");
        return;
    end

    % Old versions
    set(gcf, 'Units', 'centimeters');
    set(gcf, 'Position', [2, 2, width_cm, height_cm]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0, 0, width_cm, height_cm]);
    set(gcf, 'PaperSize', [width_cm, height_cm]);
    set(gcf, 'Color', 'w');

    [~, ~, ext] = fileparts(filename);
    ext = lower(ext);

    switch ext
        case {'.png', '.jpg', '.jpeg', '.tif', '.tiff'}
            exportgraphics(gcf, filename, ...
                'BackgroundColor', 'white', ...
                'Resolution', 600);

        case '.pdf'
            exportgraphics(gcf, filename, ...
                'ContentType', 'vector', ...
                'BackgroundColor', 'white');

        otherwise
            error('Unsupported file extension: %s', ext);
    end
end