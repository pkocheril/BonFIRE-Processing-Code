% Automatic gradient mapping and marker selection
function [LINECOLOR,MARKER,MARKSIZE] = gradientset(INDEX, LOOPLENGTH, COLORMAP)
%%% This function automatically chooses colors and marker sets for plotting
% data in a loop based on sets of pre-defined colors and markers.
% Updated version: generalized to any number of color stops.

% First color in list will be at the top of the legend
    % Previous colormaps
    presets.rgb = [1 0 0; 0 1 0; 0 0 1]; % red → green → blue
    presets.fire  = [255 134 29; 200 60 60; 100 0 100]./255; % orange → red → purple
    presets.rwb1 = [231 48 41; 235 242 250; 48 124 188]./255; % red → white → blue
    presets.rwb2 = [103 0 30; 240 240 240; 6 48 97]./255; % red → white → blue (alternate)
    presets.bwr = [48 124 188; 235 242 250; 231 48 41]./255; % blue → white → red
    presets.bwr2 = [6 48 97; 240 240 240; 103 0 30]./255; % blue → white → red (alternate)
    presets.rpb = [200 0 30; 100 0 170; 0 60 200]./255; % red → purple → blue
    presets.bpr = [0 60 200; 100 0 170; 200 0 30]./255; % blue → purple → red
    presets.parula = [253 231 36; 38 138 141; 68 1 84]./255; % yellow → green → blue
    presets.parula2 = [240 249 185; 65 182 196; 6 29 88]./255; % yellow → green → blue (alternate)
    presets.dfg = [0.5 0 0.5; 1 1 1; 0.2 0.8 0.2]; % purple → white → green
    presets.idler = [0 0.35 0.5; 1 1 1; 0.8 0.2 0.2]; % turquoise → white → red
    presets.probe = [0 0 1; 1 1 1; 1 0.5 0]; % blue → white → orange

    % Inverted viridis colormaps
    presets.fireinv = [100 0 100; 200 60 60; 255 134 29]./255; % purple → red → orange
    presets.parulainv = [68 1 84; 38 138 141; 253 231 36]./255; % blue → green → yellow
    presets.dfginv = [0.2 0.8 0.2; 1 1 1; 0.5 0 0.5]; % green → white → purple
    presets.idler = [0.8 0.2 0.2; 1 1 1; 0 0.35 0.5]; % red → white → turquoise
    presets.probe = [1 0.5 0; 1 1 1; 0 0 1]; % orange → white → blue
    presets.unionjackinv = [0 0 0; 0 0 1; 0.8 0.8 0.8; 1 0 0; 0.9 0.9 0.9];
    presets.viridisinv = [52 0 66; 54 34 107; 42 72 122; 33 108 123; 31 146 116; 63 184 90; 144 214 45; 252 229 30]./255;
    presets.plasmainv = [11 0 116; 65 0 146; 118 0 148; 168 25 118; 208 69 85; 239 115 57; 252 175 33; 237 252 27]./255;
    presets.magmainv = [0 0 5; 26 10 64; 75 0 108; 132 24 109; 198 43 91; 243 95 74; 252 172 109; 251 255 178]./255;
    presets.infernoinv = [0 0 5; 29 0 66; 81 0 91; 140 23 80; 200 50 51; 240 104 19; 247 182 31; 252 255 147]./255;
    presets.cividisinv = [1 23 60; 18 41 90; 59 66 89; 89 91 95; 123 120 102; 164 152 95; 210 189 79; 255 232 55]./255;
    presets.rocketinv = [4 5 20; 40 17 47; 92 17 69; 156 0 71; 215 27 51; 237 95 64; 242 164 124; 248 230 213]./255;
    presets.makoinv = [11 5 6; 34 21 46; 50 43 104; 43 80 140; 43 125 150; 54 170 157; 122 211 162; 215 244 223]./255;
    presets.turboinv = [36 12 45; 56 93 235; 34 199 203; 87 255 88; 201 232 41; 251 137 35; 208 37 11; 101 0 5]./255;

    % Viridis colormaps
    presets.unionjack = flipud(presets.unionjackinv);
    presets.viridis = flipud(presets.viridisinv);
    presets.plasma = flipud(presets.plasmainv);
    presets.magma = flipud(presets.magmainv);
    presets.inferno = flipud(presets.infernoinv);
    presets.cividis = flipud(presets.cividisinv);
    presets.rocket = flipud(presets.rocketinv);
    presets.mako = flipud(presets.makoinv);
    presets.turbo = flipud(presets.turboinv);

    % New colormaps
    presets.bwr = [2 24 147; 35 112 205; 14 152 227; 134 216 255; 177 230 250; 247 249 249; 255 243 241; 253 156 117; 255 126 121; 255 37 1; 179 37 1]./255;
    presets.rwb = flipud(presets.bwr);
    presets.rainbowinv = [134 112 177; 79 126 188; 79 171 200; 104 203 184; 155 211 144; 201 219 87; 248 189 84; 237 76 82]./255;
    presets.rainbow = flipud(presets.rainbowinv);
    presets.rainbow2 = [197 42 48; 247 151 32; 249 230 47; 150 201 61; 73 186 97; 43 146 154; 42 108 177; 37 43 107]./255;
    presets.rainbow2inv = flipud(presets.rainbow2);

    % If COLORMAP is missing or 'default'
    if nargin < 3 || isempty(COLORMAP) || strcmpi(COLORMAP, 'default')
        COLORMAP = 'turbo';
    end

    % If COLORMAP is a string and matches a preset
    if ischar(COLORMAP) || isstring(COLORMAP)
        cmapName = lower(string(COLORMAP));
        if isfield(presets, cmapName)
            colorStops = presets.(cmapName);
        else
            warning('Unknown colormap "%s". Using turbo.', cmapName);
            colorStops = presets.turbo;
        end
    elseif iscell(COLORMAP)
        % Cell array of RGB triplets → convert to matrix
        colorStops = cell2mat(COLORMAP(:));
    elseif isnumeric(COLORMAP)
        % Numeric array directly given
        colorStops = COLORMAP;
    else
        error('COLORMAP must be a name, a cell array of RGBs, or an n×3 numeric matrix.');
    end

    % Safety check
    if size(colorStops, 2) ~= 3
        error('Each color stop must be a 1×3 RGB vector.');
    end
    if size(colorStops, 1) < 2
        error('At least two colors are required for a gradient.');
    end

    % Interpolate along color stops
    colorPositions = linspace(0, 1, size(colorStops, 1));
    interpPositions = linspace(0, 1, LOOPLENGTH);
    gradient = zeros(LOOPLENGTH, 3);
    for c = 1:3
        gradient(:, c) = interp1(colorPositions, colorStops(:, c), interpPositions, 'linear');
    end

    % Select color for this INDEX
    INDEX = max(1, min(LOOPLENGTH, INDEX)); % clamp
    LINECOLOR = gradient(INDEX, :);

%%% Marker selection
    % Default markers
    mark = strings(2,1);
    mark(1) = "*"; mark(2) = "<"; mark(3) = "square";
    mark(4) = ">"; mark(5) = "x"; mark(6) = "diamond";
    mark(7) = "^"; mark(8) = "+"; mark(9) = "pentagram";
    mark(10) = "v"; mark(11) = "o"; mark(12) = "hexagram"; 
    mark(13) = "_"; mark(14) = "|"; mark(15) = "x"; mark(16) = ".";
    
    % Set marker type
    ind = mod(INDEX,length(mark));
    if ind == 0 % if point is a multiple of 16
        MARKER = '.'; MARKSIZE = 20; % otherwise, uses '.', size 20
    else
        MARKER = mark(ind); MARKSIZE = 8; % use unique marker and markersize 8
    end
end
