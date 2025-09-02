% Automatic color mapping and marker selection
function [LINECOLOR,MARKER,MARKSIZE] = colormarkset(INDEX,LOOPLENGTH,COLORMAP)
%%% This function automatically chooses colors and marker sets for plotting
% data in a loop based on sets of pre-defined colors and markers.
    
    % Default colors
    colorset(1,:) = [0 0.35 0.5]; % turquoise
    colorset(2,:) = [0.8 0.2 0.2]; % red
    colorset(3,:) = [0.1 0.7 0.5]; % seafoam green
    colorset(4,:) = [0.8 0.1 0.6]; % pink
    colorset(5,:) = [1 0.5 0.2]; % orange
    colorset(6,:) = [0.2 0.8 0.8]; % light blue
    colorset(7,:) = [0.8 0.8 0.2]; % yellow
    colorset(8,:) = [0.15 0.6 1]; % moderate blue
    colorset(9,:) = [0.2 0.8 0.2]; % green
    colorset(10,:) = [0.5 0.5 1]; % lavender
    colorset(11,:) = [1 0.35 0.35]; % salmon
    colorset(12,:) = [0.8 0.2 0.8]; % magenta
    colorset(13,:) = [0.5 0.8 0.5]; % pale green
    colorset(14,:) = [0.1 0.4 0.8]; % blue
    colorset(15,:) = [0.8 0.6 0.1]; % yellow-orange
    colorset(16,:) = [0.5 0 0.5]; % violet
    colorset(17,:) = [0.6 0.8 0.1]; % yellow-green
    colorset(18,:) = [0.6 0.1 0.8]; % purple
    colorset(19,:) = [0 0 1]; % dark blue
    colorset(20,:) = [0.2 0.2 0.8]; % indigo
    colorset(21,:) = [0.3 0 0.3]; % dark purple
    colorset(22,:) = [0.5 0.2 0.5]; % misc
    colorset(23,:) = [0 0.5 0.5]; % misc
    colorset(24,:) = [0.2 0.6 0.2]; % misc
    colorset(25,:) = [0.6 0.2 0.6]; % misc
    colorset(26,:) = [0.6 0.6 0.2]; % misc
    colorset(27,:) = [0.2 0.8 0.8]; % misc
    colorset(28,:) = [0.8 0 0.8]; % misc
    colorset(29,:) = [0 160 255]./255; % ATTO740
    colorset(30,:) = [0 80 255]./255; % MARS2233
    colorset(31,:) = [80 0 255]./255; % MARS2238
    colorset(32,:) = [160 0 255]./255; % O-NHS
    colorset(33,:) = [0 128 0]./255; % ATTO740
    colorset(34,:) = [0 255 0]./255; % MARS2233
    colorset(35,:) = [0 255 160]./255; % MARS2238
    colorset(36,:) = [0 128 128]./255; % O-NHS
    colorset(37,:) = [255 80 0]./255; % ATTO740
    colorset(38,:) = [255 160 0]./255; % MARS2233
    colorset(39,:) = [255 228 0]./255; % MARS2238
    colorset(40,:) = [160 255 0]./255; % O-NHS
    colorset(41,:) = [128 0 80]./255; % ATTO740
    colorset(42,:) = [255 0 160]./255; % MARS2233
    colorset(43,:) = [255 0 80]./255; % MARS2238
    colorset(44,:) = [255 128 128]./255; % O-NHS
    colorset(45,:) = [0.15 0.6 1]; % ATTO725 nitrile
    colorset(46,:) = [0.2 0.8 0.8]; % Rh800 nitrile
    colorset(47,:) = [0.8 0.2 0.2]; % 13C15N Rh800 nitrile
    
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
    % if INDEX <= length(mark) % if less than 15
    %     MARKER = mark(INDEX); MARKSIZE = 8; % use unique marker and markersize 8
    % else
    %     MARKER = '.'; MARKSIZE = 20; % otherwise, uses '.', size 20
    % end

    % Set color scheme
    if (isempty(COLORMAP) || strcmp(COLORMAP,'default')) && LOOPLENGTH <= height(colorset) % default colorset
        LINECOLOR = colorset(INDEX,:); % if total set can be described with specified number of colors
    else
        if (isempty(COLORMAP) || strcmp(COLORMAP,'default')) % not enough colors in default list
            COLORMAP = 'rainbow'; % default to rainbow
        end
        if strcmp(COLORMAP,'rainbow') % red-green-blue rainbow (optimized)
            redamt = exp((INDEX-LOOPLENGTH)*2/LOOPLENGTH);
            greenamt = ((-4/LOOPLENGTH^2)*(INDEX-0.5*LOOPLENGTH)^2+1)^2;
            blueamt = exp(-2*INDEX/LOOPLENGTH);
            LINECOLOR = [redamt greenamt blueamt];
        else % not rainbow
            startcolor = [1 0 0]; peakcolor = [0 1 0]; endcolor = [0 0 1]; % if unrecognized colormap, default to rainbow (unoptimized)
            if strcmp(COLORMAP,'fire') % orange-red-purple gradient
                startcolor = [255 134 29]./255; peakcolor = [200 60 60]./255; endcolor = [100 0 100]./255;
            end
            if strcmp(COLORMAP,'fireinv') % inverted orange-red-purple gradient
                startcolor = [100 0 100]./255; peakcolor = [200 60 60]./255; endcolor = [255 134 29]./255;
            end
            if strcmp(COLORMAP,'rwb1') || strcmp(COLORMAP,'unionjack')  % red-white-blue gradient 1
                startcolor = [231 48 41]./255; peakcolor = [235 242 250]./255; endcolor = [48 124 188]./255;
            end
            if strcmp(COLORMAP,'rwb2') % red-white-blue gradient 2
                startcolor = [103 0 30]./255; peakcolor = [240 240 240]./255; endcolor = [6 48 97]./255;
            end
            if strcmp(COLORMAP,'rpb') % red-purple-blue gradient
                startcolor = [200 0 30]./255; peakcolor = [100 0 170]./255; endcolor = [0 60 200]./255;
            end
            if strcmp(COLORMAP,'rwb1inv')  % inverted red-white-blue gradient 1
                startcolor = [48 124 188]./255; peakcolor = [235 242 250]./255; endcolor = [231 48 41]./255;
            end
            if strcmp(COLORMAP,'parula') % blue-green-yellow gradient
                startcolor = [68 1 84]./255; peakcolor = [38 138 141]./255; endcolor = [253 231 36]./255;
            end
            if strcmp(COLORMAP,'parulainv') % inverted blue-green-yellow gradient
                startcolor = [253 231 36]./255; peakcolor = [38 138 141]./255; endcolor = [68 1 84]./255;
            end
            if strcmp(COLORMAP,'parula2') % alternate yellow-green-blue gradient
                startcolor = [240 249 185]./255; peakcolor = [65 182 196]./255; endcolor = [6 29 88]./255;
            end
            if strcmp(COLORMAP,'DFG') % purple-white-green gradient (DFG contour)
                startcolor = [0.5 0 0.5]; peakcolor = [1 1 1]; endcolor = [0.2 0.8 0.2];
            end
            if strcmp(COLORMAP,'idler') % turquoise-white-red gradient (idler contour)
                startcolor = [0 0.35 0.5]; peakcolor = [1 1 1]; endcolor = [0.8 0.2 0.2];
            end
            if strcmp(COLORMAP,'probe') % blue-white-orange gradient (probe contour)
                startcolor = [0 0 1]; peakcolor = [1 1 1]; endcolor = [1 0.5 0];
            end
            if mod(LOOPLENGTH,2) == 1 % maplength is odd
                segmentlength = (LOOPLENGTH+1)./2;
                map1 = [linspace(startcolor(1),peakcolor(1),segmentlength) linspace(peakcolor(1),endcolor(1),segmentlength)]; map1(segmentlength) = [];
                map2 = [linspace(startcolor(2),peakcolor(2),segmentlength) linspace(peakcolor(2),endcolor(2),segmentlength)]; map2(segmentlength) = [];
                map3 = [linspace(startcolor(3),peakcolor(3),segmentlength) linspace(peakcolor(3),endcolor(3),segmentlength)]; map3(segmentlength) = [];
            else % maplength is even
                segmentlength = (LOOPLENGTH./2)+1;
                map1 = [linspace(startcolor(1),peakcolor(1),segmentlength) linspace(peakcolor(1),endcolor(1),segmentlength)]; map1(segmentlength) = [];
                map2 = [linspace(startcolor(2),peakcolor(2),segmentlength) linspace(peakcolor(2),endcolor(2),segmentlength)]; map2(segmentlength) = [];
                map3 = [linspace(startcolor(3),peakcolor(3),segmentlength) linspace(peakcolor(3),endcolor(3),segmentlength)]; map3(segmentlength) = [];
            end
            map = [map1.' map2.' map3.'];
            LINECOLOR = map(INDEX,:);
        end
    end
    return
end