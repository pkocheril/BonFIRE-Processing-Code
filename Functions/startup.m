% MATLAB preferences

% Get version number
versionnumber = version;
year = split(versionnumber,'.');
year = sscanf(string(year(1)),'%f');

if year >= 25 % if MATLAB R2025a or later
    set(groot, "defaultFigurePosition", [680 458 560 420]) % figure dimensions
    % set(groot, "defaultFigureWindowStyle", "normal") % separate windows for figures
    %%% can't be closed with Cmd+W
end

% Make figure defaults same as R2024b
% set(groot, "defaultFigurePosition", [680 458 560 420]) % figure dimensions
% set(groot, "defaultFigureWindowStyle", "normal") % separate windows for figures

%%% Past attempts
% set(groot, "defaultFigurePosition", [680 458 800 700]) % figure dimensions
% % Default figure manual position: drag to top left, set right edge to
% % edge of Layout icon and bottom edge to top of Command Window
% Other found preferences
% set(groot,'defaultfigureposition',[400 250 900 750])
% set(groot, 'defaultAxesFontSize', 15);
% set(groot, 'defaultFigureUnits','normalized')
% set(groot, 'defaultFigurePosition',[0 0 0.4 0.4])

clear; close all; clc;
