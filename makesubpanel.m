% Plotting time-domain spectra
function makesubpanel(T,SIGNAL,SIGSDS,TBASE,SIGBASE,BASECURVE,TFIT,FITCURVE,...
    RESID,LABEL,LEGEND,CHANNEL)
%%% This function makes subplots for time-domain spectra given raw data and
% fits.

    if isempty(SIGSDS)
        SIGSDS = zeros(length(SIGNAL),1);
    end
    % Make colors for figure
    plotcolor(1,:) = [0.5 0.5 0.5]; % Data, gray
    plotcolor(2,:) = [0.9 0.1 0.1]; % CH1 baseline data
    plotcolor(3,:) = [0.7 0.4 0.1]; % CH1 baseline fit
    plotcolor(4,:) = [0.2 0.8 0.8]; % CH2 baseline data
    plotcolor(5,:) = [0.2 0.8 0.5]; % CH2 baseline fit
    plotcolor(6,:) = [0.2 0.2 0.8]; % Lifetime fit
    plotcolor(7,:) = [0 1 1];
    plotcolor(8,:) = [174 110 180]./255; % Residuals
    % New Excel colors
    plotcolor(9,:) = [22 96 130]./255;
    plotcolor(10,:) = [233 113 49]./255;
    plotcolor(11,:) = [26 107 37]./255;
    plotcolor(12,:) = [14 158 213]./255;
    plotcolor(13,:) = [160 44 147]./255;
    plotcolor(14,:) = [77 167 46]./255;
    % New Parula colors
    plotcolor(15,:) = [68 1 84]./255;
    plotcolor(16,:) = [71 46 123]./255;
    plotcolor(17,:) = [54 96 141]./255;
    plotcolor(18,:) = [38 138 141]./255;
    plotcolor(19,:) = [53 179 124]./255;
    plotcolor(20,:) = [143 212 72]./255;
    plotcolor(21,:) = [253 231 36]./255;
    col1 = plotcolor(1,:);
    if CHANNEL == 2 % set colors by channel
        % col1 = plotcolor(15,:); % data
        col2 = plotcolor(20,:); % baseline data
        col3 = plotcolor(19,:); % baseline fit
        col4 = plotcolor(18,:); % fit
        col5 = plotcolor(15,:); % residuals
    else
        if CHANNEL == 1
            col2 = plotcolor(10,:); % baseline data
            col3 = plotcolor(9,:); % baseline fit
            col4 = plotcolor(11,:); % fit
            col5 = plotcolor(13,:); % residuals
        else
            if CHANNEL == 3
                col2 = plotcolor(4,:);
                col3 = plotcolor(5,:);
                col4 = plotcolor(6,:);
                col5 = plotcolor(8,:);
            else
                col2 = plotcolor(2,:);
                col3 = plotcolor(3,:);
                col4 = plotcolor(13,:);
                col5 = plotcolor(14,:);
            end
        end
    end
    hold on; box on; xlim([min(T) max(T)]); xlabel('Time delay (ps)'); ylabel(LABEL);
    % Make legend
    plot([min(T)-10 min(T)-10],[min(T)-10 min(T)-10],'o','Color',col1,'LineWidth',2);
    plot([min(T)-10 min(T)-10],[min(T)-10 min(T)-10],'o','Color',col2,'LineWidth',2); 
    plot([min(T)-10 min(T)-10],[min(T)-10 min(T)-10],'-','Color',col3,'LineWidth',2);
    % Plot residuals
    if sum(SIGNAL)+length(SIGNAL) ~= sum(RESID)+length(RESID) % plot fit if truly fitted (resid = signal if no fit)
        yyaxis right; % plot residuals on secondary axis
        plot(T,RESID,'Color',col5,'LineWidth',1.3); LEGEND(4) = {'Fit'};  LEGEND(5) = {'Residuals'};
        ylabel('Residuals (AU)'); yyaxis left; 
        plot([min(T)-10 min(T)-10],[min(T)-10 min(T)-10],'-','Color',col4,'LineWidth',2);
    end
    legend(LEGEND); lg = legend; lg.EdgeColor = [1 1 1]; lg.AutoUpdate = 'off';
    % Plot data
    patch([T; flip(T)],[SIGNAL-SIGSDS; flip(SIGNAL+SIGSDS)],'k','FaceColor',col1,'FaceAlpha',0.5,'EdgeColor','none')
    if ~isempty(BASECURVE)
        plot(T,BASECURVE,'-','Color',col3,'LineWidth',3);
    end
    if sum(SIGNAL)+length(SIGNAL) ~= sum(RESID)+length(RESID) % plot fit if truly fitted (resid = signal if no fit)
        plot(TFIT,FITCURVE,'-','Color',col4,'LineWidth',2.5);
    end
    plot(T,SIGNAL,'o','Color',col1,'LineWidth',2); 
    if ~isempty(SIGBASE)
        plot(TBASE,SIGBASE,'o','Color',col2,'LineWidth',2);
    end
    ax = gca; ax.FontSize = 12; ax.LineWidth = 2; hold off;
    return
end