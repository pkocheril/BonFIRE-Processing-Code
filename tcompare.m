% Time delay sweep comparison
function tcompare(TIME,SIG,NAMES,COLORMAP,OPTION)
%%% This function overlays a series of time sweeps in a new figure.
    
    % Make sure TIME is a row vector
    if length(TIME) == height(TIME)
        TIME = TIME.';
    end
    % Check if TIME and SIG match dimensions
    if width(TIME) ~= width(SIG) % if mismatched
        if width(TIME) == height(SIG) % SIG is flipped
            SIG = SIG.';
        else % trim SIG to match
            SIG = SIG(:,1:length(TIME));
            warning('Dimensions mismatched - signal matrix trimmed to match.');
        end
    end

    if isempty(COLORMAP)
        COLORMAP = 'fire';
    end
    figure; hold on; box on;
    for i=1:height(SIG)
        [LINECOLOR,MARKER,MARKSIZE] = colormarkset(i,height(SIG),COLORMAP);
        plot([-1 -1],[-1 -1],'Marker',MARKER,'MarkerSize',MARKSIZE,'Color',LINECOLOR,'LineWidth',2);
    end
    if ~isempty(NAMES)
        legend(NAMES); lg = legend; lg.Box = 'off'; lg.AutoUpdate = 'off';
    end
    for i=1:height(SIG)
        [LINECOLOR,MARKER,MARKSIZE] = colormarkset(i,height(SIG),COLORMAP);
        if height(TIME) == height(SIG)
            t = TIME(i,:);
        else
            t = TIME(1,:);
        end
        if strcmp(OPTION,'normalize')
            signal = normsweep(SIG(i,:));
            ymax = 1.05; ystring = 'Normalized BonFIRE';
        else
            signal = SIG(i,:);
            ymax = Inf; ystring = 'BonFIRE (AU)';
        end
        [TFIT,FITCURVE] = basicltfit(t,signal,2);
        [~,MAXIND] = max(FITCURVE);
        offset = TFIT(MAXIND);
        if strcmp(OPTION,'normalize')
            signal = (signal-min(FITCURVE))./(max(FITCURVE)-min(FITCURVE));
            FITCURVE = (FITCURVE-min(FITCURVE))./(max(FITCURVE)-min(FITCURVE));
        end
        plot(t-offset,signal,'o','Marker',MARKER,'MarkerSize',MARKSIZE,'Color',LINECOLOR,'LineWidth',2);
        plot(TFIT-offset,FITCURVE,'-','Color',LINECOLOR,'LineWidth',3);
    end
    xlabel('Time delay (ps)'); ylabel(ystring);
    xlim([min(TIME,[],'all') max(TIME,[],'all')]); ylim([-0.05 ymax]);
    ax = gca; ax.FontSize = 15; ax.LineWidth = 2;
    return
end