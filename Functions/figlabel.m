% Automatic labeling of figure subpanels
function figlabel(LABEL)
%%% This function creates a label at the top-left corner for a figure subpanel.

    if ~ischar(LABEL)
        LABEL = string(LABEL);
        LABEL = char(LABEL);
    end
    
    annotation('textbox',[0 0.8 0.2 0.2],'String',LABEL,...
        'FitBoxToText','on','LineStyle','none','FontName','Arial','FontSize',25);
    return
end
