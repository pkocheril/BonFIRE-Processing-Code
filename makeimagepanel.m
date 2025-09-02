% Making images
function makeimagepanel(array,XSTEPS,YSTEPS,cblabel)
%%% This function plots an image panel with my default preferences.

    image(array,'CDataMapping','Scaled'); 
    daspect([1 1 1]); xticks([]); yticks([]); cb = colorbar;
    cb.Label.Rotation=270; cb.Label.VerticalAlignment = "bottom"; cb.FontSize = 12;
    xlim([1 YSTEPS]); ylim([1 XSTEPS-2]); % crop out bottom two rows
    ax = gca; ax.FontSize = 12;
    if cblabel == 1
        cb.Label.String='BonFIRE (AU)';
        ax.CLim = [0 Inf];
    end
    if cblabel == 2
        cb.Label.String='τ_{1} (ps)';
        ax.CLim = [0 1.42];
    end
    if cblabel == 3
        cb.Label.String='τ_{2} (ps)';
    end
    if cblabel == 4
        cb.Label.String='A_{1}/A_{2}';
    end
    return
end