% All-in-one figure saving
function figsave(FIGURE,FILENAME,FILETYPE)
%%% This function saves figures (unifying the "saveas" and "print" 
% functions of MATLAB).

    FILENAME = string(FILENAME);

    if ischar(FILETYPE)
        FILETYPE = string(FILETYPE); % convert '' to ""
    end

    if isstring(FILETYPE)
        if contains(FILETYPE,"svg","IgnoreCase",true) % svg
            FILETYPE = 1;
        else
            if contains(FILETYPE,"eps","IgnoreCase",true) % eps
                FILETYPE = 2;
            else
                if contains(FILETYPE,"pdf","IgnoreCase",true) % pdf
                    FILETYPE = 3;
                else % png or other
                    FILETYPE = 0;
                end
            end
        end
    end

    if FILETYPE == 0 % png, 600x600
        print('-r600',FIGURE,FILENAME,'-dpng')
    else
        if FILETYPE == 1 % svg
            saveas(FIGURE,string(FILENAME+".svg"),'svg');
        else
            if FILETYPE == 2 % eps
                saveas(FIGURE,string(FILENAME+".eps"),'eps');
            else % PDF
                saveas(FIGURE,string(FILENAME+".pdf"),'pdf');
            end
        end
    end

    return
end