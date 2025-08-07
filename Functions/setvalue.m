% Set a definable value per-folder or globally
function current = setvalue(input,FOLDERS,II)
%%% This function is used to prepare the pre-set value for a given parameter
% (eg, pulse width) either on a per-folder basis or globally

    if ~isempty(input) % if specified
        if length(input) == length(FOLDERS) % if specified per-folder
            currentfolder = (FOLDERS == II); % get current folder position
            current = input(currentfolder);
        else
            current = input(1); % if specified as single value
        end
    else
        current = []; % unspecified
    end
    return
end
