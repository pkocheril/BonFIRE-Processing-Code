% Load spontaneous Raman data
function [X,Y,CS] = sploaddata(CURRENTFILE)
%%% This function loads spontaneous Raman data from a given input file.

    % If file doesn't exist or no file specified, use file browser GUI
    if isempty(CURRENTFILE) || ~isfile(CURRENTFILE)
        [filename,folder] = uigetfile('*.*');
        CURRENTFILE = fullfile(folder,filename);
    end

    % Read .txt file as cell
    celldata = readcell(CURRENTFILE);

    if contains(char(celldata(1,1)),'#') % file has header
        data = cell2mat(celldata(38:end,:));
        % Read header and assign values to structure
        CS.info.labels = string(celldata(1:37,1));
        CS.info.values = string(celldata(1:37,2));
    else % no header
        data = cell2mat(celldata);
    end

    CS.data.x = data(:,1); CS.data.y = data(:,2);
    X = data(:,1); Y = data(:,2);
    
    return
end
