% Load and trim data one channel at a time
function [DELAYPOS,DATA,IMAGESIZE] = ...
    loaddata(WORKINGDIRECTORY,SUBFOLDERS,II,CURRENTFILE,CURRFILETYPE,...
    TRIMLASTT,TRIMFIRSTT,XSTEPS,YSTEPS,GUESSTLIST)
%%% This function loads data from .txt, .raw, or .tif files and
% automatically trims them as specified/necessary.

    if CURRFILETYPE == 1 % load data from .txt
        IMAGESIZE = [1 1]; % only solution data
        DATA = importdata(CURRENTFILE); % load data
        DELAYPOS = DATA(:,1); % save delay position as x
        if DELAYPOS(1) == DELAYPOS(end) && TRIMLASTT == 0 % if x isn't single-valued
            TRIMLASTT = 1; % trimming is needed
        end
        DELAYPOS = DELAYPOS(TRIMFIRSTT+1:end-TRIMLASTT);
    else % load Tlist and image
        % Look for Tlist
        if isfile(GUESSTLIST) % check for per-file Tlist first
            DELAYPOS = importdata(GUESSTLIST);
            trueTlist = 1;
        else
            Tlistname = fullfile(WORKINGDIRECTORY,SUBFOLDERS{II},'Tlist.txt');
            if isfile(Tlistname)
                DELAYPOS = importdata(Tlistname); % load Tlist
                trueTlist = 1;
            else % guess Tlist - assumes 0.05 mm spacing
                warning('Tlist not found - generating Tlist with 0.05 mm spacing')
                DELAYPOS = linspace(176,180,81); DELAYPOS(end+1) = 176;
                trueTlist = 0;
            end
        end
        if CURRFILETYPE < 4 % (already ~=1) -> .tif
            DATA = double(tiffreadVolume(CURRENTFILE));
        else % .raw
            imgid = fopen(CURRENTFILE);
            imgvector = fread(imgid,'real*8'); % load raw file as a vector
            fclose('all');
            matrixarea = length(imgvector)/length(DELAYPOS);
            if matrixarea == XSTEPS*YSTEPS
                % data match expected size - do nothing
            else
                if mod(sqrt(matrixarea),1) == 0 % if data are a square
                    XSTEPS = sqrt(matrixarea); YSTEPS = sqrt(matrixarea);
                else
                    if CURRFILETYPE == 3 || CURRFILETYPE == 5
                        if trueTlist == 1
                            warning('Mismatch between dimensions of .raw file and filename.')
                        else
                            warning('Mismatch in dimensions of file - check Tlist specification and file dimensions.')
                        end
                    end
                end
            end
            tempdata = reshape(imgvector,[YSTEPS,XSTEPS,length(DELAYPOS)]);
            DATA = zeros(XSTEPS,YSTEPS,length(DELAYPOS));
            for i=1:length(DELAYPOS)
                DATA(:,:,i) = tempdata(:,:,i).';
            end
        end
        if DELAYPOS(1) == DELAYPOS(end) && TRIMLASTT == 0 % if x isn't single-valued
            TRIMLASTT = 1; % trimming is needed
        end
        DELAYPOS = DELAYPOS(TRIMFIRSTT+1:end-TRIMLASTT);
        DATA = DATA(:,:,TRIMFIRSTT+1:end-TRIMLASTT);
        IMAGESIZE = size(DATA);
        if length(DELAYPOS) ~= IMAGESIZE(3)
            warning('Tlist mismatch - trimming to correct')
            if length(DELAYPOS) > length(DATA)
                DELAYPOS = DELAYPOS(1:length(DATA));
            else
                DATA = DATA(:,:,1:length(DELAYPOS));
            end             
        end
    end
    return
end