%% Batch processing of 2DVF spectra and BonFIRE images 
%%% v2 - improved auto-lifetime fit selection
%%% v3 - bug fixes

% Initialize
cd '/Users/pkocheril/Documents/Caltech/Wei Lab/Data/2024_07_02_PK/'
clear; clc; close all;

% Main configuration options
loadprevious = 0; % [] = auto-detect, 0 = new analysis, 
% 1 = load individual .dats, 2 = read structure.xml
runtype = 0; % 0 = process all, 1 = test run with a few files,
% 2 = examine a single file, 3 = examine a single image
targetfolders = []; % indices of folders to process, [] = dialog
t0pos = 290.4; % specify t0 position (mm), [] = autofind

% Additional configuration options
ltfittype = []; % [] = auto-choose, 0 = no fitting, 1 = Gaussian*monoexp,
% 2 = Gaussian*biexp, 3 = Gauss*stretchexp, 4 = Gauss*strbiexp
basefittype = []; % [] = auto-choose, 0 = no baseline fit, 1 = linear, 
% 2 = exponential, 3 = exponential+linear
fitchannels = 2; % specify data channels to fit, [] = dialog
writeprocyn = 1; % 1 = write batch processed files, 0 = not
writefigsyn = 1; % 1 = write figure files, 0 = not
powernormtype = 1; % 0 = no normalization, 1 = normalize by IR power,
% 2 = normalize by probe and IR powers, 3 = 2 + PMT gain correction
setpulsewidth = []; % define pulse width (ps) in fit, [] = float
floatbase = 0.1; % fraction to float baseline coeffs (0.1 -> +/- 10%)
cutlow = 0.05; % lower baseline fit cutoff, (default 5th %ile)
cuthigh = 0.8; % upper baseline fit cutoff, (default 80th %ile)
verbose = 2; % 2 = all info on figure, 1 = regular figure annotation, 
% 0 = no figure annotation
showprogress = 1; % show progress as fitting proceeds over pixels in image
snrcutoff = -Inf; % SNR minimum for lifetime image fitting (-Inf = fit all)
troubleshoot = 0; % 0 = default, 1 = show everything along the way

% Even more configuration options (largely obsolete)
tempmod = 0; % 1 = temporal modulation in place, 0 = no beamsplitters
ND = 0; % ND filter strength in probe path (auto-update if possible)
pmtgain = 1; % PMT gain (norm to 1, auto-update if possible)
pmtBW = 0; % PMT bandwidth in Hz (auto-update if possible)
pinholeyn = 0; % 0 = bypassed, 1 = in place
trimlastT = 0; % how many points to remove from end of Tlist (default 0)
trimfirstT = 0; % points to remove from start of Tlist (default 0)
ltmin = [];  % define minimum lifetime (ps) in fit, [] = default (0.1)
ltmax = []; % define maximum lifetime (ps) in fit, [] = default (100)
tuningrate = []; % PicoEmerald signal t0 tuning rate (ps/nm), [] = auto
minprobe = []; maxprobe = []; % probe range for aligning tD, [] = auto-find
% pairedDCAC = []; % [] = auto-assign, 2 = SPCM + PMT CH2, 
% % 1 = PMT CH1 + CH2, 0 = only PMT CH2
setfiletype = []; % [] = auto-assign, 1 = .txt, 2 = .tif (solution), 
% 3 = .tif (image), 4 = .raw (solution), 5 = .raw (image)
setdatatype = []; % [] = auto-assign, 0 = solution, 1 = image
setfilenameconv = []; % [] = auto-assign, 0 = other,
% 1 = [pumpWL]-[pumppower]-[signalWL]-[IRpower]-[ND]-[PMTgain]-[modfreq]-[PMTBW]-[etc],
% 2 = [etc]_[size]_[idler]_[DFG]_[power]_[channel],
% 3 = [IRWN]_[probeWL]_[etc]_[size]_[idler]_[DFG]_[power]_[channel],
% 4 = [probeWL]_[IRWN]_[etc]_[size]_[idler]_[DFG]_[power]_[channel]
prbdichroic = 0; pmtcmosmirror = 0; pmtfilter = 0; spcmfilter = 0; % experimental parameters (auto-update if possible)
cmosfilter = 0; dutycycle = 0; probeFWHM = 0; levanteFWHM = 0; % FWHMs in nm
testfilesperfolder = 1; % number of files to test per folder (default 1)

workingdirectory = pwd; % get working directory path
directorycontents = dir(fullfile(workingdirectory,'*')); % search working directory
subfolders = setdiff({directorycontents([directorycontents.isdir]).name},{'.','..'}); % subfolders of working directory

if isempty(loadprevious) % auto-detect previous analysis
    guessxml = fullfile(workingdirectory,'summarystructure.xml');
    if isfile(guessxml)
        loadprevious = 2;
    else
        guessyml = fullfile(workingdirectory,'add_info.yml');
        if isfile(guessyml)
            loadprevious = 1;
        else
            loadprevious = 0;
        end
    end
end

if troubleshoot == 1
    runtype = 2;
end

if runtype > 1 % to examine a single file - choose one folder
    [indx,~] = listdlg('PromptString',{'Select one folder.'},...
    'SelectionMode','single','ListString',subfolders);
    targetfolders = indx;
end

if isempty(fitchannels) % dialog to select fit channels
    [indx,~] = listdlg('PromptString',...
        {'Select data channels for lifetime fitting.'},'SelectionMode',...
        'multiple','ListString',{'CH1','CH2','CH3','CH4','SPCM'},'InitialValue',2);
    fitchannels = indx;
end

if isempty(targetfolders) % dialog to select folders
    [indx,~] = listdlg('PromptString',{'Select folders to process.'},...
    'SelectionMode','multiple','ListString',subfolders);
    targetfolders = indx;
end

if runtype > 0  % setup test run conditions
    visibility = 'on'; % show figures
    writeprocyn = 0; % don't write batch processed files
    writefigsyn = 0; % don't write figure files
    if length(targetfolders) >= 3 % folder logic for test run
        folders = targetfolders(1:3); % limit to three folders
    else
        folders = targetfolders;
    end
    if powernormtype > 2
        powernormtype = 2;
    end
else % not a test run
    visibility = 'off'; % hide figures
    if ~isempty(targetfolders) % if target folders are specified
        folders = targetfolders;
    else
        folders = 1:numel(subfolders); % all folders
    end
end

if loadprevious == 0 % run new analysis
    for ii = folders % loop through folders
        % Set up sub-structure for current folder
        summary.('folder'+string(ii)).foldername = subfolders{ii};

        % Set current values per-folder (when specified)
        currt0 = setvalue(t0pos,folders,ii); currltfittype = setvalue(ltfittype, folders,ii);
        currbasefittype = setvalue(basefittype,folders,ii); currdatatype = setvalue(setdatatype,folders,ii);
        currpulsewidth = setvalue(setpulsewidth,folders,ii); currfiletype = setvalue(setfiletype,folders,ii);
        currfilenameconv = setvalue(setfilenameconv,folders,ii);

        % Look for files
        rawlist = dir(fullfile(workingdirectory,subfolders{ii},'*.raw')); rawlist = {rawlist(~[rawlist.isdir]).name};
        CH2raw = dir(fullfile(workingdirectory,subfolders{ii},'*CH2.raw')); CH2raw = {CH2raw(~[CH2raw.isdir]).name};
        txtlist = dir(fullfile(workingdirectory,subfolders{ii},'*.txt')); txtlist = {txtlist(~[txtlist.isdir]).name};
        tiflist = dir(fullfile(workingdirectory,subfolders{ii},'*.tif')); tiflist = {tiflist(~[tiflist.isdir]).name};
        if isempty(CH2raw) % CH2 raw files not found
            if isempty(rawlist) % no .raw files
                if isempty(tiflist) % no .tif files
                    datafiles = txtlist; % data are .txts
                else % found .tifs
                    datafiles = tiflist;
                end
            else
                datafiles = rawlist; % all .raw files
                warning('No CH2.raw files found - processing all .raw files')
            end
        else % CH2.raw files found
            datafiles = CH2raw;
        end

        infofilesraw = dir(fullfile(workingdirectory,subfolders{ii},'*IRsweep.txt'));
        if isempty(infofilesraw) % if no IRsweeps
            infofilesraw = dir(fullfile(workingdirectory,subfolders{ii},'*XYZT.txt')); % look for XYZT info
            if isempty(infofilesraw)
                infofilesraw = dir(fullfile(workingdirectory,subfolders{ii},'*Z.txt')); % look for Z stack info
                if isempty(infofilesraw)
                    infofilesraw = dir(fullfile(workingdirectory,subfolders{ii},'*T.txt')); % look for T stack info
                end
            end
        end
        infofiles = {infofilesraw(~[infofilesraw.isdir]).name}; % info files
        if numel(infofiles) >= numel(datafiles) % can have extra infofiles during scanning
            infofound = 1;
        else
            infofound = 0;
        end
        if runtype > 1 % to examine a single file - choose file
            [indx2,~] = listdlg('PromptString',{'Select a file to process.'},...
            'SelectionMode','single','ListString',datafiles);
            targetfile = indx2;
            totalfiles = targetfile; % load only chosen file
            startfile = targetfile;
        else
            startfile = 1;
            if runtype == 1 % limit test run files per folder
                totalfiles = testfilesperfolder;
            else
                totalfiles = numel(datafiles); % all files
            end
        end
        summary.('folder'+string(ii)).nfiles = totalfiles;

        if ~isempty(datafiles) % skip subfolders without data files
            for jj = startfile:totalfiles % jj = file number in subfolder ii
                currentfile = fullfile(workingdirectory,subfolders{ii},datafiles{jj}); % full file path
                
                % Filename parsing
                macseparator = count(currentfile,'/');
                winseparator = count(currentfile,'\');
                if macseparator > winseparator % MacOS/UNIX
                    delimiter = "/";
                else % Windows
                    delimiter = "\";
                end
                folderinfo = split(currentfile,delimiter); % split path into folders
                filename = char(folderinfo(end));

                % Figure out solution vs image data by filename
                if isempty(currdatatype)
                    if ~contains(currentfile,'FOV') % filename doesn't have 'FOV'
                        currdatatype = 0; % solution data
                    else % filename has 'FOV'
                        currdatatype = 1; % image data
                    end
                end

                % Figure out file naming convention
                if isempty(currfilenameconv)
                    if contains(currentfile,'.txt') %count(currentfile,'.txt') > 0
                        lookforhyphens = count(folderinfo(end),'-');
                        if numel(lookforhyphens) > 5
                            currfilenameconv = 1;
                        else
                            currfilenameconv = 0;
                        end
                    else % .tif or .raw
                        if strcmp(filename(5:9),'cm-1_') % filename ####cm-1_...
                            currfilenameconv = 3; % 2DVF probe sweep (before updated VI)
                        else
                            currfilenameconv = 2; % not a 2DVF probe sweep
                        end
                    end
                end

                % Figure out current file type
                if isempty(currfiletype)
                    if contains(currentfile,'.txt')
                        currfiletype = 1;
                    else % raw or tif
                        if contains(currentfile,'.raw')
                            if currdatatype == 0 % solution .raw
                                currfiletype = 4;
                            else
                                currfiletype = 5; % image .raw
                            end
                        else % tif
                            if currdatatype == 0 % solution .tif
                                currfiletype = 2;
                            else
                                currfiletype = 3; % image .tif
                            end
                        end
                    end
                end

                % Empty parameters for new file
                IRpower = []; prbpower = []; xsteps = []; ysteps = []; filebase = '-';
                
                if currfilenameconv == 1 % parse .txts
                    fileinfo = split(folderinfo(end),"-"); % split filename at "-"s
                    subfolderinfo = split(subfolders{ii}," "); % split subfolder by spaces
                    prbWLstrfull = char(fileinfo(1)); prbWLstr = prbWLstrfull(1:end-1); prbWL = sscanf(prbWLstr,'%f'); % extract probe wavelength and convert to double
                    prbpowerstrfull = char(fileinfo(2)); prbpowerstr = prbpowerstrfull(1:end-2); prbpower = sscanf(prbpowerstr,'%f'); % extract probe power
                    signalWLstrfull = char(fileinfo(3)); signalWLstr = signalWLstrfull(1:end-2); signalWL = sscanf(signalWLstr,'%f'); % extract IR signal WL
                    IRpowerstrfull = char(fileinfo(4)); IRpowerstr = IRpowerstrfull(1:end-2); IRpower = sscanf(IRpowerstr,'%f'); % extract IR power
                    DFGWL = 1/((2/signalWL) - (1/1031.2)); DFGWN = 1e7/DFGWL; idlerWN = 1e7*((1/1031.2)-(1/signalWL)); % calculate DFG and idler
                    NDstrfull = char(fileinfo(5)); NDstr = NDstrfull(3:end); ND = sscanf(NDstr,'%f');  % extract ND
                    PMTstrfull = char(fileinfo(6)); PMTstr = PMTstrfull(4:end); pmtgain = sscanf(PMTstr,'%f'); % extract PMT gain
                    modfreqstrfull = char(fileinfo(7)); modfreqstr = modfreqstrfull(1:end-3); % extract modulation freq
                    if modfreqstrfull(end-2) == "k" % kHz
                        modfreq = 1e3*sscanf(modfreqstr,'%f');
                    else % MHz
                        modfreq = 1e6*sscanf(modfreqstr,'%f');
                    end
                    PMTBWstrfull = char(fileinfo(8)); PMTBWstr = PMTBWstrfull(1:end-3); % extract PMT bandwidth
                    if PMTBWstrfull(end-2) == "k" % kHz
                        pmtBW = 1e3*sscanf(PMTBWstr,'%f');
                    else % MHz
                        pmtBW = 1e6*sscanf(PMTBWstr,'%f');
                    end
                else % image file (.raw or .tif)
                    fileinfo = split(folderinfo(end),"_"); % split at "_"s
                    channel = char(fileinfo(end)); % 'CH2'
                    if numel(fileinfo) > 5
                        IRpowermeterstrfull = char(fileinfo(end-1)); IRpowermeterstr = IRpowermeterstrfull(2:end); IRpowermeter = sscanf(IRpowermeterstr,'%f'); % 'P0.4409' --> 0.4409
                        DFGWNstrfull = char(fileinfo(end-2)); DFGWNstr = DFGWNstrfull(4:end); DFGWN = sscanf(DFGWNstr,'%f'); % 'DFG3797.9' --> 3797.9
                        idlerWNstrfull = char(fileinfo(end-3)); idlerWNstr = idlerWNstrfull(6:end); idlerWN = sscanf(idlerWNstr,'%f'); % 'Idler2949.8' --> 2949.8
                        imgdim = char(fileinfo(end-4)); imgdimsplit = split(string(imgdim),"X"); ysteps = imgdimsplit(1); xsteps = imgdimsplit(end); xsteps = sscanf(xsteps,'%f'); ysteps = sscanf(ysteps,'%f'); % '5X5' --> 5 and 5 (Y by X)
                        prbWL = 1; % <-- unknown from filename
                        if IRpowermeter ~= 0
                            IRpower = IRpowermeter*300; % assuming 300 mW scale
                        end
                        filebase = strjoin(fileinfo(1:end-5),'_');
                    else
                        IRpower = 0; DFGWN = 0; idlerWN = 0;
                        if numel(fileinfo) == 4
                            filebase = strjoin(fileinfo(1:end-(numel(fileinfo)-2)),'_');
                        else
                            filebase = strjoin(fileinfo(1:end-(numel(fileinfo)-1)),'_');
                        end
                    end
                end
                if infofound == 1 % parse info from .txt files
                    % Look for matching file first
                    matchIRsweep = folderinfo; matchIRsweep(end) = {strjoin({filebase,'IRsweep.txt'},'_')};
                    matchXYZTsweep = folderinfo; matchXYZTsweep(end) = {strjoin({filebase,'XYZT.txt'},'_')};
                    matchZsweep = folderinfo; matchZsweep(end) = {strjoin({filebase,'Z.txt'},'_')};
                    matchTsweep = folderinfo; matchTsweep(end) = {strjoin({filebase,'T.txt'},'_')};
                    matchIRsweep = strjoin(matchIRsweep,delimiter);
                    matchXYZTsweep = strjoin(matchXYZTsweep,delimiter);
                    matchZsweep = strjoin(matchZsweep,delimiter);
                    matchTsweep = strjoin(matchTsweep,delimiter);
                    if isfile(matchIRsweep)
                        infofilename = matchIRsweep;
                    else
                        if isfile(matchXYZTsweep)
                            infofilename = matchXYZTsweep;
                        else
                            if isfile(matchZsweep)
                                infofilename = matchZsweep;
                            else
                                if isfile(matchTsweep)
                                    infofilename = matchTsweep;
                                else
                                    warning('No matching .txt infofile found - using alphabetization to guess')
                                    infofilename = fullfile(workingdirectory,subfolders{ii},infofiles{jj});
                                end
                            end
                        end
                    end
                    infoarray = importdata(infofilename); infofile = infoarray.data;
                    xinitial = infofile(1); xfinal = infofile(2); xsteps = infofile(3);
                    yinitial = infofile(4); yfinal = infofile(5); ysteps = infofile(6);
                    zinitial = infofile(7); zfinal = infofile(8); zsteps = infofile(9);
                    tinitial = infofile(10); tfinal = infofile(11); tsteps = infofile(12);
                    if round(infofile(13)) ~= 0 % for frequency doubling
                        prbWL = round(infofile(13));
                    else
                        prbWL = 515.6;
                    end
                    if round(infofile(14)) ~= 0
                        prbpower = infofile(14); % set to prbpowerset later if this is zero
                    end
                    ND = infofile(15);
                    idlerWN = round(infofile(16));
                    DFGWN = round(infofile(17));
                    if infofile(18) > 0
                        IRpower = infofile(18);
                    end
                    dwelltime = infofile(19);
                    if length(infofile) > 19 % extended-length files with more info
                        if round(infofile(20)) ~= 0
                            prbWL = infofile(20); % manual probe WL
                        end
                        if round(infofile(21)) ~= 0
                            prbpower = infofile(21); % manual probe power
                        end
                        if length(infofile) > 21 % if using updated code
                            pmtgain = infofile(22); % PMT gain
                            pmtBW = 1e3*infofile(23); % PMT bandwidth (kHz in file)
                            modfreq = 1e3*infofile(24); % IR modulation frequency (kHz in file)
                            pinholeyn = infofile(25); % 0 = bypassed, 1 = in place
                            prbdichroic = infofile(26);
                            pmtcmosmirror = infofile(27);
                            pmtfilter = infofile(28);
                            spcmfilter = infofile(29);
                            cmosfilter = infofile(30);
                        end
                        if length(infofile) > 30
                            dutycycle = infofile(31);
                        end
                        if length(infofile) > 31
                            probeFWHM = infofile(32);
                            levanteFWHM = infofile(33);
                        end
                    else % if not extended info
                        pmtgain = 1; pmtBW = 1; modfreq = 0; pinholeyn = 0;
                        prbdichroic = 0; pmtcmosmirror = 0; pmtfilter = 0;
                        spcmfilter = 0; cmosfilter = 0; dutycycle = 0;
                        probeFWHM = 0; levanteFWHM = 0;
                        if verbose == 2
                            verbose = 1; % can't print out experimental parameters
                        end
                    end
                end
                if isempty(xsteps)
                    xsteps = 1;
                end
                if isempty(ysteps)
                    ysteps = 1;
                end
                % Assign IRWN by idler wavelength
                idlerWL = 1e7/idlerWN;
                signalWL = 1/((1/1031.2)-(1/idlerWL));
                if idlerWL < 2600 % nm
                    IRWN = DFGWN;
                else
                    IRWN = idlerWN;
                end
                if currfilenameconv == 3 % probe sweep
                    prbstrfull = char(fileinfo(2)); % probe WL - "780 (10)"
                    if length(prbstrfull) >= 6 % "780 (10)" --> 790 nm
                        prbsplit = split(prbstrfull," "); % "780" "(10)"
                        if isscalar(prbsplit) % "780finalcheck"
                            prbstr = prbstrfull(1:3); % '780'
                            prbWL = sscanf(prbstr,'%f'); % 780
                        else
                            prbWL1char = char(prbsplit(1)); % '780'
                            prbWL1 = sscanf(prbWL1char,'%f'); % 780
                            prbWL2char = char(prbsplit(2)); % '(10)' or 'finalcheck'
                            if length(prbWL2char) >= 7 % 'finalcheck'
                                prbWL2 = 0;
                            else % '(10)'
                                prbWL2str = prbWL2char(2:end-1); % '10'
                                prbWL2 = sscanf(prbWL2str,'%f'); % 10
                            end
                            if totalfiles <= 130
                                prbWL2 = 2*prbWL2; % step size 2 nm
                            end
                            prbWL = prbWL1 + prbWL2; % 790
                        end
                    else % "780" or "780nm"
                        prbstr = prbstrfull(1:3); % '780'
                        prbWL = sscanf(prbstr,'%f'); % 780
                    end
                end
                if isempty(IRpower) % if no IR power found
                    IRpower = 70; % guess power
                end
                if isempty(prbpower) % if no probe power found
                    prbpower = 300; % guess power
                end

                % Guess Tlist name (per-file)
                guessTlist = folderinfo; guessTlist(end) = {strjoin({filebase,'Tlist.txt'},'_')}; guessTlist = strjoin(guessTlist,delimiter);
                % Load data
                [delaypos,data,imagesize] = ...
                    loaddata(workingdirectory,subfolders,ii,currentfile,...
                    currfiletype,trimlastT,trimfirstT,xsteps,ysteps,guessTlist);
                channels = zeros(5,1); channels(2) = 2; % to keep track of data channels present
                if contains(currentfile,'.txt') % run again for .txt with same file
                    if width(data) > 2
                        ch1data = data(:,2);
                        signal = data(:,3);
                    else
                        signal = data(:,2);
                    end
                    channels(1) = 1;
                end
                % Check for CH1, CH3, CH4, and SPCM files
                ch1file = string(currentfile(1:end-7))+'CH1'+string(currentfile(end-3:end));
                if isfile(ch1file)
                    [~,ch1data,~] = ...
                        loaddata(workingdirectory,subfolders,ii,ch1file,...
                        currfiletype,trimlastT,trimfirstT,xsteps,ysteps,guessTlist);
                    channels(1) = 1;
                end
                ch3file = string(currentfile(1:end-7))+'CH3'+string(currentfile(end-3:end));
                if isfile(ch3file)
                    [~,ch3data,~] = ...
                        loaddata(workingdirectory,subfolders,ii,ch3file,...
                        currfiletype,trimlastT,trimfirstT,xsteps,ysteps,guessTlist);
                    channels(3) = 3;
                end
                ch4file = string(currentfile(1:end-7))+'CH4'+string(currentfile(end-3:end));
                if isfile(ch4file)
                    [~,ch4data,~] = ...
                        loaddata(workingdirectory,subfolders,ii,ch4file,...
                        currfiletype,trimlastT,trimfirstT,xsteps,ysteps,guessTlist);
                    channels(4) = 4;
                end
                spcmfile = string(currentfile(1:end-7))+'SPCM'+string(currentfile(end-3:end));
                if isfile(spcmfile)
                    [~,spcmdata,~] = ...
                        loaddata(workingdirectory,subfolders,ii,spcmfile,...
                        currfiletype,trimlastT,trimfirstT,xsteps,ysteps,guessTlist);
                    channels(5) = 5;
                end
                
                % Prepare for image test run
                if runtype > 0 && currdatatype == 1
                    preview = figure; image(data(:,:,ceil(imagesize(3)/4)),'CDataMapping','Scaled');
                    title('Preview'); colormap('hot'); cb = colorbar; cb.FontSize = 12;
                    cb.Label.String='BonFIRE preview (AU)'; cb.Label.Rotation=270;
                    cb.Label.VerticalAlignment = "bottom"; datacursormode; 
                    xlabel('X'); ylabel('Y'); ax = gca; ax.FontSize = 12;
                    if runtype == 3
                        prompt = {'Specify starting X coordinate.','Specify starting Y coordinate.',...
                            'Specify ending X coordinate.','Specify ending Y coordinate.'};
                        dlgtitle = 'Choose test area'; dims = [1 45; 1 45; 1 45; 1 45]; definput = {'1','1',string(ysteps),string(xsteps)};
                    else
                        prompt = {'Specify X coordinate.','Specify Y coordinate.'};
                        dlgtitle = 'Choose test run pixel'; dims = [1 45; 1 45]; definput = {'0','0'};
                    end
                    coordanswer = inputdlg(prompt,dlgtitle,dims,definput);
                    coordanswer = str2double(coordanswer);
                    close % closes preview image
                    if coordanswer(1) > 0 && coordanswer(1) <= imagesize(2)
                        startcolumn = coordanswer(1);
                    else
                        startcolumn = round(imagesize(2)/2);
                        warning('Specified coordinate out of bounds - using middle of image')
                    end
                    if coordanswer(2) > 0 && coordanswer(2) <= imagesize(1)
                        startrow = coordanswer(2);
                    else
                        startrow = round(imagesize(2)/2);
                        warning('Specified coordinate out of bounds - using middle of image')
                    end
                    if runtype == 3
                        if coordanswer(3) >= coordanswer(1) && coordanswer(3) <= imagesize(2)
                            endcolumn = coordanswer(3);
                        else
                            endcolumn = (5-1)+startcolumn;
                            warning('Specified coordinate out of bounds - using small area')
                        end
                        if coordanswer(4) >= coordanswer(2) && coordanswer(4) <= imagesize(1)
                            endrow = coordanswer(4);
                        else
                            endrow = (5-1)+startrow; 
                            warning('Specified coordinate out of bounds - using small area')
                        end
                    else
                        endrow = (2-1)+startrow;
                        endcolumn = (1-1)+startcolumn;
                    end
                else
                    startrow = 1; startcolumn = 1;
                    endrow = imagesize(1); endcolumn = imagesize(2);
                end

                if currdatatype == 0 % for solution data - no pixels to loop over
                    startrow = 1; startcolumn = 1; endrow = 1; endcolumn = 1;
                else % for image data - make arrays to write outputs into
                    r2array = zeros([imagesize(1) imagesize(2)]);
                    lt1array = r2array; lt2array = r2array;
                    bfarray = r2array; fwhmarray = r2array;
                    a1a2array = r2array;
                end

                % Loop over pixels in file
                for pixelrow=startrow:endrow
                    if showprogress == 1 && currdatatype == 1 % show progress in image
                        disp('Processing Row '+string(pixelrow-startrow+1)+' of '+string(endrow-startrow+1)+...
                            ' in '+filename);
                    end
                    for pixelcol=startcolumn:endcolumn
                        if currdatatype == 1 % load pixel data
                            signal = squeeze(data(pixelrow,pixelcol,:));
                            sigsds = zeros(height(signal),width(signal));
                        else % average together pixels at a given T value
                            if ~contains(currentfile,'.txt')
                                signal = squeeze(mean(data, [1 2]));
                                sigsds = squeeze(std(data, 0, [1 2]));
                            end
                        end
                        if troubleshoot == 1
                            figure; errorbar(delaypos,signal,sigsds,'o'); title('Signal vs position');
                        end
                
                        % Converting to time
                        cmmps = 299792458*1E3*1E-12; % c in mm/ps
                        t = zeros([length(delaypos),1]);
                        if isempty(currt0) % guess t0 by peak
                            indexoffset = floor(length(delaypos)/10); % using an offset to ignore sharp decay at start
                            [sigmax, t0index] = max(signal(indexoffset:end)); % estimating t0 by peak
                            t(:) = (delaypos(:)-delaypos(t0index+indexoffset-1))*4/cmmps; % time vector in ps
                        else % defined t0 position
                            t(:) = (delaypos(:)-currt0)*4/cmmps; % time vector in ps
                        end

                        % Power normalization
                        rawsig = signal;
                        [signal,sigsds] = powernorm(signal,sigsds,IRpower,...
                            prbpower,prbWL,pinholeyn,ND,tempmod,pmtgain,powernormtype);

                        % Baseline fitting
                        [basecurve,bleachrate,basefit,lbbase,ubbase,sbr,corrsig,...
                            corrsigbase,tbase,sigbase,currbasefittype] = ...
                            basefitfn(t,signal,currbasefittype,cutlow,cuthigh,troubleshoot);

                        % Calculate summary values
                        [sigpeak,peaksign,peaksd,integsig,integsd,noise,snr,basesd,snrv2,snrsweep,intsnr] = ...
                            summaryvals(corrsig,corrsigbase,sigsds);

                        % Lifetime fitting
                        if max(ismember(fitchannels,2)) && snr >= snrcutoff % if specified and SNR is high enough
                            [tint,sigint,fitvector,tfit,fitcurve,fitval,...
                                lifetime1,lifetime2,a1a2,fwhm,r2,resid,srr,fftx,ffty,currltfittype] = ...
                                ltfitfn(t,signal,currltfittype,peaksign,sigpeak,snr,...
                                currpulsewidth,ltmin,ltmax,basefit,currbasefittype,floatbase,IRWN,prbWL,troubleshoot);
                        else % don't do lifetime fitting
                            [tint,sigint,fitvector,tfit,fitcurve,fitval,...
                                lifetime1,lifetime2,a1a2,fwhm,r2,resid,srr,fftx,ffty,~] = ...
                                ltfitfn(t,signal,0,peaksign,sigpeak,snr,...
                                currpulsewidth,ltmin,ltmax,basefit,currbasefittype,floatbase,IRWN,prbWL,troubleshoot);
                        end

                        % Process CH1
                        if max(ismember(channels,1))
                            if currdatatype == 1 % load pixel data
                                ch1 = squeeze(ch1data(pixelrow,pixelcol,:));
                                ch1sds = zeros(height(ch1),width(ch1));
                            else
                                ch1 = squeeze(mean(ch1data, [1 2]));
                                ch1sds = squeeze(std(ch1data, 0, [1 2]));
                            end
                            % Power normalization
                            rawch1 = ch1;
                            [ch1,ch1sds] = powernorm(ch1,ch1sds,IRpower,...
                                prbpower,prbWL,pinholeyn,ND,tempmod,pmtgain,powernormtype);
        
                            % Baseline fitting
                            [ch1basecurve,ch1bleachrate,ch1basefit,ch1lbbase,ch1ubbase,ch1sbr,corrch1,...
                                corrch1base,~,ch1base,currbasefittype] = ...
                                basefitfn(t,ch1,currbasefittype,cutlow,cuthigh,troubleshoot);
        
                            % Calculate summary values
                            [ch1peak,ch1peaksign,ch1peaksd,ch1integsig,ch1integsd,ch1noise,ch1snr,ch1basesd,ch1snrv2,ch1snrsweep,ch1intsnr] = ...
                                summaryvals(corrch1,corrch1base,ch1sds);

                            if max(ismember(fitchannels,1)) && ch1snr >= snrcutoff
                                % Lifetime fitting
                                [tint,ch1int,ch1fitvector,~,ch1fitcurve,ch1fitval,...
                                    ch1lifetime1,ch1lifetime2,ch1a1a2,ch1fwhm,ch1r2,ch1resid,ch1srr,ch1fftx,ch1ffty,currltfittype] = ...
                                    ltfitfn(t,ch1,currltfittype,ch1peaksign,ch1peak,ch1snr,...
                                    currpulsewidth,ltmin,ltmax,ch1basefit,currbasefittype,floatbase,IRWN,prbWL,troubleshoot);
                            else % run with currltfittype = 0
                                [tint,ch1int,ch1fitvector,~,ch1fitcurve,ch1fitval,...
                                    ch1lifetime1,ch1lifetime2,ch1a1a2,ch1fwhm,ch1r2,ch1resid,ch1srr,ch1fftx,ch1ffty,~] = ...
                                    ltfitfn(t,ch1,0,ch1peaksign,ch1peak,ch1snr,...
                                    currpulsewidth,ltmin,ltmax,ch1basefit,currbasefittype,floatbase,IRWN,prbWL,troubleshoot);
                            end
                        end
                        % Process CH3
                        if max(ismember(channels,3))
                            if currdatatype == 1 % load pixel data
                                ch3 = squeeze(ch3data(pixelrow,pixelcol,:));
                                ch3sds = zeros(height(ch3),width(ch3));
                            else
                                ch3 = squeeze(mean(ch3data, [1 2]));
                                ch3sds = squeeze(std(ch3data, 0, [1 2]));
                            end
                            % Power normalization
                            rawch3 = ch3;
                            [ch3,ch3sds] = powernorm(ch3,ch3sds,IRpower,...
                                prbpower,prbWL,pinholeyn,ND,tempmod,pmtgain,powernormtype);
        
                            % Baseline fitting
                            [ch3basecurve,ch3bleachrate,ch3basefit,ch3lbbase,ch3ubbase,ch3sbr,corrch3,...
                                corrch3base,~,ch3base,currbasefittype] = ...
                                basefitfn(t,ch3,currbasefittype,cutlow,cuthigh,troubleshoot);
        
                            % Calculate summary values
                            [ch3peak,ch3peaksign,ch3peaksd,ch3integsig,ch3integsd,ch3noise,ch3snr,ch3basesd,ch3snrv2,ch3snrsweep,ch3intsnr] = ...
                                summaryvals(corrch3,corrch3base,ch3sds);
        
                            if max(ismember(fitchannels,3)) && ch3snr >= snrcutoff
                                % Lifetime fitting
                                [tint,ch3int,ch3fitvector,~,ch3fitcurve,ch3fitval,...
                                    ch3lifetime1,ch3lifetime2,ch3a1a2,ch3fwhm,ch3r2,ch3resid,ch3srr,ch3fftx,ch3ffty,currltfittype] = ...
                                    ltfitfn(t,ch3,currltfittype,ch3peaksign,ch3peak,ch3snr,...
                                    currpulsewidth,ltmin,ltmax,ch3basefit,currbasefittype,floatbase,IRWN,prbWL,troubleshoot);
                            else % run with currltfittype = 0
                                [tint,ch3int,ch3fitvector,~,ch3fitcurve,ch3fitval,...
                                    ch3lifetime1,ch3lifetime2,ch3a1a2,ch3fwhm,ch3r2,ch3resid,ch3srr,ch3fftx,ch3ffty,~] = ...
                                    ltfitfn(t,ch3,0,ch3peaksign,ch3peak,ch3snr,...
                                    currpulsewidth,ltmin,ltmax,ch3basefit,currbasefittype,floatbase,IRWN,prbWL,troubleshoot);
                            end
                        end
                        % Process CH4
                        if max(ismember(channels,4))
                            if currdatatype == 1 % load pixel data
                                ch4 = squeeze(ch4data(pixelrow,pixelcol,:));
                                ch4sds = zeros(height(ch4),width(ch4));
                            else
                                ch4 = squeeze(mean(ch4data, [1 2]));
                                ch4sds = squeeze(std(ch4data, 0, [1 2]));
                            end
                            % Power normalization
                            rawch4 = ch4;
                            [ch4,ch4sds] = powernorm(ch4,ch4sds,IRpower,...
                                prbpower,prbWL,pinholeyn,ND,tempmod,pmtgain,powernormtype);
        
                            % Baseline fitting
                            [ch4basecurve,ch4bleachrate,ch4basefit,ch4lbbase,ch4ubbase,ch4sbr,corrch4,...
                                corrch4base,~,ch4base,currbasefittype] = ...
                                basefitfn(t,ch4,currbasefittype,cutlow,cuthigh,troubleshoot);
        
                            % Calculate summary values
                            [ch4peak,ch4peaksign,ch4peaksd,ch4integsig,ch4integsd,ch4noise,ch4snr,ch4basesd,ch4snrv2,ch4snrsweep,ch4intsnr] = ...
                                summaryvals(corrch4,corrch4base,ch4sds);
        
                            if max(ismember(fitchannels,4)) && ch4snr >= snrcutoff
                                % Lifetime fitting
                                [tint,ch4int,ch4fitvector,~,ch4fitcurve,ch4fitval,...
                                    ch4lifetime1,ch4lifetime2,ch4a1a2,ch4fwhm,ch4r2,ch4resid,ch4srr,ch4fftx,ch4ffty,currltfittype] = ...
                                    ltfitfn(t,ch4,currltfittype,ch4peaksign,ch4peak,ch4snr,...
                                    currpulsewidth,ltmin,ltmax,ch4basefit,currbasefittype,floatbase,IRWN,prbWL,troubleshoot);
                            else % run with currltfittype = 0
                                [tint,ch4int,ch4fitvector,~,ch4fitcurve,ch4fitval,...
                                    ch4lifetime1,ch4lifetime2,ch4a1a2,ch4fwhm,ch4r2,ch4resid,ch4srr,ch4fftx,ch4ffty,~] = ...
                                    ltfitfn(t,ch4,0,ch4peaksign,ch4peak,ch4snr,...
                                    currpulsewidth,ltmin,ltmax,ch4basefit,currbasefittype,floatbase,IRWN,prbWL,troubleshoot);
                            end
                        end
                        % Process SPCM
                        if max(ismember(channels,5))
                            if currdatatype == 1 % load pixel data
                                spcm = squeeze(spcmdata(pixelrow,pixelcol,:));
                                spcmsds = zeros(height(spcm),width(spcm));
                            else
                                spcm = squeeze(mean(spcmdata, [1 2]));
                                spcmsds = squeeze(std(spcmdata, 0, [1 2]));
                            end
                            % Power normalization
                            rawspcm = spcm;
                            [spcm,spcmsds] = powernorm(spcm,spcmsds,IRpower,...
                                prbpower,prbWL,pinholeyn,ND,tempmod,pmtgain,powernormtype);
        
                            % Baseline fitting
                            [spcmbasecurve,spcmbleachrate,spcmbasefit,spcmlbbase,spcmubbase,spcmsbr,corrspcm,...
                                corrspcmbase,~,spcmbase,currbasefittype] = ...
                                basefitfn(t,spcm,currbasefittype,cutlow,cuthigh,troubleshoot);
        
                            % Calculate summary values
                            [spcmpeak,spcmpeaksign,spcmpeaksd,spcmintegsig,spcmintegsd,spcmnoise,spcmsnr,spcmbasesd,spcmsnrv2,spcmsnrsweep,spcmintsnr] = ...
                                summaryvals(corrspcm,corrspcmbase,spcmsds);
        
                            if max(ismember(fitchannels,5)) && spcmsnr >= snrcutoff
                                % Lifetime fitting
                                [tint,spcmint,spcmfitvector,~,spcmfitcurve,spcmfitval,...
                                    spcmlifetime1,spcmlifetime2,spcma1a2,spcmfwhm,spcmr2,spcmresid,spcmsrr,spcmfftx,spcmffty,currltfittype] = ...
                                    ltfitfn(t,spcm,currltfittype,spcmpeaksign,spcmpeak,spcmsnr,...
                                    currpulsewidth,ltmin,ltmax,spcmbasefit,currbasefittype,floatbase,IRWN,prbWL,troubleshoot);
                            else % run with currltfittype = 0
                                [tint,spcmint,spcmfitvector,~,spcmfitcurve,spcmfitval,...
                                    spcmlifetime1,spcmlifetime2,spcma1a2,spcmfwhm,spcmr2,spcmresid,spcmsrr,spcmfftx,spcmffty,~] = ...
                                    ltfitfn(t,spcm,0,spcmpeaksign,spcmpeak,spcmsnr,...
                                    currpulsewidth,ltmin,ltmax,spcmbasefit,currbasefittype,floatbase,IRWN,prbWL,troubleshoot);
                            end
                        end

                        if prbWL < 961 && prbWL > 699 && ~isempty(currt0) % interpolate for tD compensation
                            if isempty(tuningrate)
                                % if date before service visit
                                    %tuningrate = 0.033 % ps/nm
                                % else
                                tuningrate = 0.026;
                            end
                            alignlength = ceil(((t(end)-t(1))/(tuningrate/2))+1); % want talignspacing >= tuningrate/2
                            talign = linspace(t(1),t(end),alignlength).';
                            sigalign = spline(t, corrsig, talign);
                            talignspacing = talign(2)-talign(1); % evenly spaced

                            % tD compensation
                            if isempty(minprobe)
                                minprobe = 750;
                            end
                            if isempty(maxprobe)
                                maxprobe = 960;
                            end
                            if prbWL < minprobe
                                minprobe = prbWL;
                            end
                            fulltimeoffset = tuningrate*(maxprobe-minprobe); % ps
                            currentprobetimeoffset = tuningrate*(prbWL-minprobe); % ps
                            alignstart = floor(currentprobetimeoffset/talignspacing); % index
                            maxalignstart = floor(fulltimeoffset/talignspacing);
                            alignend = length(talign)-(maxalignstart-alignstart);
                            alignstart = alignstart+1;
                            talign = talign(alignstart:alignend);
                            sigalign = sigalign(alignstart:alignend);
                        else
                            talign = t; sigalign = corrsig;
                        end
                        if IRWN > 2500 && IRWN < 3500 && ~isempty(currt0) % interpolate for tD compensation
                            idlertuningrate = 0.0069; % ps/cm-1
                            alignlength = ceil(((t(end)-t(1))/(idlertuningrate/2))+1); % want talignspacing >= tuningrate/2
                            talign = linspace(t(1),t(end),alignlength).';
                            sigalign = spline(t, corrsig, talign);
                            talignspacing = talign(2)-talign(1); % evenly spaced

                            % tD compensation
                            minidler = 2500; maxidler = 3500;
                            fullidleroffset = idlertuningrate*(maxidler-minidler); % ps
                            currentidlertimeoffset = idlertuningrate*(IRWN-minidler); % ps
                            alignstart = floor(currentidlertimeoffset/talignspacing); % index
                            maxalignstart = floor(fullidleroffset/talignspacing);
                            alignend = length(talign)-(maxalignstart-alignstart);
                            talign = talign(alignstart+1:alignend);
                            sigalign = sigalign(alignstart+1:alignend);
                        end

                        if currdatatype == 1 % for images
                            bfarray(pixelrow,pixelcol) = sigpeak;
                            lt1array(pixelrow,pixelcol) = lifetime1;
                            lt2array(pixelrow,pixelcol) = lifetime2;
                            a1a2array(pixelrow,pixelcol) = a1a2;
                            fwhmarray(pixelrow,pixelcol) = fwhm;
                            r2array(pixelrow,pixelcol) = r2;
                        end

                        % Prepare output filename
                        outname = string(currentfile(1:end-8));

                        if currdatatype ~= 1 || runtype == 1 || runtype == 2 % for all solution data or test runs on partial images
                            % Make rounded strings for figure annotation
                            lt1 = string(sprintf('%0.3g',lifetime1)); % round to 3 sig figs
                            lt2 = string(sprintf('%0.3g',lifetime2));
                            sr2 = string(sprintf('%0.5g',r2));
                            slt1lt2 = string(sprintf('%0.3g',a1a2));
                            pulsewidth = string(sprintf('%0.3g',fwhm));
                            
                            if verbose > 0 % make figure annotation
                                annot = {'ω_{IR} = '+string(IRWN)+' cm^{-1}, λ_{probe} = '+...
                                    string(prbWL)+' nm, SNR = '+string(snr)};
                                if ~isempty(currltfittype)
                                    if currltfittype == 1
                                        annot(length(annot)+1) = {'SRR = '+string(srr)+', τ_{mono} = '+lt1+...
                                            ' ps, τ_{p} = '+pulsewidth+' ps, r^2 = '+sr2};
                                    end
                                    if currltfittype == 2
                                        annot(length(annot)+1) = {'SRR = '+string(srr)+', τ_1 = '+lt1+' ps, τ_2 = '+...
                                            lt2+' ps, τ_{p} = '+pulsewidth+...
                                            ' ps, A_{1}/A_{2} = '+slt1lt2+', r^2 = '+sr2};
                                    end
                                    if currltfittype == 3
                                        annot(length(annot)+1) = {'SRR = '+string(srr)+', τ_{mono} = '+lt1+...
                                            ' ps, τ_{p} = '+pulsewidth+' ps, β = '+...
                                            string(fitval(12))+', r^2 = '+sr2};
                                    end
                                    if currltfittype == 4
                                        annot(length(annot)+1) = {'SRR = '+string(srr)+', τ_1 = '+lt1+' ps, τ_2 = '+...
                                            lt2+' ps, τ_{p} = '+pulsewidth+' ps, β = '+string(fitval(12))+...
                                            ', A_{1}/A_{2} = '+slt1lt2+', r^2 = '+sr2};
                                    end
                                    if currltfittype > 0 && r2 < 0.7
                                        annot(end) = {'Poor fit, r^2 = ' + sr2};
                                    end
                                end
                                if verbose == 2 % print out all experimental parameters available
                                    if pinholeyn == 1
                                        phlabel = 'pinhole in place';
                                    else
                                        phlabel = 'bypassed pinhole';
                                    end
                                    annot(length(annot)+1) = {'PMT gain '+string(pmtgain)+...
                                        ', PMT BW '+string(pmtBW/1e3)+' kHz, '+...
                                        phlabel+', DM'+string(prbdichroic)+...
                                        ', PMT BP'+string(pmtfilter)};
                                end
                            else
                                annot = {};
                            end
    
                            % Prepare labels for figure
                            basestring = strings(4,1);
                            basestring(1) = 'Baseline'; basestring(2) = 'Linear baseline';
                            basestring(3) = 'Exponential baseline'; basestring(4) = 'Exp + linear baseline';
                            baselegend = {'Data','Baseline data',basestring(currbasefittype+1)}; ch1legend = baselegend;
                            ch2legend = baselegend; ch3legend = baselegend; ch4legend = baselegend; spcmlegend = baselegend;
                            if powernormtype > 0 % prepare y-axis label
                                normlabel = 'Norm. ';
                                if powernormtype == 1
                                    pmtunitlabel = '(V/mW_{IR})';
                                    spcmunitlabel = '(cpms/mW_{IR})';
                                else
                                    if powernormtype == 2
                                        pmtunitlabel = '(V/mW^{2})';
                                    else % powernormtype == 3
                                        pmtunitlabel = '(V_{gc}/mW^{2})';
                                    end
                                    spcmunitlabel = '(cpms/mW^{2})';
                                end
                            else
                                normlabel = 'Raw ';
                                pmtunitlabel = '(V)';
                                spcmunitlabel = '(cpms)';
                            end
    
                            ch1label = string(normlabel)+'PMT DC '+string(pmtunitlabel);
                            ch2label = string(normlabel)+'PMT AC '+string(pmtunitlabel);
                            ch3label = string(normlabel)+'SPCM DC '+string(pmtunitlabel);
                            ch4label = string(normlabel)+'SPCM AC '+string(pmtunitlabel);
                            spcmlabel = string(normlabel)+'SPCM '+string(spcmunitlabel);

                            % Figure layouts
                            %%% Type 1 - DC (CH1/SPCM) and AC (CH2)
                            % | CH1 CH1 | Info|
                            % | CH2  CH2  CH2 |
                            % | CH2  CH2  CH2 |
                            %%% Type 2 - CH1-CH4+SPCM
                            % | CH1 CH1 | Info |  SPCM  |
                            % | CH2  CH2  CH2 | CH3 CH3 |
                            % | CH2  CH2  CH2 | CH4 CH4 |
                            %%% Type 3 - CH2 only
                            % |   Info   |
                            % | CH2  CH2 |
                            % | CH2  CH2 |

                            % Make figure
                            fig = figure('visible',visibility);
                            if length(nonzeros(channels)) == 2 % DC (CH1/SPCM) and AC (CH2)
                                tiledlayout(3,3); nexttile([1 2]);  % DC
                                if max(ismember(channels,1))
                                    makesubpanel(t,ch1,ch1sds,tbase,ch1base,ch1basecurve,tfit,ch1fitcurve,...
                                        tint,ch1resid,ch1label,ch1legend,1);
                                else % SPCM
                                    makesubpanel(t,spcm,spcmsds,tbase,spcmbase,spcmbasecurve,tfit,spcmfitcurve,...
                                        tint,spcmresid,spcmlabel,spcmlegend,5);
                                end
                                nexttile([1 1]); xticks([]); yticks([]);
                                printoutdim = [0.67 0.7 0.3 0.25]; nexttile([2 3]);
                            else
                                % CH1-4 and SPCM
                                if length(nonzeros(channels)) > 2
                                    tiledlayout(3,5); nexttile([1 2]); % CH1
                                    makesubpanel(t,ch1,ch1sds,tbase,ch1base,ch1basecurve,tfit,ch1fitcurve,...
                                        tint,ch1resid,ch1label,ch1legend,1);
                                    nexttile([1 1]); xticks([]); yticks([]);
                                    printoutdim = [0.4 0.7 0.2 0.25];
                                    nexttile([1 2]); % SPCM
                                    makesubpanel(t,spcm,spcmsds,tbase,spcmbase,spcmbasecurve,tfit,spcmfitcurve,...
                                        tint,spcmresid,spcmlabel,spcmlegend,5);
                                    nexttile(9,[1 2]); % CH3
                                    makesubpanel(t,ch3,ch3sds,tbase,ch3base,ch3basecurve,tfit,ch3fitcurve,...
                                        tint,ch3resid,ch3label,ch3legend,3);
                                    nexttile(14,[1 2]); % CH4
                                    makesubpanel(t,ch4,ch4sds,tbase,ch4base,ch4basecurve,tfit,ch4fitcurve,...
                                        tint,ch4resid,ch4label,ch4legend,4);
                                    nexttile(6,[2 3]);
                                else % just CH2
                                    tiledlayout(3,2); nexttile([1 2]); xticks([]); yticks([]);
                                    printoutdim = [0.05 0.7 0.9 0.25]; nexttile([2 2]); 
                                end
                            end
                            % Plot CH2
                            makesubpanel(t,signal,sigsds,tbase,sigbase,basecurve,tfit,fitcurve,...
                                tint,resid,ch2label,ch2legend,2);
                            % Annotate with experiment info
                            annotation('rectangle',printoutdim,'Color',[1 1 1],'FaceColor',[1 1 1]); % box to cover
                            annotation('textbox',printoutdim,'String',annot);
                            if writefigsyn == 1
                                saveas(fig,outname+'_proc.png')
                            end
                        end
                    end
                end

                % Make lifetime images
                if currltfittype > 0 && currdatatype == 1 && runtype ~= 1 && runtype ~= 2
                    % Make blue-red colormap
                    startcolor = [0 0 1]; endcolor = [1 0 0];
                    map1 = [linspace(startcolor(1),1,100) linspace(1,endcolor(1),100)];
                    map2 = [linspace(startcolor(2),1,100) linspace(1,endcolor(2),100)];
                    map3 = [linspace(startcolor(3),1,100) linspace(1,endcolor(3),100)];
                    redbluemap = [map1.' map2.' map3.'];
                    % Make image
                    fitimage = figure('visible',visibility);
                    tiledlayout(1+(1-(currltfittype == 1)),2); % 1x2 monoexp, 2x2 biexp
                    nexttile; tile1 = nexttile(1); % BonFIRE
                    makeimagepanel(bfarray,xsteps,ysteps,1); 
                    colormap(tile1,'hot');
                    nexttile; tile2 = nexttile(2); % τ1
                    makeimagepanel(lt1array,xsteps,ysteps,2);
                    colormap(tile2,redbluemap);
                    if currltfittype ~= 1 % show biexponential images
                        nexttile; tile3 = nexttile(3); % τ2
                        makeimagepanel(lt2array,xsteps,ysteps,3);
                        colormap(tile3,'winter');
                        nexttile; tile4 = nexttile(4); % A1/A2
                        makeimagepanel(a1a2array,xsteps,ysteps,4);
                        colormap(tile4,'spring');
                    end
                    if writeprocyn == 1
                        saveas(fitimage,outname+'_fit.png','png');
                    end
                end
                if writeprocyn == 1 && currdatatype == 1 % write processed image tif
                    Save_tif(outname+'_proc.tif',bfarray,lt1array,...
                        lt2array,a1a2array,r2array)
                end
                if runtype == 0 % close background figures if not a test run
                    close all;
                end

                % Make structure for current file
                cs.fileName = filename;
                % Write experimental parameters to current structure
                cs.parameters.IRWN = IRWN; cs.parameters.IRpower = IRpower; cs.parameters.dutycycle = dutycycle;
                cs.parameters.prbWL = prbWL; cs.parameters.prbpower = prbpower; cs.parameters.ND = ND; cs.parameters.pinholeyn = pinholeyn;
                cs.parameters.probeFWHM = probeFWHM; cs.parameters.levanteFWHM = levanteFWHM;
                cs.parameters.modfreq = modfreq; cs.parameters.pmtgain = pmtgain; cs.parameters.pmtbandwidth = pmtBW;
                cs.parameters.pmtfilter = pmtfilter; cs.parameters.spcmfilter = spcmfilter; cs.parameters.cmosfilter = cmosfilter;
                cs.parameters.xsteps = xsteps; cs.parameters.ysteps = ysteps; cs.parameters.zsteps = zsteps; cs.parameters.tsteps = tsteps;
                cs.parameters.xinitial = xinitial; cs.parameters.yinitial = yinitial; cs.parameters.zinitial = zinitial; cs.parameters.tinitial = tinitial;
                cs.parameters.xfinal = xfinal; cs.parameters.yfinal = yfinal; cs.parameters.zfinal = zfinal; cs.parameters.tfinal = tfinal; cs.parameters.dwelltime = dwelltime;
                % Write data, summary values, and fit values to current structure
                if currdatatype == 0 % solution info
                    cs.data.t = t; cs.data.signal = signal; cs.data.tint = tint; cs.data.sigint = sigint;
                    cs.data.rawsig = rawsig; cs.data.corrsig = corrsig; cs.data.talign = talign;
                    cs.data.sigalign = sigalign; cs.data.tbase = tbase; cs.data.sigbase = sigbase;
                    cs.data.basecurve = basecurve; cs.data.snrsweep = snrsweep; cs.data.delaypos = delaypos;
                    cs.summary.peakheight = sigpeak; cs.summary.peakarea = integsig;
                    cs.summary.snr = snr; cs.summary.snrv2 = snrv2; cs.summary.sbr = sbr;
                    cs.summary.noise = noise; cs.summary.basesd = basesd; cs.summary.nfiles = totalfiles;
                    cs.fit.lifetime1 = lifetime1; cs.fit.lifetime2 = lifetime2;
                    cs.fit.a1a2 = a1a2; cs.fit.r2 = r2; cs.fit.resid = resid; cs.fit.fitval = fitval;
                    cs.fit.pulsewidth = fwhm; cs.fit.srr = srr; cs.fit.fitvector = fitvector;
                    cs.fit.tfit = tfit; cs.fit.fitcurve = fitcurve; cs.fit.fftx = fftx; cs.fit.ffty = ffty;
                    if max(ismember(channels,1))
                        cs.ch1data.signal = ch1; cs.ch1data.sigint = ch1int; cs.ch1data.rawsig = rawch1; cs.ch1data.corrsig = corrch1;
                        cs.ch1data.sigbase = ch1base; cs.ch1data.basecurve = ch1basecurve; cs.ch1data.snrsweep = ch1snrsweep;
                        cs.ch1summary.peakheight = ch1peak; cs.ch1summary.peakarea = ch1integsig;
                        cs.ch1summary.snr = ch1snr; cs.ch1summary.snrv2 = ch1snrv2; cs.ch1summary.sbr = ch1sbr;
                        cs.ch1summary.noise = ch1noise; cs.ch1summary.basesd = ch1basesd;
                        cs.ch1fit.lifetime1 = ch1lifetime1; cs.ch1fit.lifetime2 = ch1lifetime2;
                        cs.ch1fit.a1a2 = ch1a1a2; cs.ch1fit.r2 = ch1r2; cs.ch1fit.resid = ch1resid; cs.ch1fit.fitval = ch1fitval;
                        cs.ch1fit.pulsewidth = ch1fwhm; cs.ch1fit.srr = ch1srr; cs.ch1fit.fitvector = ch1fitvector;
                        cs.ch1fit.fitcurve = ch1fitcurve; cs.ch1fit.fftx = ch1fftx; cs.ch1fit.ffty = ch1ffty;
                    end
                    if max(ismember(channels,3))
                        cs.ch3data.signal = ch3; cs.ch3data.sigint = ch3int; cs.ch3data.rawsig = rawch3; cs.ch3data.corrsig = corrch3;
                        cs.ch3data.sigbase = ch3base; cs.ch3data.basecurve = ch3basecurve; cs.ch3data.snrsweep = ch3snrsweep;
                        cs.ch3summary.peakheight = ch3peak; cs.ch3summary.peakarea = ch3integsig;
                        cs.ch3summary.snr = ch3snr; cs.ch3summary.snrv2 = ch3snrv2; cs.ch3summary.sbr = ch3sbr;
                        cs.ch3summary.noise = ch3noise; cs.ch3summary.basesd = ch3basesd;
                        cs.ch3fit.lifetime1 = ch3lifetime1; cs.ch3fit.lifetime2 = ch3lifetime2;
                        cs.ch3fit.a1a2 = ch3a1a2; cs.ch3fit.r2 = ch3r2; cs.ch3fit.resid = ch3resid; cs.ch3fit.fitval = ch3fitval;
                        cs.ch3fit.pulsewidth = ch3fwhm; cs.ch3fit.srr = ch3srr; cs.ch3fit.fitvector = ch3fitvector;
                        cs.ch3fit.fitcurve = ch3fitcurve; cs.ch3fit.fftx = ch3fftx; cs.ch3fit.ffty = ch3ffty;
                    end
                    if max(ismember(channels,4))
                        cs.ch4data.signal = ch4; cs.ch4data.sigint = ch4int; cs.ch4data.rawsig = rawch4; cs.ch4data.corrsig = corrch4;
                        cs.ch4data.sigbase = ch4base; cs.ch4data.basecurve = ch4basecurve; cs.ch4data.snrsweep = ch4snrsweep;
                        cs.ch4summary.peakheight = ch4peak; cs.ch4summary.peakarea = ch4integsig;
                        cs.ch4summary.snr = ch4snr; cs.ch4summary.snrv2 = ch4snrv2; cs.ch4summary.sbr = ch4sbr;
                        cs.ch4summary.noise = ch4noise; cs.ch4summary.basesd = ch4basesd;
                        cs.ch4fit.lifetime1 = ch4lifetime1; cs.ch4fit.lifetime2 = ch4lifetime2;
                        cs.ch4fit.a1a2 = ch4a1a2; cs.ch4fit.r2 = ch4r2; cs.ch4fit.resid = ch4resid; cs.ch4fit.fitval = ch4fitval;
                        cs.ch4fit.pulsewidth = ch4fwhm; cs.ch4fit.srr = ch4srr; cs.ch4fit.fitvector = ch4fitvector;
                        cs.ch4fit.fitcurve = ch4fitcurve; cs.ch4fit.fftx = ch4fftx; cs.ch4fit.ffty = ch4ffty;
                    end
                    if max(ismember(channels,5))
                        cs.spcmdata.signal = spcm; cs.spcmdata.sigint = spcmint; cs.spcmdata.rawsig = rawspcm; cs.spcmdata.corrsig = corrspcm;
                        cs.spcmdata.sigbase = spcmbase; cs.spcmdata.basecurve = spcmbasecurve; cs.spcmdata.snrsweep = spcmsnrsweep;
                        cs.spcmsummary.peakheight = spcmpeak; cs.spcmsummary.peakarea = spcmintegsig;
                        cs.spcmsummary.snr = spcmsnr; cs.spcmsummary.snrv2 = spcmsnrv2; cs.spcmsummary.sbr = spcmsbr;
                        cs.spcmsummary.noise = spcmnoise; cs.spcmsummary.basesd = spcmbasesd;
                        cs.spcmfit.lifetime1 = spcmlifetime1; cs.spcmfit.lifetime2 = spcmlifetime2;
                        cs.spcmfit.a1a2 = spcma1a2; cs.spcmfit.r2 = spcmr2; cs.spcmfit.resid = spcmresid; cs.spcmfit.fitval = spcmfitval;
                        cs.spcmfit.pulsewidth = spcmfwhm; cs.spcmfit.srr = spcmsrr; cs.spcmfit.fitvector = spcmfitvector;
                        cs.spcmfit.fitcurve = spcmfitcurve; cs.spcmfit.fftx = spcmfftx; cs.spcmfit.ffty = spcmffty;
                    end
                else % image info - reshaped into vectors to allow saving structure
                    cs.fit.lt1vector = nonzeros(lt1array); cs.fit.lt2vector = nonzeros(lt2array); cs.fit.a1a2vector = nonzeros(lt2array);
                end
                % Write values to summary structure
                summary.('folder'+string(ii)).('file'+string(jj)) = cs;
                % Clear cs for next file
                cs = structfun(@(x) [], cs, 'UniformOutput', false);
            end
        end
    end
    % Write out final results
    if writeprocyn == 1 % write summary file
        writestruct(summary,'summarystructure.xml',FileType="xml");
    end
else % reload previously processed data
    if loadprevious == 2
        summary = readstruct('summarystructure.xml','FileType','xml');
    else
        addinfo = importdata('add_info.yml'); % load info from previous analysis
        master_size = addinfo(1:4);
        master = NaN(master_size(1),master_size(2),master_size(3),master_size(4));
        targetfolders = addinfo(5:end);
        emptydatcount = 0;
        folders = targetfolders;
        for ii = folders
            datafilesraw = dir(fullfile(workingdirectory,subfolders{ii},'*_proc.dat'));
            datafiles = {datafilesraw(~[datafilesraw.isdir]).name}; % all proc files
            if ~isempty(datafilesraw)
                for jj = 1:numel(datafiles)
                    currentfile = fullfile(workingdirectory,subfolders{ii},datafiles{jj});
                    tempimport = importdata(currentfile);
                    master(1:height(tempimport),1:width(tempimport),ii,jj) = tempimport;
                end
            else
                emptydatcount = emptydatcount+1;
            end
        end
        if emptydatcount == length(folders) % no proc.dat files found
            checkdat = dir('batch.dat'); % try to load from batch.dat
            warning('Loading from batch.dat - master array may have errors.');
            if ~isempty(checkdat) % found batch.dat
                batchmaster = importdata('batch.dat');
                if width(batchmaster) == master_size(2)*master_size(3)*master_size(4)
                    for ii=1:master_size(3)
                        for jj=1:master_size(4)
                            startcol = 1+master_size(2)*(ii+2*jj-3);
                            endcol = startcol+master_size(2)-1;
                            master(:,:,ii,jj) = batchmaster(1:master_size(1),startcol:endcol);
                        end
                    end
                else
                    error('Batch.dat size is mismatched with addinfo.yml - please rerun analysis.');
                end
            else
                error('No batch processed files found - please rerun analysis.');
            end
        end
    end
end


%% Post-batch analysis
% Master array has been replaced by the summary structure! Use the GUI to 
% inspect its contents and summary.(folder).(file).(value) to pull from it.
% Legacy post-batch analysis (pulling from the master array) is available 
% in sweep_proc_align_v14.m

clc; close all;
clearvars -except targetfolders summary workingdirectory subfolders
figvis = 'on';

% Analysis options
analysistype = [1 2]; % 0 = nothing, 1 = load data, 2 = contours, 3 = lifetimes, [] = dialog
savecontours = 0; % 1 = save contours, 0 = not
xaxischoice = 3; % 0 = probeWL, 1 = sumfreq, 2 = detuning (specify),
% 3 = detuning (Rh800 max), 4 = detuning (Rh800 0-0)
logcontour = 0; % 0 = linear scale, 1 = log scale
subset = [2 3 4 5 6 8 9]; % [] = dialog, targetfolders = run all

if isempty(analysistype) % dialog to select analysis type
    [indx,~] = listdlg('PromptString',{'Select post-batch analysis type.'},...
    'SelectionMode','multiple','ListString',{'Load data','Contour maps','Lifetime comparisons'});
    analysistype = indx;
end

if isempty(subset) % dialog to select folders
    [indx,~] = listdlg('PromptString',{'Select folders to process.'},...
    'SelectionMode','multiple','ListString',subfolders(targetfolders));
    subset = indx;
end

sf = fieldnames(summary); % summary fields
maxfiles = 0; maxTlength = 0; maxalignlength = 0;
if ismember(analysistype,1)
    for i=1:length(subset)
        nfiles = summary.(sf{i}).nfiles;
        if maxfiles < nfiles
            maxfiles = nfiles;
        end
        for j=1:nfiles
            Tlength = length(summary.(sf{i}).('file'+string(j)).data.t);
            if maxTlength < Tlength
                maxTlength = Tlength;
            end
            alignlength = length(summary.(sf{i}).('file'+string(j)).data.talign);
            if maxalignlength < alignlength
                maxalignlength = alignlength;
            end
        end
    end
end

% Set up arrays
time = NaN(length(subset),maxTlength);
timealign = NaN(length(subset),maxalignlength);
prb = NaN(maxfiles,length(subset));
wIR = prb; lt1 = prb; lt2 = prb; ltr = prb; pulsewidth = prb;
IRpowers = prb; prbpowers = prb; noise = prb; dcsbr = prb;
pkht = NaN(maxfiles,length(subset));
csig = NaN(height(wIR),width(time),length(subset)); rawsig = csig;
csiga = NaN(height(prb),width(timealign),length(subset));

% Pulling data
for i=subset % folder
    nfiles = summary.(sf{i}).nfiles;
    for j=1:nfiles % file
        tlength = length(summary.(sf{i}).('file'+string(j)).data.t);
        alength = length(summary.(sf{i}).('file'+string(j)).data.talign);
        time(i,1:tlength) = summary.(sf{i}).('file'+string(j)).data.t;
        timealign(i,1:alength) = summary.(sf{i}).('file'+string(j)).data.talign;
        prb(j,i) = summary.(sf{i}).('file'+string(j)).parameters.prbWL;
        wIR(j,i) = summary.(sf{i}).('file'+string(j)).parameters.IRWN;
        lt1(j,i) = summary.(sf{i}).('file'+string(j)).fit.lifetime1;
        lt2(j,i) = summary.(sf{i}).('file'+string(j)).fit.lifetime2;
        ltr(j,i) = summary.(sf{i}).('file'+string(j)).fit.a1a2;
        csig(j,1:tlength,i) = summary.(sf{i}).('file'+string(j)).data.corrsig;
        rawsig(j,1:tlength,i) = summary.(sf{i}).('file'+string(j)).data.rawsig;
        csiga(j,1:alength,i) = summary.(sf{i}).('file'+string(j)).data.sigalign;
        pkht(j,i) = summary.(sf{i}).('file'+string(j)).summary.peakheight;
        pulsewidth(j,i) = summary.(sf{i}).('file'+string(j)).fit.pulsewidth;
        IRpowers(j,i) = summary.(sf{i}).('file'+string(j)).parameters.IRpower;
        prbpowers(j,i) = summary.(sf{i}).('file'+string(j)).parameters.prbpower;
        noise(j,i) = summary.(sf{i}).('file'+string(j)).summary.noise;
        dcsbr(j,i) = summary.(sf{i}).('file'+string(j)).ch1summary.sbr;
    end
end

pkht = abs(pkht);
ltavg = (ltr.*lt1+lt2)./(ltr+1);
photothermal = squeeze(rawsig(:,end,:));
defaultfont = 25;
tD = time(1,:);

if max(ismember(analysistype,2))
    % Plotting contours
    for i=subset
        nfiles = summary.(sf{i}).nfiles;
        if isscalar((unique(rmmissing(nonzeros(prb(:,i)))))) % no probe tuning --> IR sweep
            xstring = 'Time delay (ps)';
            ystring = 'ω_{IR} (cm^{-1})';
            w1 = wIR(1:nfiles,i);
            w2 = prb(1:nfiles,i);
            nonswept = string(round(w2(1)))+' nm';
            sumfreq = 1e7./w2+w1;
            % Gradients for contour maps
            if min(w1) > 2000 % idler sweep
                startcolor = [0 0.35 0.5]; % turquoise
                endcolor = [0.8 0.2 0.2]; % red
            else % DFG sweep
                startcolor = [0.5 0 0.5]; % purple
                endcolor = [0.2 0.8 0.2]; % green
            end
        else % probe sweep
            xstring = 'Compensated time delay (ps)';
            ystring = 'λ_{probe} (nm)';
            w1 = prb(1:nfiles,i);
            w2 = wIR(1:nfiles,i);
            nonswept = string(round(w2(1)))+' cm^{-1}';
            sumfreq = 1e7./w1+w2;
            % Blue-white-orange gradient for contour map
            startcolor = [0 0 1]; % blue
            endcolor = [1 0.5 0]; % orange
        end
        t1 = timealign(i,:);
        sigmatrix = csiga(:,:,i);

        if ismatrix(rmmissing(w1))
            % Auto-sorting
            [t1,tind] = sort(t1);
            [w1,wind] = sort(w1);
            [x1,x2] = meshgrid(t1,w1);
            sigmatrix = sigmatrix(wind,tind);
            if logcontour == 1
                sigmatrix = log(abs(sigmatrix));
                contourstring = 'log|Corrected signal| at '+nonswept+' (AU)';
            else
                contourstring = 'Corrected signal at '+nonswept+' (AU)';
            end

            % Color mapping
            map1 = [linspace(startcolor(1),1,100) linspace(1,endcolor(1),100)];
            map2 = [linspace(startcolor(2),1,100) linspace(1,endcolor(2),100)];
            map3 = [linspace(startcolor(3),1,100) linspace(1,endcolor(3),100)];
            map = [map1.' map2.' map3.'];

            % Plotting
            contour = figure('visible',figvis);
            tiledlayout(4,4,'TileSpacing','compact','Padding','compact');
            % Time 1D cross-section
            nexttile([1 3]);
            %plot(t1,max(sigmatrix,[],1),'Color',startcolor,'LineWidth',2); ylabel('Peak (AU)');
            plot(t1,sum(sigmatrix(1:length(rmmissing(w1)),:)/length(rmmissing(w1)),1),'Color',startcolor,'LineWidth',2); ylabel('Mean (AU)');
            xlim([min(t1) max(t1)]); xlabel(xstring);

            % Normal contour
            nexttile([3 3]);
            contourf(x1,x2,sigmatrix); colormap(map); cb=colorbar;
            cb.Label.String=contourstring;
            cb.Label.Rotation=270; cb.Label.VerticalAlignment = "bottom";
            xlabel(xstring); ylabel(ystring);
            xlim([min(t1) max(t1)]); ylim([min(w1) max(w1)]);

            % Blank square (top-right)
            nexttile([1 1]);
            xticks([]); yticks([]);
            annotation('rectangle',[0.78 0.78 0.4 0.4],'Color',[1 1 1],'FaceColor',[1 1 1]); % box to cover

            % Frequency 1D cross-section
            nexttile([3 1]);
            plot(max(sigmatrix,[],2),w1,'Color',endcolor,'LineWidth',2); xlabel('Peak (AU)');
            %plot(sum(sigmatrix(:,1:length(rmmissing(t1))),2)/length(rmmissing(t1)),w1,'Color',endcolor,'LineWidth',2); xlabel('Mean (AU)');
            ylim([min(w1) max(w1)]); ylabel(ystring);
            if savecontours == 1
                saveas(contour,'Probe contours/'+nonswept+'_sweep.png')
            end
        end
    end
end

% Set up x-axis for lifetime plots
if xaxischoice == 0
    xdata = w1;
    label = 'λ_{probe} (nm)';
else
    if xaxischoice == 1
        xdata = sumfreq;
        label = 'ω_{IR}+ω_{probe} (cm^{-1})';
    else % some detuning
        if xaxischoice == 3
            detuning = sumfreq-(1e7/696);
            label = 'ω_{IR}+ω_{probe}-ω_{max} (cm^{-1})';
        else
            if xaxischoice == 4
                detuning = sumfreq-(1e7/713);
                label = 'ω_{IR}+ω_{probe}-ω_{0-0} (cm^{-1})';
            else % xaxischoice = 2 or other
                prompt = {'Specify absorption maximum (nm).'};
                dlgtitle = 'Set up detuning'; dims = [1 45]; definput = {'696'};
                coordanswer = inputdlg(prompt,dlgtitle,dims,definput);
                detuning = sumfreq-(1e7/str2double(coordanswer));
                label = 'ω_{IR}+ω_{probe}-ω_{max} (cm^{-1})';
            end
        end
        xdata = detuning;
    end
end

if max(ismember(analysistype,3))
    % Default lifetime comparison
    ltlegend = strings(length(subset),1);
    figure('visible',figvis);
    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
    % τ1
    nexttile([1 2]); hold on;
    for i=subset
        redamt = exp((i-length(subset))*2/length(subset));
        greenamt = ((-4/length(subset)^2)*(i-0.5*length(subset))^2+1)^2;
        blueamt = exp(-2*i/length(subset));
        linecolor = [redamt greenamt blueamt];
        plot(xdata,lt1(:,i),'Color',linecolor,'LineWidth',2);
        ltlegend(i) = string(wIR(1,i))+' cm{-1}';
    end
    xticks([]); ylabel('τ_{1} (ps)');
    hold off; xlim([min(xdata) max(xdata)]); ylim([0 5]); 
    legend(ltlegend);
    % τ2
    nexttile([1 2]); hold on;
    for i=1:length(subset)
        redamt = exp((i-length(subset))*2/length(subset));
        greenamt = ((-4/length(subset)^2)*(i-0.5*length(subset))^2+1)^2;
        blueamt = exp(-2*i/length(subset));
        linecolor = [redamt greenamt blueamt];
        plot(xdata,lt2(:,i),'Color',linecolor,'LineWidth',2);
    end
    xlabel('ω_{IR}+ω_{probe} (cm^{-1})'); ylabel('τ_{2} (ps)');
    hold off; xlim([min(xdata) max(xdata)]); ylim([0 10]); 
    legend(ltlegend);
    % A1/A2
    nexttile([2 1]); hold on;
    for i=1:length(subset)
        redamt = exp((i-length(subset))*2/length(subset));
        greenamt = ((-4/length(subset)^2)*(i-0.5*length(subset))^2+1)^2;
        blueamt = exp(-2*i/length(subset));
        linecolor = [redamt greenamt blueamt];
        plot(xdata,ltr(:,i),'Color',linecolor,'LineWidth',2);
    end
    set(gca,'Yscale','log');
    xlabel('ω_{IR}+ω_{probe} (cm^{-1})'); ylabel('A_{1}/A_{2}');
    xlim([min(xdata) max(xdata)]); ylim([1e-2 1e2]); 
    legend(ltlegend);
end

%% Pairwise comparisons
% somehow sort lifetimes into condensate or buffer (by BonFIRE intensity?)
condensate = [1 2 3]; buffer = [4 5 6];

% Run comparison and make summary plot
figure;
[condmean,buffermean,pval,cohend,cohenlower,cohenupper] = ...
    statcompare(condensate,buffer);

%% Lifetime histograms



%% Functions
% Most of the data processing has now been moved to dedicated functions.
% When possible, the variable names in the functions are capitalized
% versions of the variable names in the main script.
% An attempt has been made to write descriptions of each of the functions.

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

% Power normalization
function [SIGNAL,SIGSDS] = powernorm(SIGNAL,SIGSDS,IRPOWER,PRBPOWER,...
    PRBWL,PINHOLEYN,NDF,TEMPMOD,PMTGAIN,POWERNORMTYPE)
%%% This function performs power normalization, corrects for the
% power losses from the pinhole, objective, pellicle beamsplitters
% (temporal modulation), and scales for increased PMT gain.

    if POWERNORMTYPE > 0 % power normalization
        % Power corrections
        if PINHOLEYN == 1
            PRBPOWER = 0.27*PRBPOWER; % pinhole transmission is 27%
        end
        PRBPOWER = PRBPOWER/(10^(NDF)); % power loss from ND filter
        % Correct for objective transmission (400-1600 nm)
        xlpln = [84.4573697852295	84.5194986387765	84.5937140550633	84.6790481672078	84.7745331083276	84.8792010115406	84.9920840099646	85.1122142367173	85.2386238249166	85.3703449076801	85.5064096181257	85.6458500893712	85.7876984545343	85.9309868467329	86.0747473990846	86.2180122447073	86.3598135167188	86.4991833482368	86.6351538723792	86.7667572222636	86.8930255310079	87.0129909317299	87.1256855575473	87.2301415415780	87.3253910169396	87.4104661167500	87.4843989741270	87.5462217221883	87.5949717389603	87.6303181338952	87.6531553748221	87.6645184982279	87.6654425405991	87.6569625384224	87.6401135281844	87.6159305463718	87.5854486294712	87.5497028139692	87.5097281363526	87.4665596331079	87.4212323407218	87.3747812956809	87.3282415344719	87.2826480935815	87.2390360094961	87.1982925804875	87.1605878117055	87.1258743253634	87.0941046954252	87.0652314958549	87.0392073006166	87.0159846836741	86.9955162189915	86.9777544805327	86.9626520422619	86.9501614781428	86.9402353621396	86.9328262682162	86.9278867703366	86.9253694424648	86.9252268585648	86.9274163080248	86.9319601661879	86.9389282532680	86.9483914657739	86.9604207002147	86.9750868530993	86.9924608209367	87.0126135002358	87.0356157875056	87.0615385792550	87.0904527719930	87.1224292622286	87.1575389464707	87.1958527212282	87.2374414830102	87.2823761283256	87.3307261311845	87.3824558903725	87.4373544527600	87.4951942212391	87.5557475987017	87.6187869880397	87.6840847921451	87.7514134139099	87.8205452562261	87.8912527219857	87.9633082140805	88.0364841354026	88.1105528888440	88.1852868772966	88.2604585036524	88.3358401708034	88.4112042816415	88.4863352053220	88.5610665010406	88.6352443658681	88.7087149968848	88.7813245911706	88.8529193458059	88.9233454578707	88.9924491244450	89.0600765426092	89.1260739094433	89.1902874220275	89.2525632774419	89.3127476727666	89.3706868050818	89.4262268714676	89.4792140690041	89.5295018196827	89.5770247505171	89.6217685862474	89.6637198470075	89.7028650529315	89.7391907241533	89.7726833808070	89.8033295430265	89.8311157309458	89.8560284646989	89.8780542644198	89.8971796502426	89.9133911423011	89.9266752607294	89.9370185256615	89.9444074572313	89.9488291691859	89.9503001876356	89.9488793786822	89.9446288956254	89.9376108917646	89.9278875203995	89.9155209348297	89.9005732883546	89.8831067342740	89.8631834258874	89.8408655164943	89.8162151593944	89.7892945078872	89.7601657152723	89.7288909348493	89.6955323199177	89.6601520243808	89.6228191979486	89.5836273577648	89.5426753230852	89.5000619131655	89.4558859472613	89.4102462446285	89.3632416245226	89.3149709061994	89.2655329089146	89.2150264519239	89.1635503544830	89.1112034358475	89.0580845152733	89.0042924120159	88.9499259453312	88.8950839344747	88.8398680119149	88.7844058820563	88.7288393771492	88.6733104689172	88.6179611290836	88.5629333293719	88.5083690415056	88.4544102372081	88.4011988882029	88.3488769662135	88.2975864429631	88.2474692901754	88.1986674795738	88.1513229828817	88.1055777718226	88.0615738181198	88.0194509867361	87.9792746715435	87.9410176169093	87.9046467703020	87.8701290791901	87.8374314910419	87.8065209533259	87.7773644135107	87.7499288190645	87.7241811174560	87.7000882561534	87.6776171826254	87.6567348443403	87.6374081887665	87.6196041633726	87.6032897156270	87.5884317866010	87.5749891784533	87.5628966192670	87.5520844072430	87.5424828405822	87.5340222174854	87.5266328361537	87.5202449947878	87.5147889915888	87.5101951247576	87.5063936924950	87.5033149930020	87.5008893244795	87.4990469851284	87.4977182731497	87.4968334867441	87.4963229241127	87.4961165777058	87.4961420809350	87.4963259697164	87.4965947738678	87.4968750232070	87.4970932475519	87.4971759767202	87.4970497405296	87.4966410687981	87.4958764913434	87.4946825379833	87.4929857385356	87.4907126228180	87.4877897206485	87.4841435618447	87.4797006762245	87.4743886486124	87.4681628757545	87.4610086922646	87.4529129136754	87.4438623555199	87.4338438333308	87.4228441626409	87.4108501589832	87.3978486378905	87.3838264148956	87.3687703055313	87.3526671253306	87.3355036898262	87.3172668145510	87.2979433150379	87.2775200068196	87.2559837069647	87.2333218422006	87.2095233760081	87.1845775102349	87.1584734467282	87.1312003873357	87.1027475339047	87.0731040882827	87.0422592523174	87.0102022278560	86.9769222167461	86.9424084208352	86.9066500419707	86.8696362820002	86.8313563427711	86.7917994261309	86.7509547339271	86.7088138915968	86.6653842049501	86.6206792206438	86.5747125013528	86.5274976097523	86.4790481085170	86.4293775603221	86.3784995278424	86.3264275737529	86.2731752607286	86.2187561514444	86.1631838085753	86.1064717947963	86.0486336727823	85.9896830052082	85.9296333547491	85.8684991914614	85.8063134215665	85.7431261289567	85.6789880528130	85.6139499323161	85.5480625066467	85.4813765149857	85.4139426965139	85.3458117904120	85.2770345358609	85.2076616720413	85.1377439381341	85.0673320733199	84.9964768167797	84.9252289076941	84.8536390852441	84.7817580448157	84.7096284934900	84.6372759133962	84.5647235459648	84.4919946326263	84.4191124148114	84.3461001339506	84.2729810314745	84.1997783488137	84.1265153273987	84.0532152086602	83.9799012340288	83.9065966449349	83.8333246828091	83.7601085890822	83.6869716051845	83.6139369725468	83.5410346384077	83.4683311187708	83.3959053337275	83.3238362134189	83.2522026879858	83.1810836875690	83.1105581423095	83.0407049823481	82.9716031378257	82.9033315388831	82.8359691156614	82.7695947983012	82.7042875169436	82.6401262017294	82.5771897827995	82.5155571902947	82.4553045903829	82.3964637701095	82.3390307387157	82.2830004868022	82.2283680049696	82.1751282838187	82.1232763139500	82.0728070859642	82.0237155904621	81.9759968180442	81.9296457593113	81.8846574048639	81.8410267453028	81.7987487712286	81.7578184732420	81.7182308419437	81.6799809419201	81.6430714103353	81.6075189022342	81.5733415910352	81.5405576501567	81.5091852530171	81.4792425730347	81.4507477836279	81.4237190582152	81.3981745702150	81.3741324930455	81.3516110001253	81.3306282648727	81.3112024607061	81.2933517610438	81.2770943393044	81.2624483689061	81.2494270549409	81.2380206910871	81.2282129710027	81.2199875877353	81.2133282343327	81.2082186038425	81.2046423893125	81.2025832837901	81.2020249803232	81.2029511719594	81.2053455517464	81.2091918127318	81.2144736479633	81.2211747504886	81.2292788133554	81.2387695296113	81.2496303700177	81.2618419300221	81.2753828016687	81.2902315362660	81.3063666851222	81.3237667995457	81.3424104308448	81.3622761303277	81.3833424493028	81.4055879390785	81.4289911509629	81.4535306362644	81.4791849462913	81.5059326323519	81.5337522457546	81.5626223378076	81.5925213843446	81.6234229715139	81.6552928911261	81.6880962386969	81.7217981097420	81.7563635997771	81.7917578043181	81.8279458188806	81.8648927389805	81.9025636601333	81.9409236778550	81.9799378876613	82.0195713850678	82.0597892655904	82.1005566247447	82.1418385580466	82.1836001610131	82.2258092857760	82.2684445504549	82.3114872010805	82.3549184836831	82.3987196442931	82.4428719289410	82.4873565836572	82.5321548544722	82.5772479874164	82.6226172285202	82.6682438238140	82.7141090193284	82.7601940610937	82.8064801951403	82.8529486674988	82.8995807241994	82.9463566738246	82.9932469090548	83.0402158625736	83.0872278856533	83.1342473295661	83.1812385455845	83.2281658849807	83.2749936990270	83.3216863389958	83.3682081561594	83.4145235017901	83.4605967271601	83.5063921835419	83.5518742222078	83.5970071944299	83.6417554514808	83.6860834994845	83.7299627192817	83.7733739508586	83.8162987224979	83.8587185624824	83.9006149990946	83.9419695606174	83.9827637753333	84.0229791715251	84.0625972774755	84.1015996214671	84.1399677317827	84.1776831367050	84.2147273645165	84.2510819435001	84.2867284019384	84.3216482680956	84.3558229829549	84.3892336985313	84.4218615071180	84.4536875010086	84.4846927724963	84.5148584138745	84.5441655174366	84.5725951754760	84.6001284802861	84.6267465241601	84.6524303993916	84.6771611982738	84.7009200131001	84.7236879361640	84.7454460597587	84.7661754761777	84.7858579942200	84.8044816922316	84.8220378908752	84.8385179380353	84.8539131815962	84.8682149694423	84.8814146494580	84.8935035695278	84.9044730775360	84.9143145213671	84.9230192489054	84.9305786080353	84.9369839466413	84.9422266126077	84.9462979538190	84.9491893181595	84.9508922000413	84.9514028058121	84.9507229473945	84.9488547636578	84.9458003934714	84.9415619757047	84.9361416492271	84.9295415529078	84.9217638256163	84.9128106062220	84.9026840335942	84.8913862466023	84.8789193841157	84.8652855850037	84.8504869881357	84.8345257323811	84.8174039560411	84.7991233232559	84.7796841637181	84.7590865740776	84.7373306509839	84.7144164910867	84.6903441910357	84.6651138474804	84.6387255570706	84.6111794164559	84.5824755222859	84.5526139712104	84.5215948598788	84.4894182849410	84.4560843430465	84.4215931308451	84.3859447449862	84.3491412781378	84.3111993992910	84.2721422425426	84.2319929710545	84.1907747479889	84.1485107365079	84.1052240997736	84.0609380009482	84.0156756031938	83.9694600696724	83.9223145635462	83.8742622479774	83.8253262861279	83.7755298411600	83.7248960762358	83.6734481545173	83.6212086495648	83.5681858361838	83.5143732680497	83.4597638247378	83.4043503858237	83.3481258308830	83.2910830394912	83.2332148912238	83.1745142656563	83.1149740423642	83.0545871009231	82.9933463209085	82.9312445818959	82.8682747634608	82.8044297451788	82.7397024066253	82.6740856683352	82.6075848770240	82.5402352165138	82.4720762569085	82.4031475683121	82.3334887208284	82.2631392845614	82.1921388296151	82.1205269260933	82.0483431441000	81.9756270537392	81.9024182251147	81.8287562283305	81.7546806334906	81.6802310106987	81.6054469300590	81.5303679616753	81.4550307133975	81.3794536205548	81.3036482310642	81.2276260797609	81.1513987014803	81.0749776310576	80.9983744033280	80.9216005531270	80.8446676152896	80.7675871246513	80.6903706160472	80.6130296243126	80.5355756842829	80.4580203307932	80.3803750986789	80.3026515227752	80.2248605960010	80.1470030961903	80.0690706983956	79.9910547591411	79.9129466349509	79.8347376823491	79.7564192578599	79.6779827180072	79.5994194193153	79.5207207183083	79.4418779715102	79.3628825354452	79.2837257666374	79.2043990216109	79.1248936568898	79.0452010289983	78.9653126269066	78.8852398228520	78.8050348764548	78.7247550783306	78.6444577190950	78.5642000893635	78.4840394797518	78.4040331808754	78.3242384833498	78.2447126777908	78.1655130548138	78.0866969050344	78.0083215190682	77.9304441875308	77.8531222010377	77.7764128502046	77.7003734256469	77.6250476411489	77.5504089234982	77.4764080226095	77.4029956772535	77.3301226262007	77.2577396082219	77.1857973620877	77.1142466265687	77.0430381404355	76.9721226424588	76.9014508714093	76.8309735660575	76.7606414651742	76.6904053075298	76.6202158318952	76.5500237770409	76.4797839017536	76.4095112366381	76.3392672713546	76.2691146942198	76.1991161935509	76.1293344576648	76.0598321748785	75.9906720335090	75.9219167218732	75.8536289282882	75.7858713410708	75.7187066485382	75.6521975390073	75.5864067007950	75.5213968222183	75.4572305915943	75.3939703073196	75.3316439715703	75.2702189800172	75.2096565354554	75.1499178406800	75.0909640984862	75.0327565116690	74.9752562830234	74.9184246153447	74.8622227114279	74.8066117740681	74.7515530060605	74.6970076102000	74.6429367892818	74.5893017461011	74.5360636834529	74.4831838041322	74.4306341556452	74.3784342926939	74.3266167752349	74.2752141635549	74.2242590179410	74.1737838986801	74.1238213660591	74.0744039803650	74.0255643018847	73.9773348909050	73.9297483077130	73.8828371125955	73.8366338658396	73.7911711277320	73.7464814585597	73.7025974186098	73.6595493715981	73.6173410164582	73.5759582969886	73.5353868352480	73.4956122532949	73.4566201731882	73.4183962169863	73.3809260067481	73.3441951645321	73.3081893123970	73.2728940724016	73.2382950666043	73.2043779170640	73.1711282458393	73.1385316749888	73.1065738265713	73.0752404119470	73.0445222546065	73.0144179635545	72.9849268016115	72.9560480315979	72.9277809163343	72.9001247186410	72.8730787013386	72.8466421272475	72.8208142591882	72.7955943599812	72.7709816924469	72.7469755194059	72.7235751036784	72.7007797080852	72.6785885954465	72.6570010286526	72.6360210498812	72.6156704359557	72.5959750763938	72.5769608607134	72.5586536784323	72.5410794190683	72.5242639721392	72.5082332271627	72.4930130736568	72.4786294011392	72.4651080991277	72.4524750571402	72.4407561646944	72.4299773113081	72.4201643864992	72.4113432797855	72.4035347536071	72.3967084770204	72.3908047922058	72.3857636927294	72.3815251721568	72.3780292240541	72.3752158419869	72.3730250195212	72.3713967502228	72.3702710276576	72.3695878453913	72.3692871969899	72.3693090760192	72.3695934760450	72.3700803906332	72.3707098133496	72.3714229728845	72.3722104864981	72.3731279348148	72.3742353220920	72.3755926525869	72.3772599305569	72.3792971602592	72.3817643459513	72.3847214918903	72.3882286023336	72.3923456815386	72.3971327337625	72.4026497632627	72.4089567742964	72.4161137711210	72.4241807579939	72.4332177302606	72.4432630699808	72.4542871253461	72.4662468919371	72.4790993653345	72.4928015411189	72.5073104148708	72.5225829821709	72.5385762385998	72.5552471797382	72.5725528011665	72.5904500984656	72.6088960672158	72.6278477029980	72.6472620013927	72.6670959579804	72.6873065683419	72.7078547350230	72.7287336652290	72.7499525053896	72.7715205142934	72.7934469507287	72.8157410734840	72.8384121413478	72.8614694131084	72.8849221475544	72.9087796034741	72.9330510396560	72.9577457148885	72.9828728879602	73.0084418176593	73.0344617627745	73.0609419820941	73.0878912925815	73.1153055367305	73.1431657961149	73.1714523513285	73.2001454829652	73.2292254716187	73.2586725978830	73.2884671423519	73.3185893856193	73.3490196082790	73.3797380909249	73.4107251141508	73.4419609585506	73.4734259047181	73.5051002332473	73.5369642247318	73.5689981755734	73.6011915865978	73.6335586143870	73.6661175004888	73.6988864864511	73.7318838138218	73.7651277241487	73.7986364589797	73.8324282598628	73.8665213683458	73.9009340259766	73.9356844743030	73.9707909548730	74.0062717092344	74.0421449789352	74.0784290055231	74.1151420305461	74.1522988636232	74.1898905784169	74.2278982153560	74.2663027791578	74.3050852745393	74.3442267062176	74.3837080789098	74.4235103973330	74.4636146662044	74.5040018902409	74.5446530741598	74.5855492226781	74.6266713405129	74.6680004323813	74.7095175030004	74.7512035570873	74.7930399618071	74.8350161936897	74.8771297145201	74.9193783237080	74.9617598206627	75.0042720047937	75.0469126755106	75.0896796322227	75.1325706743396	75.1755836012708	75.2187162124257	75.2619663072139	75.3053316850447	75.3488101453278	75.3923994874724	75.4360975108883	75.4799020100437	75.5238096015103	75.5678142066969	75.6119093717044	75.6560886426340	75.7003455655867	75.7446736866634	75.7890665519654	75.8335177075936	75.8780206996490	75.9225690742328	75.9671563774460	76.0117761553896	76.0564219541646	76.1010873198722	76.1457657986134	76.1904509364892	76.2351382236587	76.2798344634523	76.3245505406988	76.3692973457211	76.4140857688424	76.4589267003854	76.5038310306734	76.5488096500292	76.5938734487758	76.6390333172362	76.6843001457334	76.7296848245904	76.7751982441302	76.8208512946758	76.8666548665501	76.9126198500763	76.9587557213825	77.0050471712264	77.0514577680656	77.0979504046563	77.1444879737545	77.1910333681164	77.2375494804980	77.2839992036554	77.3303454303447	77.3765510533221	77.4225789653435	77.4683920591651	77.5139532275430	77.5592253632332	77.6041713589918	77.6487541075750	77.6929366595321	77.7367018687116	77.7800714380074	77.8230715884651	77.8657285411302	77.9080685170481	77.9501177372644	77.9919024228246	78.0334487947742	78.0747830741587	78.1159314820236	78.1569202394144	78.1977755673767	78.2385236869559	78.2791908191977	78.3198031851473	78.3603870058505	78.4009627617711	78.4415227129847	78.4820504621707	78.5225296097204	78.5629437560256	78.6032765014777	78.6435114464685	78.6836321913893	78.7236223366319	78.7634654825878	78.8031452296486	78.8426451782059	78.8819489286512	78.9210400813760	78.9599022367722	78.9985189952310	79.0368737945090	79.0749477916219	79.1127204628580	79.1501712453898	79.1872795763898	79.2240248930305	79.2603866324843	79.2963442319238	79.3318771285215	79.3669647594497	79.4015865618811	79.4357219729880	79.4693504299430	79.5024513699186	79.5350042300872	79.5669884476213	79.5983836917328	79.6291873438528	79.6594266736011	79.6891318288081	79.7183329573047	79.7470602069213	79.7753437254887	79.8032136608374	79.8307001607981	79.8578333732014	79.8846434458779	79.9111605266584	79.9374147633733	79.9634363038534	79.9892552959292	80.0149018874314	80.0404062261907	80.0657872120670	80.0910169330352	80.1160553008431	80.1408622272114	80.1653976238607	80.1896214025118	80.2134934748854	80.2369737527021	80.2600221476826	80.2825985715476	80.3046629360177	80.3261751528138	80.3470951336564	80.3673827902662	80.3869980343640	80.4059007776704	80.4240547790827	80.4414676833565	80.4581750560021	80.4742129105466	80.4896172605170	80.5044241194402	80.5186695008432	80.5323894182532	80.5456198851970	80.5583969152016	80.5707565217942	80.5827347185016	80.5943675188510	80.6056909363692	80.6167409845834	80.6275536770204	80.6381644851829	80.6485812805709	80.6587717673960	80.6687004829171	80.6783319643936	80.6876307490844	80.6965613742488	80.7050883771460	80.7131762950350	80.7207896651751	80.7278930248253	80.7344509112449	80.7404278616930	80.7457884134287	80.7504971037112	80.7545184697996	80.7578170496460	80.7603678791959	80.7621830119658	80.7632826565322	80.7636870214717	80.7634163153610	80.7624907467765	80.7609305242948	80.7587558564925	80.7559869519461	80.7526440192322	80.7487472669275	80.7443169036084	80.7393731378515	80.7339361782334	80.7280262333306	80.7216635117198	80.7148667627410	80.7076410215345	80.6999838079792	80.6918925649395	80.6833647352800	80.6743977618652	80.6649890875595	80.6551361552274	80.6448364077335	80.6340872879422	80.6228862387180	80.6112307029255	80.5991181234291	80.5865459430932	80.5735116047825	80.5600125513614	80.5460460389226	80.5316025649766	80.5166641273419	80.5012121830979	80.4852281893242	80.4686936031002	80.4515898815055	80.4338984816194	80.4156008605214	80.3966784752911	80.3771127830078	80.3568852407512	80.3359773056006	80.3143704346355	80.2920460849354	80.2689857135798	80.2451707844157	80.2205924035611	80.1952705436525	80.1692305561915	80.1424977926799	80.1150976046194	80.0870553435117	80.0583963608584	80.0291460081614	79.9993296369224	79.9689725986429	79.9381002448248	79.9067379269698	79.8749109965795	79.8426448051557	79.8099647042001	79.7768960452143	79.7434596401627	79.7096408048988	79.6754081510578	79.6407301928694	79.6055754445630	79.5699124203683	79.5337096345149	79.4969356012323	79.4595588347501	79.4215478492979	79.3828711591053	79.3434972784019	79.3033947214172	79.2625320023808	79.2208776355223	79.1784001350713	79.1350692449158	79.0908877987161	79.0458946340594	79.0001304023236	78.9536357548870	78.9064513431276	78.8586178184235	78.8101758321527	78.7611660356934	78.7116290804237	78.6616056177217	78.6111362989654	78.5602617755330	78.5090226988025	78.4574597201520	78.4056134909596	78.3535246165080	78.3012140795740	78.2486528155226	78.1958038963578	78.1426303940832	78.0890953807029	78.0351619282206	77.9807931086402	77.9259519939656	77.8706016562006	77.8147051673491	77.7582255994149	77.7011260244019	77.6433695143139	77.5849191411548	77.5257379769285	77.4657890936388	77.4050429635479	77.3435185592784	77.2812543847212	77.2182889973811	77.1546609547629	77.0904088143714	77.0255711337114	76.9601864702878	76.8942933816053	76.8279304251688	76.7611361584830	76.6939491390529	76.6264079243832	76.5585510719786	76.4904171393441	76.4220446839845	76.3534705890956	76.2846970866156	76.2156937739237	76.1464289778841	76.0768710253612	76.0069882432194	75.9367489583229	75.8661214975362	75.7950741877236	75.7235753557495	75.6515933284781	75.5790964327738	75.5060529955011	75.4324313435241	75.3581998037074	75.2833267029151	75.2077804836400	75.1315517516193	75.0546794546648	74.9772089133843	74.8991854483858	74.8206543802772	74.7416610296663	74.6622507171610	74.5824687633692	74.5023604888989	74.4219712143578	74.3413462603539	74.2605309474950	74.1795705963891	74.0985105276440	74.0173960618676	73.9362725196679	73.8551852216526	73.7741794884297	73.6933006406071	73.6125939987926	73.5321048835941	73.4518786156196	73.3719605154769	73.2923959037739	73.2132301011184	73.1345084281185	73.0562762053819	72.9785787535165	72.9014613931303	72.8249694448310];
        PRBPOWER = PRBPOWER*xlpln(PRBWL-399)/100;
        % Pellicle beamsplitter power correction (1-2500 nm)
        BStrans = [0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.0316300000000000	0.805310000000000	0.118730000000000	0.0373600000000000	0.139960000000000	0.152130000000000	0.472740000000000	1.13583000000000	0.131530000000000	0.725900000000000	0.619680000000000	0.655050000000000	0.955960000000000	0.235500000000000	1.11524000000000	1.20048000000000	0.350210000000000	0.377190000000000	0.255320000000000	1.00664000000000	1.45784000000000	1.62808000000000	1.33865000000000	2.16181000000000	2.11515000000000	1.10618000000000	0.528770000000000	2.72680000000000	0.803080000000000	1.77604000000000	0.233680000000000	0.632840000000000	0.365060000000000	0.566720000000000	0.698830000000000	1.29796000000000	0.953840000000000	0.574910000000000	0.0752000000000000	0.127720000000000	0.208090000000000	0.413980000000000	1.59303000000000	0.524520000000000	0.425960000000000	2.73864000000000	0.346680000000000	0.697960000000000	0.439850000000000	1.39002000000000	0.233760000000000	1.10112000000000	1.39920000000000	1.06977000000000	3.07187000000000	0.386960000000000	0.713600000000000	3.05694000000000	2.65171000000000	0.752580000000000	0.797670000000000	1.45916000000000	1.10274000000000	1.79768000000000	2.90613000000000	1.20756000000000	1.92496000000000	0.300060000000000	1.06546000000000	0.209610000000000	3.57832000000000	2.59083000000000	0.328390000000000	3.13160000000000	0.765140000000000	1.89232000000000	1.63980000000000	1.41732000000000	1.56093000000000	2.99158000000000	1.55941000000000	2.28477000000000	2.71539000000000	1.71160000000000	2.10356000000000	2.16121000000000	2.97121000000000	2.58634000000000	1.87521000000000	2.53224000000000	1.46936000000000	2.28168000000000	2.84737000000000	1.16128000000000	2.44232000000000	2.55651000000000	2.82912000000000	3.50579000000000	3.38702000000000	3.29702000000000	2.79215000000000	2.95893000000000	3.80206000000000	3.53211000000000	3.75801000000000	3.54704000000000	3.78942000000000	3.73001000000000	3.14379000000000	3.37963000000000	4.21475000000000	4.09041000000000	4.32873000000000	4.63199000000000	4.69174000000000	4.96309000000000	5.23248000000000	5.21629000000000	5.14598000000000	5.27854000000000	5.68675000000000	5.35653000000000	5.68165000000000	5.74316000000000	5.64435000000000	5.89168000000000	6.09394000000000	6.16011000000000	6.10264000000000	6.26718000000000	6.37395000000000	6.55372000000000	6.71860000000000	6.95152000000000	7.53432000000000	7.67449000000000	7.80717000000000	8.21823000000000	8.43010000000000	8.54389000000000	8.62407000000000	8.87726000000000	8.91881000000000	9.30336000000000	9.20459000000000	9.31761000000000	9.46774000000000	9.62059000000000	9.45815000000000	9.70941000000000	9.89090000000000	9.99207000000000	10.2195300000000	10.4203400000000	10.7602300000000	10.9549500000000	11.4074300000000	11.8140800000000	12.2212600000000	12.6464600000000	12.9797500000000	13.4344800000000	13.8005700000000	14.2208500000000	14.5315000000000	14.8044600000000	15.0291400000000	15.3423700000000	15.5084700000000	15.6750100000000	15.8700600000000	15.9151600000000	16.0290700000000	16.2408300000000	16.3319500000000	16.4517200000000	16.6469300000000	16.8835500000000	17.1340000000000	17.4305800000000	17.8193500000000	18.1811400000000	18.5587000000000	19.0380700000000	19.4999800000000	20.1102300000000	20.5887800000000	21.2014200000000	21.8644000000000	22.4015700000000	23.0323000000000	23.5552100000000	24.1055700000000	24.5849800000000	25.0831600000000	25.4955900000000	25.8744000000000	26.1880200000000	26.4478100000000	26.6497100000000	26.9136600000000	27.0945500000000	27.2155100000000	27.3580600000000	27.4732500000000	27.6355400000000	27.7650400000000	27.9251600000000	28.1383000000000	28.3860200000000	28.6518700000000	28.9637400000000	29.4033100000000	29.7927700000000	30.3007800000000	30.8350300000000	31.4006500000000	32.0432900000000	32.7095500000000	33.4496000000000	34.2009300000000	34.9854900000000	35.7744300000000	36.5456800000000	37.3178000000000	38.0703200000000	38.8128600000000	39.4867600000000	40.0723300000000	40.6353200000000	41.0789500000000	41.4648700000000	41.7648400000000	41.9811500000000	42.1324700000000	42.1651300000000	42.1651800000000	42.1271200000000	42.0427500000000	41.9053600000000	41.7568700000000	41.6183900000000	41.4953900000000	41.4019300000000	41.3190300000000	41.3131500000000	41.3651700000000	41.4497800000000	41.6075100000000	41.8199200000000	42.1149200000000	42.4859800000000	42.9300200000000	43.4430000000000	44.0337600000000	44.6801600000000	45.4184300000000	46.1636300000000	46.9929100000000	47.8399300000000	48.7142300000000	49.5853700000000	50.4656400000000	51.2979700000000	52.1373600000000	52.8789400000000	53.5729200000000	54.1486900000000	54.6338000000000	55.0393900000000	55.3120600000000	55.4608400000000	55.4839600000000	55.4123900000000	55.2262300000000	54.9812600000000	54.5690800000000	54.1543700000000	53.6800700000000	53.1645500000000	52.6270100000000	52.0769500000000	51.5319900000000	51.0502200000000	50.5771500000000	50.1439600000000	49.7586700000000	49.4733200000000	49.2306400000000	49.0983000000000	48.9929000000000	49.0318500000000	49.1072600000000	49.2562700000000	49.5091700000000	49.8485200000000	50.2562200000000	50.7005200000000	51.2846400000000	51.8789900000000	52.5552700000000	53.3106700000000	54.0455500000000	54.8715600000000	55.7319900000000	56.5612900000000	57.4187700000000	58.2384500000000	59.0384600000000	59.8261600000000	60.5151200000000	61.1648100000000	61.7216200000000	62.2043400000000	62.5943100000000	62.8406800000000	62.9597500000000	63.0261800000000	62.9597700000000	62.7498500000000	62.4835100000000	62.1461900000000	61.6738700000000	61.1734500000000	60.6013900000000	59.9792800000000	59.3066100000000	58.6289300000000	57.9993100000000	57.3145300000000	56.6539600000000	56.0303300000000	55.4365800000000	54.8827100000000	54.3884900000000	53.9197800000000	53.5303700000000	53.2085000000000	52.9554300000000	52.7343300000000	52.6082700000000	52.5424100000000	52.5459400000000	52.6049400000000	52.7388900000000	52.9468700000000	53.1821300000000	53.5200300000000	53.9015700000000	54.3329800000000	54.8370700000000	55.3925100000000	55.9826400000000	56.6071000000000	57.2777100000000	57.9950800000000	58.6963100000000	59.4378100000000	60.1956400000000	60.9437800000000	61.6912000000000	62.4112700000000	63.1039900000000	63.7844900000000	64.3995700000000	64.9614700000000	65.4753000000000	65.9319500000000	66.3158300000000	66.6307800000000	66.8604700000000	67.0074300000000	67.0758000000000	67.0470200000000	66.9580600000000	66.7890400000000	66.5259800000000	66.2303900000000	65.8534900000000	65.4143800000000	64.9372600000000	64.4042000000000	63.8612700000000	63.2852200000000	62.6992700000000	62.0762900000000	61.4773400000000	60.8726800000000	60.2879300000000	59.7029600000000	59.1499500000000	58.6098300000000	58.1087200000000	57.6543200000000	57.2060900000000	56.8105600000000	56.4764900000000	56.1591000000000	55.8913500000000	55.6894500000000	55.5120400000000	55.3809900000000	55.3109800000000	55.2796400000000	55.2834700000000	55.3470200000000	55.4639500000000	55.6127400000000	55.8016300000000	56.0559700000000	56.3368700000000	56.6625000000000	57.0349800000000	57.4230500000000	57.8584200000000	58.3456300000000	58.8340800000000	59.3637800000000	59.9076100000000	60.4893300000000	61.0836900000000	61.6904800000000	62.3101700000000	62.9205600000000	63.5673400000000	64.1903000000000	64.7973900000000	65.4056100000000	66.0130500000000	66.5755000000000	67.0990400000000	67.6325000000000	68.1211300000000	68.5556600000000	68.9600100000000	69.3189900000000	69.6312600000000	69.8908800000000	70.0802500000000	70.2262400000000	70.3318300000000	70.3915000000000	70.3584600000000	70.3034100000000	70.1878900000000	70.0262800000000	69.7983900000000	69.5601700000000	69.2657900000000	68.9260000000000	68.5426900000000	68.1591200000000	67.7350300000000	67.3043000000000	66.8586100000000	66.3652300000000	65.8908500000000	65.4010800000000	64.9080900000000	64.4214200000000	63.9374200000000	63.4722700000000	63.0077300000000	62.5556100000000	62.1424500000000	61.7269100000000	61.3357500000000	60.9633600000000	60.6074200000000	60.2894500000000	59.8285700000000	59.7337000000000	59.5278300000000	59.2201700000000	59.0416000000000	58.8560800000000	58.7805300000000	58.6816600000000	58.5369800000000	58.5326700000000	58.5056400000000	58.5143200000000	58.6497800000000	58.6830100000000	58.7960100000000	58.9550700000000	59.0177200000000	59.2726200000000	59.4560000000000	59.7647800000000	60.0421500000000	60.2989400000000	60.6761300000000	61.0014800000000	61.2969200000000	61.7732500000000	62.1680500000000	62.5946600000000	63.0636900000000	63.4786600000000	63.9452300000000	64.4365700000000	64.8737900000000	65.3899200000000	65.9317000000000	66.4414500000000	66.9271200000000	67.4019700000000	67.9191200000000	68.4187900000000	68.9263300000000	69.4078400000000	69.9438700000000	70.3174100000000	70.7315100000000	71.1699400000000	71.6641000000000	71.9967100000000	72.3735300000000	72.7750400000000	73.0256100000000	73.3261000000000	73.5485200000000	73.8297700000000	73.9296600000000	74.0600400000000	74.1816600000000	74.2702200000000	74.3850200000000	74.3566200000000	74.2523300000000	74.1702000000000	74.2203400000000	74.1391400000000	73.9375200000000	73.7000900000000	73.5623000000000	73.3478600000000	72.9575100000000	72.7903400000000	72.3787500000000	72.2166200000000	71.8501500000000	71.4184000000000	71.0751100000000	70.7836800000000	70.3535200000000	70.0246300000000	69.6448100000000	69.2006900000000	68.9215300000000	68.5279700000000	68.1364200000000	67.7739700000000	67.3762400000000	66.9377700000000	66.6366000000000	66.3682400000000	66.0111900000000	65.6616600000000	65.5021700000000	65.0728700000000	64.9130900000000	64.5432400000000	64.2681300000000	64.0390900000000	63.8369000000000	63.7470900000000	63.5123400000000	63.2602800000000	63.2048100000000	63.0392500000000	62.8867600000000	62.8816500000000	62.8218300000000	62.7282400000000	62.6017900000000	62.6328700000000	62.6375100000000	62.6121700000000	62.7235600000000	62.7708400000000	62.7979900000000	62.9004800000000	63.0106000000000	63.1344800000000	63.2813100000000	63.4021100000000	63.6024900000000	63.7928400000000	63.9152700000000	64.1572300000000	64.4038500000000	64.6617800000000	64.9003600000000	65.2075800000000	65.4333300000000	65.7587800000000	66.0430800000000	66.3889500000000	66.7427000000000	67.0366700000000	67.4128000000000	67.7781100000000	68.1070800000000	68.4868700000000	68.8803100000000	69.2281200000000	69.6544500000000	70.0348700000000	70.4306000000000	70.8292800000000	71.2449600000000	71.6390200000000	72.0232400000000	72.4065100000000	72.8034400000000	73.1466700000000	73.5678700000000	73.9422200000000	74.2763900000000	74.6321400000000	75.0075600000000	75.3090600000000	75.6562200000000	75.9754600000000	76.2669500000000	76.5563800000000	76.8446400000000	77.0889900000000	77.3298600000000	77.5725600000000	77.7654300000000	77.9431900000000	78.1296700000000	78.2832000000000	78.4065100000000	78.5386300000000	78.6217100000000	78.6890300000000	78.7334800000000	78.7741700000000	78.7771800000000	78.7893100000000	78.7594100000000	78.7123900000000	78.6457700000000	78.5798000000000	78.4678600000000	78.3715700000000	78.2360600000000	78.0988300000000	77.9262800000000	77.7728800000000	77.5741600000000	77.3807700000000	77.1990100000000	76.9673100000000	76.7391700000000	76.5087500000000	76.2761200000000	76.0247200000000	75.7738500000000	75.5114100000000	75.2511900000000	74.9710500000000	74.7055600000000	74.4354400000000	74.1621900000000	73.8710600000000	73.5875900000000	73.3127600000000	73.0427100000000	72.7782400000000	72.5011400000000	72.2345400000000	71.9712400000000	71.7098200000000	71.4615000000000	71.2138200000000	70.9717100000000	70.7231900000000	70.4922800000000	70.2738500000000	70.0590100000000	69.8505800000000	69.6510800000000	69.4377600000000	69.2624400000000	69.0874500000000	68.9205900000000	68.7660600000000	68.6066400000000	68.4645100000000	68.3226000000000	68.2093100000000	68.0919800000000	67.9815700000000	67.8971000000000	67.8025300000000	67.7217700000000	67.6643800000000	67.6194500000000	67.5595600000000	67.5356600000000	67.5232800000000	67.5064800000000	67.5238800000000	67.5041300000000	67.5346300000000	67.5417600000000	67.5781600000000	67.6186500000000	67.6909700000000	67.7618000000000	67.8030100000000	67.8851700000000	68.0221200000000	68.0886100000000	68.2042900000000	68.3405600000000	68.4725300000000	68.5984500000000	68.7543000000000	68.9077400000000	69.0763500000000	69.2388500000000	69.4172700000000	69.5960500000000	69.7980900000000	70.0050500000000	70.2008100000000	70.4276100000000	70.6453200000000	70.8656200000000	71.1166500000000	71.3590500000000	71.5946600000000	71.8511400000000	72.1110600000000	72.3717300000000	72.6242400000000	72.9142000000000	73.1785600000000	73.4498000000000	73.7216400000000	74.0337500000000	74.2898800000000	74.5650300000000	74.8575900000000	75.1347000000000	75.4323200000000	75.7371700000000	76.0007300000000	76.2877200000000	76.5830100000000	76.8694400000000	77.1433600000000	77.4377800000000	77.7156100000000	77.9950600000000	78.2694100000000	78.5323700000000	78.7942400000000	79.0717400000000	79.3485800000000	79.5943800000000	79.8315400000000	80.0715300000000	80.3066600000000	80.5454700000000	80.7557400000000	80.9632600000000	81.1827900000000	81.3851900000000	81.5795600000000	81.7823700000000	81.9537500000000	82.1189400000000	82.2680800000000	82.4291300000000	82.5803400000000	82.7130100000000	82.8304500000000	82.9582500000000	83.0716900000000	83.1497600000000	83.2406000000000	83.3262200000000	83.3846700000000	83.4410400000000	83.5040100000000	83.5105500000000	83.5408100000000	83.5628400000000	83.5841800000000	83.5555500000000	83.5565900000000	83.5357300000000	83.4892700000000	83.5542400000000	83.4033300000000	83.5002800000000	83.3018200000000	83.2995800000000	83.2153200000000	83.1046000000000	83.0074200000000	82.8490900000000	82.6786400000000	82.6315000000000	82.6258700000000	82.4539800000000	82.2716400000000	82.0721500000000	82.0479700000000	81.7958300000000	81.6185800000000	81.4838800000000	81.3010600000000	81.2186100000000	81.0824800000000	80.8095300000000	80.5952000000000	80.4976700000000	80.1553000000000	80.1486100000000	79.9055800000000	79.8081600000000	79.5559700000000	79.3588300000000	79.1483300000000	78.9702800000000	78.8638900000000	78.4725600000000	78.5061800000000	78.2178400000000	78.1082500000000	77.8408600000000	77.6479300000000	77.4766500000000	77.1428900000000	77.0814600000000	76.9172900000000	76.6884000000000	76.5940700000000	76.3782100000000	76.2570100000000	76.0793800000000	75.9256400000000	75.7993000000000	75.6477400000000	75.4841400000000	75.2828800000000	75.1442000000000	75.0491000000000	74.8671300000000	74.6724600000000	74.6333500000000	74.5600400000000	74.3413900000000	74.2860900000000	74.1018400000000	74.0635200000000	73.9547500000000	73.8129000000000	73.7734700000000	73.6545700000000	73.6224600000000	73.5343200000000	73.4324000000000	73.3878900000000	73.2873600000000	73.3745000000000	73.2121600000000	73.1592800000000	73.1691800000000	73.0272400000000	72.9974600000000	72.9938400000000	73.0028100000000	72.9673000000000	72.9941100000000	72.9362300000000	72.9623800000000	73.0004000000000	72.9341000000000	72.9262000000000	73.0302300000000	73.0378300000000	73.0334800000000	73.0969500000000	73.0818800000000	73.1065300000000	73.1545300000000	73.1981300000000	73.2272900000000	73.3711700000000	73.3985100000000	73.4929900000000	73.5709800000000	73.6112500000000	73.7055500000000	73.7298000000000	73.8705400000000	73.9231300000000	74.0895600000000	74.2056300000000	74.2766900000000	74.3844700000000	74.4418200000000	74.6003100000000	74.6500500000000	74.8510300000000	74.9920000000000	75.1298500000000	75.2890000000000	75.3285100000000	75.5353200000000	75.6969800000000	75.7949600000000	75.9471200000000	76.1240700000000	76.2656600000000	76.3890900000000	76.6264100000000	76.7327900000000	76.9115300000000	77.0582900000000	77.2494400000000	77.4404700000000	77.5303200000000	77.8336900000000	77.9712900000000	78.1167000000000	78.3677300000000	78.5355800000000	78.6897100000000	78.8991400000000	79.0726600000000	79.2223700000000	79.4412800000000	79.6719400000000	79.8336700000000	80.0321400000000	80.1652600000000	80.4921300000000	80.5827600000000	80.7399600000000	80.9377200000000	81.1405900000000	81.3213300000000	81.5190700000000	81.7476700000000	82.0311300000000	82.1034000000000	82.2981400000000	82.4951100000000	82.6200300000000	82.9078800000000	83.0256700000000	83.2452700000000	83.4666300000000	83.5731700000000	83.7844000000000	83.9620300000000	84.1168100000000	84.2092600000000	84.5180800000000	84.6348700000000	84.8063500000000	84.9630700000000	85.0365400000000	85.3106500000000	85.3977100000000	85.5977900000000	85.7358600000000	85.8990200000000	85.9763300000000	86.1910000000000	86.2890900000000	86.3975400000000	86.4997500000000	86.6196000000000	86.8308800000000	86.8985100000000	86.9917800000000	87.0927900000000	87.1754000000000	87.2985700000000	87.3623400000000	87.4524000000000	87.5649300000000	87.6429800000000	87.8037900000000	87.8549400000000	87.8844900000000	87.9634400000000	88.0319800000000	88.0971600000000	88.1348300000000	88.1743500000000	88.2187400000000	88.2057400000000	88.2862200000000	88.3163000000000	88.3135800000000	88.3202300000000	88.3642000000000	88.3692000000000	88.4435100000000	88.3769100000000	88.4002500000000	88.4329500000000	88.3677800000000	88.3298600000000	88.4189400000000	88.3374600000000	88.3747900000000	88.2562600000000	88.3109900000000	88.1957700000000	88.1528500000000	88.1445900000000	88.0638000000000	88.0721800000000	88.0040300000000	87.9102800000000	87.8462100000000	87.8521400000000	87.7338800000000	87.7176800000000	87.6390600000000	87.5968000000000	87.4665000000000	87.3660900000000	87.3076400000000	87.2000100000000	87.0726900000000	87.0424400000000	86.9079000000000	86.8288300000000	86.8273800000000	86.6488900000000	86.5819700000000	86.4155700000000	86.3354300000000	86.2528800000000	86.1142600000000	86.0634800000000	85.9559600000000	85.8808700000000	85.7473500000000	85.5605900000000	85.4696300000000	85.4477200000000	85.2838500000000	85.1907000000000	85.0072800000000	84.8837600000000	84.8313100000000	84.6395500000000	84.6243400000000	84.4698700000000	84.2919800000000	84.2529400000000	84.1033200000000	84.0015800000000	83.9397400000000	83.7178200000000	83.6133300000000	83.5490500000000	83.3769800000000	83.2776200000000	83.1775300000000	83.0398300000000	82.9432900000000	82.8521000000000	82.7398000000000	82.6231800000000	82.5384800000000	82.3682700000000	82.2934300000000	82.1647300000000	82.0373200000000	81.9860600000000	81.8918200000000	81.7906200000000	81.6123600000000	81.5495600000000	81.4868600000000	81.3570300000000	81.2687500000000	81.1790300000000	81.0521600000000	80.9575200000000	80.8844100000000	80.8210900000000	80.6850000000000	80.6581600000000	80.5217700000000	80.4473200000000	80.3426000000000	80.3104200000000	80.2073400000000	80.1340700000000	80.1133400000000	80.0440500000000	79.9165500000000	79.7660100000000	79.7604800000000	79.7119600000000	79.6538000000000	79.5944400000000	79.4669500000000	79.4491800000000	79.4364900000000	79.3743500000000	79.3410100000000	79.2492500000000	79.1997100000000	79.1666200000000	79.1273800000000	79.0773700000000	79.0990600000000	78.9676100000000	78.9805700000000	79.0142500000000	78.8646500000000	78.8065600000000	78.8565100000000	78.8284600000000	78.8162800000000	78.7894300000000	78.8027200000000	78.7858700000000	78.6898700000000	78.7074700000000	78.7510500000000	78.7465000000000	78.7117000000000	78.7449300000000	78.6896400000000	78.6422800000000	78.7104100000000	78.7417400000000	78.6837800000000	78.7320900000000	78.7319100000000	78.7369800000000	78.7382200000000	78.7529700000000	78.7877000000000	78.8283800000000	78.8282200000000	78.8733200000000	78.8056400000000	78.8615000000000	78.8998300000000	78.9399300000000	78.9793300000000	79.0102500000000	79.0410400000000	79.0927200000000	79.1106700000000	79.1693600000000	79.2075000000000	79.2319300000000	79.2809900000000	79.3313800000000	79.4099800000000	79.4571200000000	79.4630100000000	79.5410100000000	79.5505800000000	79.5868200000000	79.7484500000000	79.7101000000000	79.8143200000000	79.8641500000000	79.9499400000000	80.0799700000000	80.0946100000000	80.1916000000000	80.2112700000000	80.3403700000000	80.4000700000000	80.4340000000000	80.5790700000000	80.6528500000000	80.7281300000000	80.8062400000000	80.8920000000000	80.9503600000000	81.0716900000000	81.1375000000000	81.2718400000000	81.3440300000000	81.4298500000000	81.5668900000000	81.6356500000000	81.7161200000000	81.7978100000000	81.8678200000000	81.9803100000000	82.1020000000000	82.1498100000000	82.3352000000000	82.3624000000000	82.4678800000000	82.6052400000000	82.7639500000000	82.8462400000000	82.8961600000000	83.0465800000000	83.1106300000000	83.2149700000000	83.3509600000000	83.4790900000000	83.6155500000000	83.6660100000000	83.8251000000000	83.8639100000000	83.9690200000000	84.0945400000000	84.2244600000000	84.3503000000000	84.4332100000000	84.5851200000000	84.6561700000000	84.7781400000000	84.9106900000000	84.9616600000000	85.1115300000000	85.2407200000000	85.3592300000000	85.4884900000000	85.5538200000000	85.6863400000000	85.8112400000000	85.9041100000000	85.9614700000000	86.1448800000000	86.2362200000000	86.3786600000000	86.4640100000000	86.5883800000000	86.7111100000000	86.8090600000000	86.9232300000000	87.0430800000000	87.1584500000000	87.3021600000000	87.3341400000000	87.4790800000000	87.6139100000000	87.6873900000000	87.8019200000000	87.9116400000000	87.9931000000000	88.1805000000000	88.2086600000000	88.3592100000000	88.4537200000000	88.5551900000000	88.6374500000000	88.7384800000000	88.8887300000000	88.9198700000000	89.0320500000000	89.1652900000000	89.2450000000000	89.3420300000000	89.4459500000000	89.5522400000000	89.6335900000000	89.6678300000000	89.8102700000000	89.9378800000000	89.9978300000000	90.0449400000000	90.1640200000000	90.2261300000000	90.3599500000000	90.4428600000000	90.5277900000000	90.5788300000000	90.6555600000000	90.7239500000000	90.8233000000000	90.8902500000000	90.9708200000000	91.0272200000000	91.0690700000000	91.1505200000000	91.2541700000000	91.2820800000000	91.4006200000000	91.4325700000000	91.4661000000000	91.6065000000000	91.5628300000000	91.6459100000000	91.6764800000000	91.8398100000000	91.8160900000000	91.8715400000000	91.8916000000000	91.9855100000000	92.0317900000000	92.0969200000000	92.1490600000000	92.1723300000000	92.2159700000000	92.2258800000000	92.2776800000000	92.4126500000000	92.3947400000000	92.3994100000000	92.4416900000000	92.4478100000000	92.4485900000000	92.5086400000000	92.5393000000000	92.6452100000000	92.6305500000000	92.6440600000000	92.6722400000000	92.6852000000000	92.6779000000000	92.7335300000000	92.7399400000000	92.7850300000000	92.7909600000000	92.7798800000000	92.8366700000000	92.8042400000000	92.7630900000000	92.8312900000000	92.7959000000000	92.8086600000000	92.8600300000000	92.8354600000000	92.8010800000000	92.8086900000000	92.8380800000000	92.7418900000000	92.8065000000000	92.8128700000000	92.7808100000000	92.7490500000000	92.7426700000000	92.7969200000000	92.7249500000000	92.7248500000000	92.6453900000000	92.6568300000000	92.5830700000000	92.6059000000000	92.5395600000000	92.5220600000000	92.4631000000000	92.4702600000000	92.4349800000000	92.4438200000000	92.3943600000000	92.3484900000000	92.2915200000000	92.2884400000000	92.2473400000000	92.1845100000000	92.1523500000000	92.1602700000000	92.0732700000000	92.0524600000000	92.0774000000000	91.9785100000000	91.9275100000000	91.8265800000000	91.8502000000000	91.8356500000000	91.7655100000000	91.6610600000000	91.6281100000000	91.5847300000000	91.5678500000000	91.5006600000000	91.4717700000000	91.4189300000000	91.3459200000000	91.3908200000000	91.2439000000000	91.2373800000000	91.1801900000000	91.0923700000000	91.0382800000000	90.9656400000000	90.9865100000000	90.9401900000000	90.8497000000000	90.7519300000000	90.6727500000000	90.6133700000000	90.5489300000000	90.5177300000000	90.4305800000000	90.4442600000000	90.3731100000000	90.2653100000000	90.2489400000000	90.1626300000000	90.0758900000000	90.0361100000000	89.9906800000000	89.8944500000000	89.8161000000000	89.7949300000000	89.6584600000000	89.5611600000000	89.6264000000000	89.5030700000000	89.4628400000000	89.4274600000000	89.2933700000000	89.2764700000000	89.2399000000000	89.1827600000000	89.0651600000000	89.0432600000000	88.9947700000000	88.9361600000000	88.8100100000000	88.7307600000000	88.6340900000000	88.6849700000000	88.6256000000000	88.5510000000000	88.4840300000000	88.3881600000000	88.3092000000000	88.2745600000000	88.1702200000000	88.1558700000000	88.0343700000000	88.0641300000000	88.0481900000000	87.8593900000000	87.9126800000000	87.7720900000000	87.7092300000000	87.7098000000000	87.6020300000000	87.6088900000000	87.4952500000000	87.3989300000000	87.4295500000000	87.2857100000000	87.2378500000000	87.2077400000000	87.1439300000000	87.1014600000000	87.0140500000000	86.9129900000000	86.9010400000000	86.8213800000000	86.7848100000000	86.7643400000000	86.7078700000000	86.6643200000000	86.5981900000000	86.5574200000000	86.4458500000000	86.4461000000000	86.4386400000000	86.4173300000000	86.2970400000000	86.2509400000000	86.2540100000000	86.1629300000000	86.1536600000000	86.0909900000000	86.0149400000000	85.9964800000000	85.9879800000000	85.9289500000000	85.8533600000000	85.8578500000000	85.7779800000000	85.8163500000000	85.6662800000000	85.6995300000000	85.5979000000000	85.5559400000000	85.6237600000000	85.5344900000000	85.5247400000000	85.4224200000000	85.4319300000000	85.3659200000000	85.2743800000000	85.3265400000000	85.2612800000000	85.3357700000000	85.1733600000000	85.1732500000000	85.1623700000000	85.1378600000000	85.1068000000000	85.0542200000000	85.0293200000000	85.0241700000000	84.8877100000000	84.9658400000000	84.8456300000000	84.8459500000000	84.8918600000000	84.7950400000000	84.7044900000000	84.7643900000000	84.7659800000000	84.7517200000000	84.7671100000000	84.7027900000000	84.6819100000000	84.7061200000000	84.6364900000000	84.5799600000000	84.6452600000000	84.5078600000000	84.4802700000000	84.5759000000000	84.5826300000000	84.5567200000000	84.5105800000000	84.4767800000000	84.4832700000000	84.4648100000000	84.4046900000000	84.4243900000000	84.5063500000000	84.4402500000000	84.3939100000000	84.4044600000000	84.4336200000000	84.4281500000000	84.4215200000000	84.3251100000000	84.3452200000000	84.3389300000000	84.3632400000000	84.3696400000000	84.4097100000000	84.3427000000000	84.3948100000000	84.4034500000000	84.4298000000000	84.4219900000000	84.5113200000000	84.4158900000000	84.3997300000000	84.4271200000000	84.4231800000000	84.4090000000000	84.4506900000000	84.4713300000000	84.4648500000000	84.4802700000000	84.4987500000000	84.5946900000000	84.4539300000000	84.5636500000000	84.4220700000000	84.4418800000000	84.5268700000000	84.6556200000000	84.4433200000000	84.6548500000000	84.6032700000000	84.5975300000000	84.6653400000000	84.6683600000000	84.6569300000000	84.7500800000000	84.6352300000000	84.5877100000000	84.7150400000000	84.7959500000000	84.7978600000000	84.7340700000000	84.8033800000000	84.7872200000000	84.7977400000000	84.9133500000000	84.8568600000000	84.9516600000000	84.9912900000000	84.9656100000000	84.9159200000000	85.1057100000000	85.1260600000000	85.0165600000000	85.1180300000000	85.1591200000000	85.1521600000000	85.1648000000000	85.1797900000000	85.2842700000000	85.4006400000000	85.3263000000000	85.4085700000000	85.4267400000000	85.4947700000000	85.4685200000000	85.5332900000000	85.5180800000000	85.5260000000000	85.6823000000000	85.6903100000000	85.6617900000000	85.7580300000000	85.7425500000000	85.8746200000000	85.7839300000000	85.8115500000000	85.8679200000000	86.0572100000000	86.0268000000000	86.0132400000000	85.9280500000000	86.1130500000000	86.1444400000000	86.1912800000000	86.1568800000000	86.3279500000000	86.3057900000000	86.3560700000000	86.3125000000000	86.5318500000000	86.3930400000000	86.4527400000000	86.3444400000000	86.5925700000000	86.5192500000000	86.6460000000000	86.6202600000000	86.7534500000000	86.7494400000000	86.5693600000000	86.7695400000000	86.7100100000000	86.7022300000000	86.8690100000000	86.9965500000000	86.9766600000000	87.1341300000000	87.1253300000000	87.2444200000000	87.3041300000000	87.3340400000000	87.3282700000000	87.2743300000000	87.4933500000000	87.3738300000000	87.5199700000000	87.8115600000000	87.6434200000000	87.7006100000000	87.7552800000000	87.9572400000000	87.8617300000000	87.9007700000000	88.1211500000000	87.9979200000000	88.1737000000000	88.2172000000000	88.3103900000000	88.4443800000000	88.3516300000000	88.3205400000000	88.5012700000000	88.3982100000000	88.6968800000000	88.5847900000000	88.7182800000000	88.7309500000000	88.6854200000000	88.7701700000000	88.9035900000000	88.8616200000000	89.1121300000000	89.0622300000000	89.1406400000000	89.1976500000000	89.2231900000000	89.2623700000000	89.2456400000000	89.4251700000000	89.5194600000000	89.4858900000000	89.6339700000000	89.5267000000000	89.6050000000000	89.6476100000000	89.8502000000000	89.7357500000000	89.8417000000000	89.9231600000000	90.0881200000000	90.0997800000000	90.1007700000000	90.0533800000000	90.1259500000000	90.2581800000000	90.4073000000000	90.4377100000000	90.4693700000000	90.5277400000000	90.4906500000000	90.6233600000000	90.6948200000000	90.7440900000000	90.6161200000000	90.9085000000000	90.8768800000000	91.1223500000000	91.0969300000000	91.0067700000000	91.1381400000000	91.2088800000000	91.2390400000000	91.3985500000000	91.4490400000000	91.3494900000000	91.4082000000000	91.5070200000000	91.6213600000000	91.7514000000000	91.5548700000000	91.6571000000000	91.6887300000000	91.8265200000000	91.9421500000000	91.7378800000000	92.1016200000000	92.0283600000000	92.1060100000000	92.1937100000000	92.2692600000000	92.2788700000000	92.4186000000000	92.3643500000000	92.5810400000000	92.3132600000000	92.3938200000000	92.5299900000000	92.4739300000000	92.5361800000000	92.7276200000000	92.6985800000000	92.7388200000000	92.7939100000000	92.8436700000000	92.8269400000000	92.8799600000000	92.9848900000000	92.9362100000000	93.0098400000000	93.2184100000000	92.9133900000000	93.2075500000000	93.3084100000000	93.4873400000000	93.3121100000000	93.4637500000000	93.3802000000000	93.4309800000000	93.5155200000000	93.6279800000000	93.7298700000000	93.7179600000000	93.7894700000000	93.5967300000000	93.8424500000000	93.8914000000000	93.9106700000000	93.7086600000000	93.9993100000000	93.9385000000000	94.1419300000000	94.1729600000000	94.0983000000000	94.0577500000000	93.8949300000000	94.2774400000000	94.4227100000000	94.3669400000000	94.4522400000000	94.4461100000000	94.4826000000000	94.3750900000000	94.5050200000000	94.5100000000000	94.5652500000000	94.7032000000000	94.6241200000000	94.6442600000000	94.6104500000000	94.6719700000000	94.7332100000000	94.8208000000000	94.9155700000000	94.8767200000000	94.7355100000000	94.8263800000000	94.9460700000000	95.0013000000000	94.9107700000000	95.1277200000000	95.0777200000000	95.1236300000000	95.1771700000000	95.2015800000000	95.1438200000000	95.2718400000000	95.3016800000000	95.1440700000000	95.2857800000000	95.3604100000000	95.3027700000000	95.0430800000000	95.5567900000000	95.4342000000000	95.4810600000000	95.4028700000000	95.4121400000000	95.5573500000000	95.5931400000000	95.5727500000000	95.6897000000000	95.6767500000000	95.6631100000000	95.6524000000000	95.6648900000000	95.8687600000000	95.7126200000000	95.8553200000000	95.8323100000000	95.7560000000000	95.8552200000000	95.9190200000000	95.7872200000000	95.8508000000000	95.9666400000000	96.0822400000000	95.9370100000000	95.9169300000000	95.7954800000000	95.7185100000000	96.1343300000000	96.0591300000000	96.0212100000000	95.8725100000000	96.0403600000000	95.9495200000000	96.1336300000000	96.0696000000000	96.1536900000000	96.1771800000000	96.0965700000000	96.2021100000000	96.1030800000000	95.8495100000000	96.0488400000000	96.3746000000000	96.1940000000000	96.5109000000000	96.2363000000000	96.2945600000000	96.2582300000000	96.1972800000000	96.1012500000000	96.1977800000000	96.1771700000000	96.2688500000000	96.2134900000000	96.4011400000000	96.2567900000000	96.2939500000000	96.1386100000000	96.2826200000000	96.2392400000000	96.3840100000000	96.1135200000000	96.2759400000000	96.1191900000000	95.9944500000000	96.3328700000000	96.1936900000000	96.2365500000000	96.1925600000000	96.2558400000000	96.3860200000000	96.2261400000000	96.3669000000000	96.2316900000000	96.1114200000000	96.2170600000000	96.2696400000000	96.2813100000000	96.3484100000000	95.9716700000000	96.2487000000000	96.0151400000000	96.2384800000000	96.2946100000000	96.0454200000000	96.1324400000000	96.2971800000000	96.2283100000000	96.0812000000000	96.1798900000000	96.0042700000000	96.0521200000000	96.1201100000000	95.9834000000000	96.2049300000000	96.1082400000000	96.1897400000000	95.9657400000000	96.0489700000000	96.1428500000000	96.0934400000000	96.0244800000000	95.9911100000000	96.0435000000000	96.0494800000000	95.7285800000000	95.8462400000000	95.7960700000000	95.8452000000000	95.9276000000000	95.8741500000000	95.9089400000000	95.7770200000000	96.1609900000000	95.8094200000000	96.1091600000000	95.7969300000000	95.8244000000000	95.6538300000000	95.7354100000000	95.6572000000000	95.7396300000000	95.6526000000000	95.7140700000000	95.6530500000000	95.7663800000000	95.7833600000000	95.5131600000000	95.5768000000000	95.6681400000000	95.5085400000000	95.4855300000000	95.6511900000000	95.5241200000000	95.6151100000000	95.4209600000000	95.5487300000000	95.3369700000000	95.3205500000000	95.2950100000000	95.3161200000000	95.4305300000000	95.1585500000000	94.9449800000000	95.1589100000000	95.4731000000000	95.2760000000000	95.2037800000000	95.0891700000000	95.2051800000000	95.1488900000000	95.1710100000000	94.9248200000000	94.9211300000000	95.0754500000000	94.9651000000000	94.7547800000000	94.8121500000000	94.7204900000000	94.8088100000000	94.6752500000000	94.3471700000000	94.8583600000000	94.8774700000000	94.4959900000000	94.5998600000000	94.5610200000000	94.6785100000000	94.6621700000000	94.5546100000000	94.4685400000000	94.5185700000000	94.5717500000000	94.2494700000000	94.2950500000000	94.2633600000000	94.2354300000000	94.1289900000000	94.6396000000000	94.1129500000000	94.1710700000000	93.9914000000000	93.8625300000000	94.0079900000000	94.3154600000000	93.7877500000000	93.6828000000000	93.9718200000000	94.4337200000000	93.7845100000000	94.4710300000000	94.0344900000000	93.6987700000000	93.5494700000000	94.2173500000000	93.8378800000000	93.8737000000000	93.8140900000000	93.5826000000000	93.9200500000000	93.6200600000000	93.9101900000000	93.3416300000000	93.6411600000000	93.7811400000000	93.2447600000000	93.5460800000000	93.2610200000000	93.5301700000000	93.4854000000000	93.6601600000000	93.2393300000000	93.2928700000000	93.0603400000000	93.2741500000000	93.1185700000000	93.1016200000000	92.9314300000000	93.0502300000000	92.9839700000000	92.9866600000000	92.9529100000000	92.9893900000000	92.8573500000000	92.8240900000000	92.8031000000000	92.9205800000000	92.6946600000000	92.9572800000000	92.7476700000000	92.7254400000000	92.6662300000000	92.7348700000000	92.6384700000000	92.5802400000000	92.5814400000000	92.6058500000000	92.5272400000000	92.6192000000000	92.5296400000000	92.4160400000000	92.4591200000000	92.4816800000000	92.3549700000000	92.3728900000000	92.3003800000000	92.2378400000000	92.1591800000000	92.2422600000000	92.2384100000000	92.2910000000000	92.2805400000000	92.0833900000000	92.2010000000000	92.1152900000000	92.1209500000000	92.0213200000000	91.9537300000000	91.9525100000000	91.9060700000000	91.9363700000000	91.8821700000000	91.8672800000000	91.7849500000000	91.8140000000000	91.7664200000000	91.7199500000000	91.6622500000000	91.5594000000000	91.6986100000000	91.6465700000000	91.6750300000000	91.5669800000000	91.6132800000000	91.5378400000000	91.5335300000000	91.4733600000000	91.4841600000000	91.5236700000000	91.5552200000000	91.5295600000000	91.4868100000000	91.3899400000000	91.6097400000000	91.4843200000000	91.4753900000000	91.4251200000000	91.3814200000000	91.3580200000000	91.3979000000000	91.3979800000000	91.4092400000000	91.3131200000000	91.2973200000000	91.2028000000000	91.2844600000000	91.2893200000000	91.1803000000000	91.1506800000000	91.1338400000000	91.1913800000000	91.1287300000000	91.0468900000000	91.0677100000000	91.1436200000000	90.9859900000000	90.9763900000000	90.9901000000000	90.9853700000000	90.9407300000000	90.9539400000000	90.9653100000000	90.9334300000000	90.8008200000000	90.9486200000000	90.7737200000000	90.8490500000000	90.7938900000000	90.8091900000000	90.6567500000000	90.6500100000000	90.7652000000000	90.7568700000000	90.7563200000000	90.7470500000000	90.6800200000000	90.7123000000000	90.6943700000000	90.5535900000000	90.5409200000000	90.5712700000000	90.6577700000000	90.6235800000000	90.3629500000000	90.6147100000000	90.4302300000000	90.5963900000000	90.3898800000000	90.3898200000000	90.4619100000000	90.4115800000000	90.3395700000000	90.3372200000000	90.2689400000000	90.1382700000000	90.2753100000000	90.3037600000000	90.3116000000000	90.2014300000000	90.1875400000000	90.2053600000000	90.1494000000000	90.1391700000000	90.0507600000000	90.1369300000000	90.1088300000000	90.0046700000000	90.1794300000000	90.0036200000000	89.9799400000000	90.1236400000000	89.9851200000000	89.9089000000000	89.9618400000000	89.9938300000000	89.8693600000000	89.9954800000000	89.9886600000000	89.9159700000000	90.0371600000000	89.8083200000000	89.8111900000000	89.7599000000000	89.8781500000000	89.8394500000000	89.8595900000000	89.7913300000000	89.8114100000000	89.7060100000000	89.9371400000000	89.7348300000000	89.7851700000000	89.7314500000000	89.6612300000000	89.6584200000000	89.4354100000000	89.7822100000000	89.7089200000000	89.7078600000000	89.6366500000000	89.7856100000000	89.5686400000000	89.7355700000000	89.5941500000000	89.4623000000000	89.4520000000000	89.5643100000000	89.5645500000000	89.3775300000000	89.4671400000000	89.4349700000000	89.4814800000000	89.3946900000000	89.3148200000000	89.5879400000000	89.3359900000000	89.3050800000000	89.4592700000000	89.2665300000000	89.3324400000000	89.3233000000000	89.3769800000000	89.4765300000000	89.1013400000000	89.4418300000000	89.3112000000000	89.2604700000000	89.3734100000000	89.5398600000000	89.4178200000000];
        BSref = [0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	9.55232000000000	9.39092000000000	5.41764000000000	5.29360000000000	7.28716000000000	9.96871000000000	10.6633900000000	12.2526200000000	10.2393700000000	8.78400000000000	5.79283000000000	6.20194000000000	8.50825000000000	12.3557400000000	12.8999300000000	15.2487400000000	14.6960200000000	10.6187800000000	7.93268000000000	7.92405000000000	9.10017000000000	8.67855000000000	11.4495400000000	12.8809900000000	16.2445700000000	13.8376600000000	13.5221600000000	10.4644500000000	9.47150000000000	9.31693000000000	10.6454100000000	11.2038000000000	13.8457000000000	15.6186100000000	16.4196600000000	18.9209100000000	16.1149900000000	14.4317800000000	10.7327200000000	10.5637600000000	8.56792000000000	11.4125000000000	13.0389200000000	15.2756500000000	17.6898400000000	20.1903000000000	20.0823100000000	20.7111600000000	18.4032400000000	16.3474400000000	12.6801700000000	12.9871700000000	11.4736800000000	10.9725000000000	13.8404700000000	15.9820700000000	18.1234400000000	21.0921900000000	21.3619800000000	21.2881300000000	21.5114700000000	17.6580300000000	16.9971000000000	15.7498100000000	13.6271700000000	11.6341400000000	11.1611500000000	11.2115300000000	13.8128400000000	16.1580200000000	19.4843800000000	20.9506900000000	19.5592300000000	21.8721800000000	20.4329700000000	19.5984700000000	18.3043500000000	16.8381900000000	14.0390100000000	13.7971700000000	11.4424400000000	10.0218100000000	11.3226500000000	13.7370200000000	14.8619400000000	16.2958700000000	18.7446600000000	21.1421600000000	21.4866700000000	21.8122800000000	20.5056200000000	20.2397900000000	19.5750500000000	16.9367000000000	15.7546800000000	13.2434000000000	11.5911600000000	11.6092100000000	10.3883400000000	11.1398800000000	12.3656000000000	13.3104900000000	15.3514500000000	17.3633500000000	18.7150000000000	20.4682600000000	20.6407700000000	21.1186100000000	21.1514000000000	20.4501800000000	19.1785800000000	17.5218600000000	15.8047800000000	14.2897300000000	12.6851900000000	11.2324200000000	10.0957800000000	10.0183200000000	9.98029000000000	10.4378800000000	11.7819600000000	13.0332500000000	14.3559100000000	16.1435500000000	17.3246800000000	18.7949000000000	19.8235400000000	20.0944000000000	20.2452100000000	19.8155100000000	18.9959800000000	17.8220100000000	16.4699700000000	14.7818200000000	13.5414000000000	11.8876700000000	10.6317400000000	9.51846000000000	8.62878000000000	8.32590000000000	8.32393000000000	8.75975000000000	9.58412000000000	10.8579500000000	11.8338400000000	13.2019700000000	14.4926500000000	15.6624700000000	16.4346800000000	17.2631300000000	17.7594000000000	18.0045600000000	17.8754400000000	17.4644100000000	16.7780700000000	15.8799900000000	14.6838700000000	13.4561900000000	12.2768000000000	10.9765400000000	9.52680000000000	8.70149000000000	7.70356000000000	7.14439000000000	6.73476000000000	6.60779000000000	6.82936000000000	7.31764000000000	8.04111000000000	8.82619000000000	9.80264000000000	10.7810800000000	11.9006500000000	12.9123100000000	13.8414600000000	14.5822800000000	15.1986700000000	15.6503600000000	16.0023500000000	16.0149500000000	15.8254300000000	15.5248700000000	14.9953400000000	14.3217600000000	13.4756200000000	12.5998600000000	11.6739900000000	10.6441500000000	9.73905000000000	8.76672000000000	7.95658000000000	7.22418000000000	6.65561000000000	6.26617000000000	6.11787000000000	6.08521000000000	6.32614000000000	6.69447000000000	7.27556000000000	8.00795000000000	8.91167000000000	9.88981000000000	10.8466600000000	11.8477100000000	12.9026300000000	13.8535900000000	14.7738300000000	15.5338300000000	16.2532700000000	16.8336400000000	17.2957300000000	17.5682400000000	17.7482000000000	17.7676700000000	17.6290500000000	17.3803100000000	17.0061600000000	16.4239500000000	15.8698300000000	15.2236300000000	14.5213300000000	13.7427300000000	12.9855700000000	12.2445500000000	11.6061100000000	11.0123100000000	10.5716500000000	10.2353400000000	10.0390600000000	10.0300900000000	10.2090800000000	10.5604400000000	11.1180500000000	11.8258200000000	12.6674900000000	13.6630000000000	14.7727400000000	15.9590200000000	17.2029500000000	18.4800400000000	19.7508600000000	20.9621700000000	22.1154800000000	23.2388100000000	24.2406700000000	25.1508500000000	25.9592300000000	26.5848500000000	27.0839100000000	27.4706700000000	27.7347400000000	27.8251100000000	27.8251300000000	27.6589400000000	27.4005700000000	27.0175200000000	26.5708200000000	26.0130900000000	25.3781400000000	24.7416700000000	24.0504500000000	23.3811200000000	22.7556500000000	22.1267700000000	21.5838800000000	21.1656700000000	20.8411300000000	20.6204800000000	20.5727300000000	20.6950600000000	20.9504200000000	21.4177400000000	22.0247200000000	22.7365300000000	23.6010100000000	24.6162800000000	25.6606800000000	26.8134300000000	27.9968700000000	29.2384400000000	30.4651800000000	31.6503500000000	32.8534500000000	33.9572700000000	35.0090100000000	35.9778100000000	36.8359200000000	37.6282600000000	38.2788700000000	38.8234100000000	39.2626200000000	39.6340600000000	39.8363800000000	39.9430300000000	39.9648500000000	39.8716200000000	39.6968300000000	39.4114200000000	39.0239500000000	38.5741700000000	38.0639700000000	37.4791800000000	36.8550500000000	36.1299600000000	35.4679900000000	34.7468700000000	34.0406100000000	33.3568300000000	32.6740500000000	32.0553300000000	31.5056100000000	31.0445300000000	30.6834800000000	30.4045600000000	30.2455900000000	30.1804200000000	30.2693300000000	30.4580200000000	30.7541600000000	31.1519500000000	31.6962100000000	32.3011200000000	32.9993600000000	33.7733300000000	34.5700600000000	35.4085500000000	36.2548100000000	37.1549600000000	38.0536500000000	38.9119500000000	39.7154800000000	40.5520900000000	41.2813300000000	41.9715400000000	42.6464000000000	43.1936400000000	43.7040400000000	44.1438900000000	44.5202300000000	44.8016700000000	45.0186500000000	45.1769600000000	45.2505400000000	45.2511700000000	45.1962800000000	45.0554100000000	44.8614900000000	44.5788000000000	44.2379100000000	43.8490800000000	43.3786300000000	42.8910900000000	42.3252400000000	41.7210800000000	41.0910700000000	40.4265500000000	39.7288200000000	39.0248200000000	38.3091200000000	37.5914500000000	36.8749200000000	36.1868400000000	35.5236400000000	34.8755000000000	34.3049400000000	33.7594300000000	33.2798600000000	32.8605000000000	32.5252200000000	32.2542400000000	32.0695800000000	31.9589000000000	31.9300900000000	31.9860800000000	32.1237900000000	32.3403800000000	32.6167600000000	32.9558600000000	33.3734800000000	33.8514300000000	34.3472500000000	34.8984100000000	35.4789400000000	36.0661300000000	36.6997800000000	37.3215200000000	37.9406700000000	38.5593600000000	39.1648200000000	39.7666200000000	40.3193800000000	40.8687900000000	41.3826300000000	41.8591300000000	42.2998800000000	42.6934200000000	43.0553400000000	43.3739400000000	43.6564800000000	43.8708200000000	44.0529900000000	44.1907500000000	44.2672400000000	44.2943400000000	44.2848700000000	44.2298800000000	44.1314200000000	43.9700900000000	43.7815600000000	43.5515400000000	43.2611000000000	42.9453400000000	42.5851500000000	42.1694300000000	41.7474500000000	41.2749300000000	40.7733400000000	40.2562100000000	39.6914600000000	39.1252100000000	38.5299000000000	37.9336400000000	37.3146700000000	36.6882500000000	36.0690300000000	35.4473500000000	34.8354400000000	34.2431500000000	33.6533200000000	33.0930400000000	32.5511600000000	32.0516000000000	31.5776800000000	31.1284800000000	30.7540100000000	30.4064100000000	30.1118500000000	29.8538500000000	29.6501700000000	29.5093200000000	29.4241300000000	29.3877600000000	29.4057100000000	29.4731900000000	29.5913500000000	29.7751000000000	29.9793500000000	30.2480000000000	30.5469300000000	30.8879300000000	31.2388200000000	31.6470600000000	32.0607000000000	32.5291300000000	32.9783800000000	33.4438600000000	33.9306900000000	34.4217900000000	34.8889100000000	35.4003400000000	35.8704400000000	36.3467700000000	36.8189100000000	37.2440500000000	37.6890100000000	38.0913400000000	38.4922700000000	38.8552600000000	39.1756000000000	39.5209800000000	39.9283200000000	40.1108100000000	40.4188200000000	40.6146200000000	40.9421500000000	41.0372700000000	41.0302900000000	41.2058700000000	41.3054600000000	41.3425900000000	41.3102600000000	41.3050000000000	41.2874300000000	41.1545000000000	41.1213700000000	41.0314600000000	40.6840000000000	40.6180300000000	40.3244700000000	40.2014000000000	39.9038500000000	39.4740400000000	39.2257000000000	38.8805700000000	38.5805600000000	38.1330200000000	37.7576200000000	37.3540000000000	36.8941100000000	36.3990700000000	35.9629100000000	35.4805900000000	34.9491000000000	34.4828400000000	33.9335900000000	33.4038200000000	32.8292000000000	32.3993500000000	31.9296700000000	31.4717600000000	30.9604600000000	30.4476300000000	30.0180000000000	29.5143400000000	29.1507800000000	28.7034200000000	28.2432100000000	27.9468700000000	27.5631100000000	27.1507200000000	26.8635400000000	26.6037900000000	26.3863400000000	26.1376900000000	25.9699900000000	25.9078300000000	25.6579100000000	25.6322700000000	25.5928300000000	25.6620200000000	25.7291100000000	25.5831800000000	25.7589300000000	25.8845100000000	25.9795800000000	26.3240100000000	26.5418800000000	26.6929700000000	26.9890000000000	27.2519300000000	27.5155900000000	27.8326400000000	28.1309000000000	28.5608200000000	28.8149600000000	29.3355800000000	29.5269200000000	29.9008400000000	30.4010900000000	30.7754400000000	31.1239200000000	31.4964600000000	31.8873000000000	32.1915900000000	32.5903400000000	32.9774100000000	33.2319700000000	33.6858000000000	33.9580000000000	34.2490700000000	34.5704600000000	34.8616100000000	35.1344200000000	35.4652500000000	35.6724700000000	35.9184100000000	36.1072500000000	36.3687200000000	36.5029500000000	36.6064400000000	36.8231800000000	36.9994100000000	36.9681600000000	37.1771800000000	37.2371500000000	37.2888800000000	37.3365100000000	37.3757700000000	37.3984300000000	37.3569300000000	37.3415000000000	37.2370600000000	37.1568300000000	37.0657100000000	37.0535400000000	36.8677600000000	36.7898100000000	36.5490200000000	36.4432400000000	36.2697100000000	36.0069000000000	35.8376300000000	35.5701400000000	35.3388000000000	35.1086700000000	34.8454000000000	34.5318900000000	34.2554100000000	33.9283900000000	33.5990100000000	33.3136600000000	32.9292000000000	32.6109100000000	32.2539600000000	31.9178300000000	31.4786100000000	31.1173400000000	30.7347100000000	30.3581100000000	29.9689500000000	29.5729200000000	29.2139100000000	28.7685500000000	28.3794800000000	27.9981300000000	27.5859600000000	27.2331800000000	26.8342300000000	26.4817400000000	26.0875600000000	25.7383600000000	25.3695200000000	25.0353300000000	24.7025100000000	24.3486600000000	24.0299800000000	23.7569200000000	23.4841900000000	23.2005400000000	22.9705400000000	22.7287600000000	22.4986200000000	22.3105500000000	22.1222500000000	21.9638000000000	21.8124900000000	21.6758100000000	21.5701300000000	21.4653400000000	21.4165700000000	21.3584600000000	21.3231700000000	21.3107500000000	21.3204000000000	21.3501800000000	21.4014800000000	21.4615600000000	21.5358400000000	21.6391900000000	21.7604800000000	21.8967100000000	22.0332000000000	22.1967600000000	22.3672900000000	22.5467900000000	22.7509000000000	22.9659300000000	23.1891600000000	23.4102600000000	23.6613300000000	23.9045600000000	24.1459800000000	24.4126500000000	24.6962400000000	24.9613500000000	25.2378700000000	25.5183600000000	25.7830600000000	26.0760400000000	26.3643800000000	26.6363200000000	26.9226800000000	27.1952000000000	27.4800600000000	27.7473900000000	28.0438300000000	28.3137400000000	28.5669100000000	28.8374000000000	29.0877100000000	29.3320800000000	29.5865200000000	29.8265800000000	30.0444800000000	30.2706400000000	30.4813500000000	30.6768000000000	30.8774600000000	31.0828000000000	31.2610500000000	31.4304600000000	31.5981900000000	31.7514500000000	31.8911500000000	32.0315500000000	32.1608000000000	32.2691500000000	32.3817600000000	32.4762800000000	32.5481100000000	32.6383200000000	32.6997300000000	32.7492600000000	32.7960600000000	32.8288300000000	32.8561000000000	32.8661900000000	32.8758400000000	32.8579800000000	32.8538400000000	32.8185900000000	32.7768500000000	32.7328000000000	32.6805200000000	32.6108200000000	32.5214900000000	32.4442300000000	32.3692200000000	32.2513400000000	32.1417600000000	32.0195000000000	31.8898300000000	31.7468400000000	31.6014700000000	31.4550100000000	31.2890500000000	31.1091300000000	30.9257100000000	30.7374900000000	30.5426200000000	30.3432600000000	30.1349500000000	29.9123000000000	29.6834900000000	29.4543700000000	29.2137900000000	28.9752200000000	28.7279200000000	28.4697600000000	28.2005000000000	27.9401700000000	27.6761300000000	27.4036100000000	27.1414500000000	26.8589000000000	26.5611600000000	26.2785600000000	25.9892000000000	25.7033700000000	25.4098900000000	25.1310500000000	24.8250000000000	24.5232800000000	24.2373500000000	23.9371200000000	23.6529300000000	23.3713700000000	23.0725400000000	22.7792300000000	22.4948900000000	22.2066200000000	21.9258000000000	21.6564600000000	21.3744800000000	21.0977400000000	20.8302600000000	20.5697200000000	20.3166900000000	20.0687800000000	19.8360100000000	19.5870700000000	19.3507900000000	19.1373600000000	18.9133400000000	18.7051200000000	18.5093600000000	18.3095300000000	18.1180800000000	17.9431000000000	17.7755900000000	17.6196900000000	17.4722600000000	17.3315400000000	17.1910800000000	17.0710000000000	16.9619500000000	16.8561900000000	16.7636400000000	16.6865200000000	16.6076600000000	16.5441700000000	16.4931100000000	16.4467500000000	16.4134800000000	16.3877100000000	16.3738800000000	16.3690600000000	16.3758100000000	16.3899800000000	16.4129800000000	16.1276700000000	16.2335400000000	16.5202900000000	16.3081300000000	16.6473900000000	16.6947900000000	16.6187100000000	16.6993700000000	16.9713900000000	17.0187800000000	17.2032500000000	17.3131900000000	17.3835100000000	17.6352000000000	17.8412200000000	17.9730100000000	17.9938300000000	18.1633000000000	18.3226900000000	18.5853500000000	18.7101900000000	19.0007700000000	19.0973200000000	19.3006800000000	19.5444500000000	19.6823200000000	19.8950700000000	20.0130300000000	20.2985000000000	20.4087700000000	20.6849800000000	20.8255000000000	20.9544600000000	21.2183500000000	21.3489300000000	21.5905400000000	21.8307400000000	22.0615200000000	22.2824000000000	22.2926300000000	22.4715900000000	22.6768300000000	22.9412800000000	23.1169600000000	23.2668800000000	23.3939800000000	23.7353200000000	23.8768300000000	24.0122900000000	24.1905300000000	24.3348800000000	24.5762300000000	24.6380000000000	24.7727500000000	25.0508400000000	25.0783800000000	25.2169200000000	25.3809900000000	25.5651800000000	25.7004100000000	25.7422600000000	25.9023400000000	25.9448800000000	26.1714200000000	26.2024300000000	26.3436300000000	26.4341600000000	26.4896900000000	26.6576100000000	26.7480600000000	26.7786500000000	26.8743200000000	26.9689200000000	26.9575500000000	27.0335200000000	27.1177900000000	27.1879300000000	27.1972000000000	27.2536500000000	27.2327700000000	27.3878500000000	27.2630200000000	27.4157100000000	27.3901500000000	27.4000300000000	27.4087600000000	27.3536800000000	27.3678100000000	27.3257700000000	27.4247600000000	27.3490000000000	27.3722900000000	27.2750800000000	27.2615800000000	27.2689000000000	27.1589100000000	27.0992600000000	27.0456800000000	27.0247400000000	26.9522000000000	26.8421700000000	26.8526900000000	26.7992200000000	26.6695500000000	26.5847300000000	26.4691600000000	26.4429000000000	26.3309900000000	26.1763400000000	26.1079300000000	25.9574900000000	25.7868100000000	25.7062500000000	25.5913800000000	25.5277300000000	25.3902400000000	25.2698800000000	25.1136400000000	24.9765700000000	24.8164900000000	24.6198900000000	24.5832300000000	24.4365100000000	24.2085200000000	24.0006400000000	23.9752300000000	23.7541100000000	23.5531500000000	23.3936400000000	23.2353000000000	23.0514500000000	22.8625400000000	22.7657000000000	22.5628000000000	22.3823000000000	22.2581300000000	22.0010400000000	21.8460800000000	21.6113400000000	21.4863400000000	21.2585200000000	21.1191300000000	20.9222200000000	20.6902100000000	20.4917900000000	20.3320900000000	20.1428200000000	19.9050400000000	19.6770400000000	19.5147200000000	19.3816000000000	19.0961500000000	19.0050700000000	18.7623000000000	18.6548900000000	18.3102700000000	18.1441000000000	18.0394200000000	17.9001800000000	17.6337100000000	17.4429400000000	17.3242300000000	17.1468600000000	16.9012300000000	16.6424500000000	16.4884100000000	16.4137800000000	16.2147400000000	16.0349200000000	15.8807700000000	15.7068200000000	15.5242000000000	15.3667500000000	15.2303000000000	15.0554500000000	14.8681600000000	14.7777100000000	14.6581300000000	14.4671900000000	14.3777500000000	14.1799900000000	14.0692600000000	13.9659200000000	13.7669100000000	13.7292500000000	13.6032300000000	13.4788600000000	13.3747000000000	13.2202100000000	13.1869700000000	13.0670300000000	13.0089200000000	12.9334400000000	12.8105700000000	12.7153100000000	12.6769900000000	12.5684800000000	12.5261500000000	12.4621300000000	12.4783700000000	12.3454000000000	12.3473900000000	12.2480500000000	12.2515700000000	12.2293400000000	12.1542900000000	12.1289100000000	12.1472700000000	12.1070800000000	12.1095300000000	12.1532900000000	12.2334200000000	12.1846700000000	12.1306200000000	12.1690800000000	12.1574200000000	12.1730600000000	12.2123800000000	12.1561400000000	12.2085700000000	12.2987300000000	12.3484200000000	12.3851600000000	12.4051100000000	12.4379700000000	12.4493800000000	12.5747900000000	12.6060900000000	12.6849300000000	12.7406500000000	12.8441900000000	12.8822100000000	12.9761400000000	13.0000200000000	13.1148200000000	13.1755800000000	13.2876000000000	13.4097300000000	13.4490400000000	13.5508000000000	13.6831100000000	13.7842700000000	13.8909800000000	13.9217900000000	14.0794000000000	14.1594900000000	14.2871200000000	14.4160300000000	14.5030000000000	14.6153900000000	14.7531600000000	14.9079300000000	14.9210500000000	15.0134200000000	15.1750700000000	15.3002600000000	15.4198700000000	15.5411600000000	15.5923300000000	15.7459800000000	15.8823000000000	15.9301600000000	16.0941800000000	16.2352800000000	16.3615700000000	16.4724400000000	16.5474000000000	16.7271800000000	16.8668800000000	16.9637500000000	17.0941200000000	17.1895800000000	17.2963800000000	17.4141200000000	17.5682100000000	17.6603900000000	17.7351000000000	17.9234600000000	18.0416700000000	18.1387500000000	18.2640600000000	18.3112700000000	18.3944900000000	18.5530700000000	18.7225000000000	18.7536700000000	18.8888600000000	18.9467200000000	19.0494200000000	19.1933500000000	19.2316900000000	19.4187300000000	19.4515100000000	19.5228800000000	19.6268800000000	19.7332200000000	19.8208400000000	19.9157200000000	20.0204800000000	20.0980400000000	20.2011700000000	20.2057200000000	20.3503600000000	20.3921300000000	20.3973800000000	20.5344800000000	20.5825300000000	20.6426400000000	20.7357900000000	20.8283500000000	20.8620900000000	20.9745300000000	21.0408400000000	21.1053900000000	21.1453600000000	21.1700400000000	21.1923900000000	21.2552400000000	21.3261400000000	21.3301000000000	21.3609100000000	21.4506400000000	21.4560600000000	21.5018100000000	21.5638700000000	21.5434500000000	21.6220500000000	21.6164400000000	21.7016500000000	21.6926300000000	21.7498300000000	21.7495700000000	21.7587900000000	21.7202500000000	21.8144500000000	21.7343000000000	21.7693700000000	21.8504100000000	21.8172500000000	21.7822600000000	21.7852500000000	21.8374400000000	21.8441000000000	21.8149700000000	21.8220300000000	21.7646200000000	21.7634000000000	21.7262600000000	21.7234800000000	21.6906500000000	21.6704100000000	21.6736100000000	21.6226900000000	21.6579100000000	21.6807900000000	21.5365000000000	21.5389500000000	21.4903200000000	21.4575500000000	21.4311600000000	21.3793300000000	21.3559100000000	21.3213700000000	21.2617600000000	21.2113800000000	21.1467700000000	21.0993300000000	21.0492100000000	21.0337000000000	20.9750300000000	20.8732700000000	20.8412700000000	20.7606300000000	20.7638800000000	20.6464700000000	20.5887600000000	20.5091700000000	20.4814800000000	20.4231000000000	20.3187200000000	20.2480600000000	20.1793000000000	20.0717200000000	20.0058800000000	19.9403600000000	19.8622500000000	19.7875100000000	19.6708900000000	19.6016100000000	19.5646300000000	19.4596700000000	19.3748900000000	19.2746200000000	19.1580700000000	19.0775700000000	18.9753800000000	18.8882900000000	18.7946600000000	18.7127600000000	18.6245600000000	18.5412900000000	18.4116300000000	18.3500500000000	18.2245800000000	18.1433700000000	18.0559300000000	17.9261900000000	17.8456300000000	17.7253800000000	17.6132100000000	17.5096900000000	17.4296700000000	17.3101200000000	17.1856000000000	17.0917100000000	16.9820800000000	16.8312200000000	16.7784500000000	16.6003100000000	16.5309700000000	16.3862200000000	16.3062100000000	16.1836200000000	16.0830900000000	15.9424800000000	15.8502700000000	15.7227800000000	15.6198500000000	15.5014100000000	15.3514100000000	15.2497600000000	15.1558200000000	15.0457100000000	14.9549300000000	14.7757800000000	14.6881700000000	14.5759900000000	14.5104500000000	14.3324600000000	14.2042300000000	14.1239100000000	14.0063300000000	13.8737500000000	13.8115900000000	13.6045700000000	13.5108800000000	13.4538800000000	13.3125200000000	13.2360000000000	13.0378500000000	12.9867300000000	12.8783300000000	12.7316400000000	12.6012000000000	12.5350100000000	12.4032400000000	12.3219200000000	12.1798100000000	12.0879600000000	11.9524800000000	11.8955400000000	11.7372700000000	11.6522000000000	11.5737200000000	11.4736700000000	11.3328200000000	11.2554500000000	11.1808400000000	11.0300000000000	10.9578500000000	10.8436500000000	10.7792200000000	10.6640600000000	10.5591700000000	10.4601800000000	10.4049000000000	10.3006200000000	10.2199700000000	10.1023100000000	9.97630000000000	9.92741000000000	9.83921000000000	9.76028000000000	9.65054000000000	9.63687000000000	9.45379000000000	9.40698000000000	9.34640000000000	9.20574000000000	9.18230000000000	9.11345000000000	9.01409000000000	8.97946000000000	8.84171000000000	8.80488000000000	8.69982000000000	8.67083000000000	8.61643000000000	8.54683000000000	8.46050000000000	8.40984000000000	8.40398000000000	8.27267000000000	8.23544000000000	8.15720000000000	8.14651000000000	8.06167000000000	8.03601000000000	8.00526000000000	7.97113000000000	7.82221000000000	7.84665000000000	7.71988000000000	7.76380000000000	7.67634000000000	7.63108000000000	7.64979000000000	7.57952000000000	7.53109000000000	7.51541000000000	7.54610000000000	7.50653000000000	7.40838000000000	7.38836000000000	7.42465000000000	7.37921000000000	7.37467000000000	7.31236000000000	7.28093000000000	7.30569000000000	7.25641000000000	7.28921000000000	7.26327000000000	7.20379000000000	7.22173000000000	7.20520000000000	7.19907000000000	7.23347000000000	7.20702000000000	7.18582000000000	7.20122000000000	7.21348000000000	7.24831000000000	7.25133000000000	7.14827000000000	7.23198000000000	7.24379000000000	7.23439000000000	7.26383000000000	7.23549000000000	7.24493000000000	7.26955000000000	7.32164000000000	7.34294000000000	7.40255000000000	7.36211000000000	7.39406000000000	7.39935000000000	7.40142000000000	7.45483000000000	7.47392000000000	7.50026000000000	7.53385000000000	7.59808000000000	7.55976000000000	7.64334000000000	7.70204000000000	7.68055000000000	7.69428000000000	7.77867000000000	7.79140000000000	7.84168000000000	7.85676000000000	7.95402000000000	7.98272000000000	7.98250000000000	8.06323000000000	8.06667000000000	8.15448000000000	8.20912000000000	8.27631000000000	8.25311000000000	8.38019000000000	8.34955000000000	8.42581000000000	8.42383000000000	8.53330000000000	8.61176000000000	8.67487000000000	8.66394000000000	8.69607000000000	8.81373000000000	8.83383000000000	8.85312000000000	8.94724000000000	9.01466000000000	9.09353000000000	9.10262000000000	9.20432000000000	9.24938000000000	9.26906000000000	9.36035000000000	9.46434000000000	9.47937000000000	9.62690000000000	9.55900000000000	9.71344000000000	9.79188000000000	9.78627000000000	9.82893000000000	9.94994000000000	9.94472000000000	9.99708000000000	10.1291900000000	10.1100500000000	10.2067800000000	10.2249700000000	10.2924800000000	10.4078800000000	10.4319500000000	10.5427700000000	10.5911300000000	10.7108200000000	10.7125400000000	10.7783300000000	10.8795200000000	10.9173600000000	10.9462500000000	11.0035800000000	11.0962100000000	11.1300600000000	11.2072000000000	11.3370500000000	11.3298000000000	11.4355100000000	11.4856200000000	11.5405500000000	11.6719200000000	11.7671500000000	11.8421200000000	11.8110700000000	11.8141900000000	11.9436900000000	12.0565200000000	12.0409700000000	12.1744300000000	12.1494500000000	12.1945400000000	12.3108800000000	12.3770500000000	12.4227200000000	12.4668700000000	12.5571200000000	12.5809700000000	12.7188700000000	12.8115600000000	12.8060800000000	12.8724300000000	12.8688400000000	12.9987700000000	13.0077700000000	13.0656400000000	13.1497100000000	13.1693300000000	13.2594700000000	13.3319100000000	13.3849700000000	13.4301500000000	13.4843300000000	13.5152300000000	13.5858700000000	13.7076800000000	13.6479400000000	13.7980200000000	13.7183800000000	13.8944400000000	13.8071300000000	13.9629900000000	14.0204600000000	14.0712300000000	14.0673200000000	14.0810000000000	14.0892000000000	14.1761400000000	14.2450700000000	14.3852700000000	14.3196900000000	14.4539400000000	14.4701200000000	14.4509300000000	14.6341500000000	14.6042100000000	14.5936000000000	14.7329500000000	14.7612500000000	14.8288000000000	14.7677700000000	14.8101300000000	14.9051800000000	14.8839500000000	14.9670000000000	14.9969400000000	15.0391800000000	15.0741100000000	15.0789200000000	15.1561200000000	15.1560800000000	15.3090600000000	15.2583500000000	15.2453600000000	15.3450100000000	15.2330500000000	15.3375000000000	15.3079700000000	15.4351500000000	15.4306200000000	15.4678000000000	15.5179100000000	15.4706000000000	15.5276800000000	15.5520000000000	15.6001200000000	15.5905700000000	15.6301000000000	15.6234200000000	15.6209300000000	15.7085900000000	15.6866400000000	15.6675900000000	15.7057900000000	15.7564000000000	15.6739400000000	15.8549700000000	15.7195100000000	15.8098300000000	15.8020000000000	15.7826300000000	15.8638700000000	15.8448000000000	15.9611300000000	15.8821100000000	15.8958300000000	15.9140300000000	15.9765800000000	15.8302800000000	15.9454300000000	15.9647900000000	15.9557800000000	15.9320500000000	15.9970700000000	15.9511800000000	15.9371400000000	15.9042100000000	15.9103900000000	15.9646300000000	16.0406500000000	15.8893100000000	15.9603600000000	15.9680800000000	16.0834300000000	16.0422800000000	16.0792600000000	15.9780500000000	15.9696400000000	15.7965200000000	15.8662600000000	15.8302000000000	15.9037500000000	15.8839900000000	15.8179500000000	15.9233000000000	15.7860800000000	15.8369800000000	15.8691100000000	15.8389500000000	15.7590800000000	15.7894200000000	15.8377500000000	15.8310900000000	15.7559500000000	15.7725300000000	15.8664900000000	15.7224600000000	15.7361800000000	15.6749200000000	15.6622200000000	15.6819100000000	15.5943600000000	15.6483400000000	15.5969500000000	15.5604300000000	15.5199500000000	15.4297000000000	15.4898100000000	15.4580500000000	15.5680400000000	15.5683200000000	15.4474100000000	15.2642500000000	15.3774000000000	15.2942400000000	15.1569300000000	15.2281700000000	15.1934300000000	15.3345800000000	15.1647400000000	15.1749200000000	15.0924500000000	15.0476700000000	15.0697100000000	15.0138800000000	14.9964500000000	14.9127900000000	14.9042300000000	14.8997700000000	14.8202700000000	14.8260100000000	14.7640900000000	14.8736400000000	14.6607700000000	14.6225700000000	14.7011300000000	14.6562700000000	14.5390400000000	14.4935300000000	14.5389800000000	14.4106200000000	14.3962100000000	14.2643000000000	14.3794200000000	14.2780600000000	14.2729100000000	14.2924900000000	14.2500700000000	14.1814800000000	14.1577500000000	13.9981500000000	13.9386700000000	13.9056700000000	13.9939000000000	13.8954800000000	13.8874500000000	13.7490800000000	13.6851900000000	13.7565000000000	13.6985700000000	13.6761000000000	13.4983300000000	13.4814400000000	13.3934900000000	13.3821700000000	13.2026000000000	13.2167700000000	13.2068300000000	13.1913900000000	13.1761400000000	13.1000400000000	13.0414300000000	13.0738400000000	12.9611400000000	12.8909200000000	12.8839100000000	12.6457200000000	12.7179700000000	12.6783100000000	12.7408700000000	12.4990600000000	12.6374000000000	12.5050300000000	12.5102000000000	12.3398800000000	12.3381500000000	12.3132400000000	12.2057200000000	12.1284100000000	12.1292100000000	12.0656300000000	12.0328400000000	12.1187000000000	11.9803600000000	11.8189100000000	11.8508300000000	11.7261700000000	11.6002200000000	11.5926000000000	11.6188900000000	11.4889800000000	11.5983200000000	11.4578300000000	11.3212100000000	11.3506200000000	11.1753400000000	11.2287400000000	11.2576900000000	11.2271600000000	11.0962900000000	10.9437400000000	10.9918200000000	10.8666300000000	10.8369800000000	10.7823400000000	10.6788300000000	10.5431700000000	10.6596200000000	10.5229900000000	10.5552900000000	10.4711600000000	10.4059300000000	10.2921600000000	10.2287900000000	10.2648400000000	10.1060400000000	10.1803600000000	9.96138000000000	9.88087000000000	9.89172000000000	10.0090900000000	9.78840000000000	9.78652000000000	9.58635000000000	9.48667000000000	9.61598000000000	9.44767000000000	9.49498000000000	9.45597000000000	9.40811000000000	9.18694000000000	9.29270000000000	9.41237000000000	9.06259000000000	9.09988000000000	8.98734000000000	8.91117000000000	8.78318000000000	8.86076000000000	8.85996000000000	8.83512000000000	8.65596000000000	8.64555000000000	8.49175000000000	8.56205000000000	8.53196000000000	8.30412000000000	8.43654000000000	8.18189000000000	8.26469000000000	8.16274000000000	8.18542000000000	8.02881000000000	8.09038000000000	7.95267000000000	7.86746000000000	7.94185000000000	7.80496000000000	7.73175000000000	7.61856000000000	7.60907000000000	7.69778000000000	7.54911000000000	7.50830000000000	7.50622000000000	7.41533000000000	7.25051000000000	7.31388000000000	7.18978000000000	7.22048000000000	7.10399000000000	7.11980000000000	6.92316000000000	7.04799000000000	7.08691000000000	6.92547000000000	6.84555000000000	6.88853000000000	6.78315000000000	6.62203000000000	6.70637000000000	6.65336000000000	6.70351000000000	6.49385000000000	6.33344000000000	6.45291000000000	6.39706000000000	6.32252000000000	6.23981000000000	6.16983000000000	6.29037000000000	6.15019000000000	6.10986000000000	6.25077000000000	6.03984000000000	6.03551000000000	6.02428000000000	5.87559000000000	5.92285000000000	5.88932000000000	5.72588000000000	5.76363000000000	5.72420000000000	5.68451000000000	5.64334000000000	5.62132000000000	5.69612000000000	5.62993000000000	5.53687000000000	5.49209000000000	5.32041000000000	5.38367000000000	5.38297000000000	5.32666000000000	5.34584000000000	5.17947000000000	5.28501000000000	5.31560000000000	5.09062000000000	5.18713000000000	4.94444000000000	5.05718000000000	4.99661000000000	5.08993000000000	4.92249000000000	4.91793000000000	5.01076000000000	5.06538000000000	4.81314000000000	5.03861000000000	4.79817000000000	4.67352000000000	4.61118000000000	4.73800000000000	4.61498000000000	4.63909000000000	4.81342000000000	4.62477000000000	4.56442000000000	4.45838000000000	4.54229000000000	4.47606000000000	4.34804000000000	4.44003000000000	4.61718000000000	4.47283000000000	4.24978000000000	4.36399000000000	4.42482000000000	4.39587000000000	4.38990000000000	4.14732000000000	4.37212000000000	4.32028000000000	4.13452000000000	4.41485000000000	4.07521000000000	4.41365000000000	4.20255000000000	4.27142000000000	4.24217000000000	4.12369000000000	4.24562000000000	4.13837000000000	4.09511000000000	4.15071000000000	4.19257000000000	4.10707000000000	4.05723000000000	4.10150000000000	3.96851000000000	4.04512000000000	4.03759000000000	3.86313000000000	3.83361000000000	4.07813000000000	3.87629000000000	3.67449000000000	3.99722000000000	3.87206000000000	4.03101000000000	3.86122000000000	3.96209000000000	3.85849000000000	4.01981000000000	3.95239000000000	3.98592000000000	3.78436000000000	3.64915000000000	3.88977000000000	3.80074000000000	3.67306000000000	3.76473000000000	3.69148000000000	3.64133000000000	3.92840000000000	3.96332000000000	3.82927000000000	3.81623000000000	3.86069000000000	3.83209000000000	3.83458000000000	3.53912000000000	3.69169000000000	3.74984000000000	3.88119000000000	3.86299000000000	3.88793000000000	3.73023000000000	4.06675000000000	3.82090000000000	3.88531000000000	3.77727000000000	3.87035000000000	3.83032000000000	3.72809000000000	3.63731000000000	3.88678000000000	3.78578000000000	3.81261000000000	3.98652000000000	3.75331000000000	4.00195000000000	3.78431000000000	3.87412000000000	3.71599000000000	4.13462000000000	3.80532000000000	4.10182000000000	3.95715000000000	4.03436000000000	3.88724000000000	4.02672000000000	4.01935000000000	4.01921000000000	3.96372000000000	4.11305000000000	4.11940000000000	4.04311000000000	4.35181000000000	4.13511000000000	4.23236000000000	4.31851000000000	4.25243000000000	4.17402000000000	4.14070000000000	4.22744000000000	4.43210000000000	4.28074000000000	4.42725000000000	4.45216000000000	4.50700000000000	4.14011000000000	4.22407000000000	4.29360000000000	4.34308000000000	4.38938000000000	4.28218000000000	4.59074000000000	4.60133000000000	4.32938000000000	4.71841000000000	4.47698000000000	4.40460000000000	4.60875000000000	4.54999000000000	4.48363000000000	4.53994000000000	4.58667000000000	4.50622000000000	4.56864000000000	4.44814000000000	4.65676000000000	4.64959000000000	4.76146000000000	4.72607000000000	4.79894000000000	4.83301000000000	4.83668000000000	4.86738000000000	5.01546000000000	4.95047000000000	5.04289000000000	4.94468000000000	5.13672000000000	4.89765000000000	5.13028000000000	5.01700000000000	5.10228000000000	5.11716000000000	5.17764000000000	5.05855000000000	5.34030000000000	5.32780000000000	5.08000000000000	5.30428000000000	5.31498000000000	5.35937000000000	5.11362000000000	5.44137000000000	5.21175000000000	5.57371000000000	5.68176000000000	4.99065000000000	5.36265000000000	5.40643000000000	5.82747000000000	5.20752000000000	5.36386000000000	5.68368000000000	5.49109000000000	5.82309000000000	5.86944000000000	5.62570000000000	5.66385000000000	5.72621000000000	6.27788000000000	5.73641000000000	5.79633000000000	6.02302000000000	5.93436000000000	5.67989000000000	6.21493000000000	5.63794000000000	6.36169000000000	5.95980000000000	5.86788000000000	6.34611000000000	5.86577000000000	6.24315000000000	6.50187000000000	6.48968000000000	5.94750000000000	5.83210000000000	6.28200000000000	6.26393000000000	6.36403000000000	6.40709000000000	6.20848000000000	6.38120000000000	6.57575000000000	6.43661000000000	6.50658000000000	6.63227000000000	6.47898000000000	6.70431000000000	6.59091000000000	6.64636000000000	6.60423000000000	6.69779000000000	6.84153000000000	6.59573000000000	6.68461000000000	6.78660000000000	6.78386000000000	6.73120000000000	6.96043000000000	6.87847000000000	6.93157000000000	6.82062000000000	7.12875000000000	6.96306000000000	7.10531000000000	7.08344000000000	7.14288000000000	7.19758000000000	7.21018000000000	7.25861000000000	7.23237000000000	7.18923000000000	7.15536000000000	7.34550000000000	7.31538000000000	7.45249000000000	7.36740000000000	7.41617000000000	7.53259000000000	7.53404000000000	7.47999000000000	7.50309000000000	7.52838000000000	7.58061000000000	7.58316000000000	7.74689000000000	7.66748000000000	7.72963000000000	7.63604000000000	7.69890000000000	7.77640000000000	7.73123000000000	7.82885000000000	7.86312000000000	7.85215000000000	7.95741000000000	7.94975000000000	7.90991000000000	7.96613000000000	8.11376000000000	7.97395000000000	8.05254000000000	8.06058000000000	8.11254000000000	8.10488000000000	8.11502000000000	8.21571000000000	8.25523000000000	8.25329000000000	8.27850000000000	8.32835000000000	8.44962000000000	8.34821000000000	8.33425000000000	8.37276000000000	8.38781000000000	8.45814000000000	8.47408000000000	8.54214000000000	8.55168000000000	8.65214000000000	8.62617000000000	8.62945000000000	8.61778000000000	8.67731000000000	8.66187000000000	8.68404000000000	8.70514000000000	8.80199000000000	8.76613000000000	8.86659000000000	8.86808000000000	8.82583000000000	8.94215000000000	8.92580000000000	8.91337000000000	8.88906000000000	8.96325000000000	9.07433000000000	9.12554000000000	9.01202000000000	9.13008000000000	9.14287000000000	9.13011000000000	9.06824000000000	9.26359000000000	9.17776000000000	9.28026000000000	9.27179000000000	9.34839000000000	9.15067000000000	9.29290000000000	9.32971000000000	9.45525000000000	9.37459000000000	9.42749000000000	9.42442000000000	9.28696000000000	9.41986000000000	9.40623000000000	9.62998000000000	9.54339000000000	9.46335000000000	9.53782000000000	9.62518000000000	9.39550000000000	9.56475000000000	9.55811000000000	9.57936000000000	9.60567000000000	9.50044000000000	9.63257000000000	9.61223000000000	9.73863000000000	9.73938000000000	9.76631000000000	9.77494000000000	9.77942000000000	9.85213000000000	9.77745000000000	9.87319000000000	9.75706000000000	9.85857000000000	9.97425000000000	9.74408000000000	9.91550000000000	9.85888000000000	10.0678300000000	9.98671000000000	9.88626000000000	9.87178000000000	9.92735000000000	10.0666500000000	10.1901300000000	10.0452900000000	9.98188000000000	10.1531700000000	10.0745800000000	10.0170500000000	10.2438200000000	10.2619800000000	10.0695200000000	10.0898400000000	10.1983600000000	10.2937800000000	10.3215100000000	10.2914600000000	10.2266500000000	10.3275200000000	10.2196000000000	10.4844300000000	10.1664900000000	10.3559500000000	10.4161800000000	10.0723400000000	10.2030700000000	10.3393600000000	10.3933500000000	10.2031100000000	10.4279600000000	10.4769100000000	10.3525800000000	10.3496200000000	10.5330800000000	10.5523100000000	10.6902600000000	10.5193400000000	10.4807800000000	10.3527800000000	10.5320800000000	10.5630800000000	10.2028300000000	10.5636700000000	10.3941300000000	10.6328300000000	10.4177300000000	10.6388900000000	10.6964100000000	10.4155600000000	10.3131400000000	10.5327400000000	10.6359500000000	10.6261900000000	10.7764600000000	10.4167200000000	10.4096700000000	10.5768900000000	10.9492500000000	10.7608000000000	10.5275400000000	10.3270700000000	10.5906200000000	10.6402100000000	10.4713800000000	10.6251500000000	10.3216100000000	10.3730700000000];
        if TEMPMOD == 1
            PRBPOWER = (BStrans(PRBWL)./100).*(BSref(PRBWL)./100).*PRBPOWER;
        end

        % Correct powers
        SIGNAL = SIGNAL./IRPOWER;
        SIGSDS = SIGSDS./IRPOWER;
        if POWERNORMTYPE > 1 % 2 = IR and probe
            SIGNAL = SIGNAL./PRBPOWER;
            SIGSDS = SIGSDS./PRBPOWER;
            if POWERNORMTYPE > 2 % 3 = PMT gain correction
                gaincoef1 = [0.0327055169810707	-0.0126807096058019]; gaincoef2 = [-1.43914965365560	-0.0691139784840697];
                gaincurve1 = polyval(gaincoef1,1:20); gaincurve2 = gaincoef2(1)*exp(gaincoef2(2)*(21:100))+1;
                gaincurve = [gaincurve1 gaincurve2]; gaincurve = gaincurve.'/gaincurve(1); % scaled to PMT gain 1
                gainfactor = gaincurve(PMTGAIN); SIGNAL = SIGNAL/gainfactor; SIGSDS = SIGSDS/gainfactor;
            end
        end
    end
    return
end

% Baseline fitting
function [BASECURVE,BLEACHRATE,BASEFIT,LBBASE,UBBASE,SBR,CORRSIG,CORRSIGBASE,TBASE,SIGBASE,CURRBASEFITTYPE] = ...
    basefitfn(T,SIGNAL,CURRBASEFITTYPE,CUTLOW,CUTHIGH,TROUBLESHOOT)
%%% This function performs baseline fitting and can guess which function
% to use based on the variation of the data.

    % Setup baseline fitting
    ilow = ceil(length(T)*CUTLOW); % using ceiling to round up
    ihigh = ceil(length(T)*CUTHIGH);
    if isempty(CURRBASEFITTYPE)
        frontmean = mean(SIGNAL(1:ilow));
        tailmean = mean(SIGNAL(ihigh:end));
        fronttaildiff = 2*(frontmean-tailmean)/(frontmean+tailmean);
        if fronttaildiff > 0.03
            CURRBASEFITTYPE = 3;
        else
            CURRBASEFITTYPE = 1;
        end
    end
    if CURRBASEFITTYPE >= 2 % if exponential baseline
        startbase = 1; % start at first point
    else % if linear baseline
        startbase = 2; % start at second point
    end
    TBASE = [T(startbase:ilow); T(ihigh:end)]; % trimmed t
    SIGBASE = [SIGNAL(startbase:ilow); SIGNAL(ihigh:end)]; % trimmed signal
    SBR = max(SIGNAL(ilow:ihigh))/mean(SIGNAL(ihigh:end));
    % Baseline fitting
    if CURRBASEFITTYPE == 0 % no baseline fitting, set basecurve to 0
        BASECURVE = zeros(height(t),width(t));
        BLEACHRATE = 0;
        BASEFIT = zeros(1,5); LBBASE = BASEFIT; UBBASE = BASEFIT;
    else
        basefn = @(b) b(1)+b(2)*exp(-b(3)*(TBASE-b(4)))+b(5)*TBASE-SIGBASE; % error function
        LBBASE = [-Inf -Inf 0 -Inf -Inf];
        UBBASE = [Inf Inf Inf Inf Inf];
        if CURRBASEFITTYPE == 1
            LBBASE(2:4) = 0; UBBASE(2:4) = 0;
        end
        if CURRBASEFITTYPE == 2
            LBBASE(5) = 0; UBBASE(5) = 0;
        end
        baseopt=optimoptions(@lsqnonlin);
        baseopt.MaxFunctionEvaluations = 1e6; % let the fit run longer
        baseopt.MaxIterations = 1e6; % let the fit run longer
        baseopt.FunctionTolerance = 1e-12; % make the fit more accurate
        baseopt.OptimalityTolerance = 1e-12; % make the fit more accurate
        baseopt.StepTolerance = 1e-12; % make the fit more precise
        if TROUBLESHOOT == 0
            baseopt.Display = 'off'; % silence console output
        end
        basegs = [min(SIGNAL) 0.1*(max(SIGNAL)-min(SIGNAL)) 0.1 0 0];
        BASEFIT = lsqnonlin(basefn,basegs,LBBASE,UBBASE,baseopt);
        BASECURVE = BASEFIT(1)+BASEFIT(2)*exp(-BASEFIT(3)*(T-BASEFIT(4)))+BASEFIT(5)*T;
        if BASEFIT(2) > 0 && BASEFIT(3) > 0 
            BLEACHRATE = BASEFIT(3); % ps-1
        else % bleach rate doesn't have physical meaning
            BLEACHRATE = 0;
        end
        if TROUBLESHOOT == 1 % Check baseline fitting
            figure; plot(T,SIGNAL,'o',T,BASECURVE,'-'); title('Baseline fitting');
        end
    end
    CORRSIG = SIGNAL-BASECURVE;
    CORRSIGBASE = [CORRSIG(2:ilow); CORRSIG(ihigh:end)];
    return
end

% Summary values
function [SIGPEAK,PEAKSIGN,PEAKSD,INTEGSIG,INTEGSD,NOISE,SNR,BASESD,SNRV2,SNRSWEEP,INTSNR] = ...
    summaryvals(CORRSIG,CORRSIGBASE,SIGSDS)
%%% This function calculates summary statistics/values for data after 
% baseline correction.

    % Determine if signal peak is negative or positive
    abscorrsig = abs(CORRSIG);
    flipcorrsig = -1*CORRSIG;
    corrabsdist = sum((CORRSIG-abscorrsig).^2);
    flipabsdist = sum((flipcorrsig-abscorrsig).^2);
    if corrabsdist < flipabsdist
        [SIGPEAK,maxindex] = max(CORRSIG); % corrsig is closer to abscorrsig
        PEAKSIGN = 1; % peak is positive
    else
        [SIGPEAK,maxindex] = min(CORRSIG); % flipcorrsig is closer to abscorrsig
        PEAKSIGN = -1; % peak is negative
    end
    % Calculate summary values
    PEAKSD = SIGSDS(maxindex);
    INTEGSIG = sum(CORRSIG);
    INTEGSD = mean(SIGSDS,'all');
    NOISE = range(CORRSIGBASE);
    BASESD = std(CORRSIGBASE);
    SNR = SIGPEAK/NOISE;
    SNRV2 = SIGPEAK/BASESD;
    SNRSWEEP = CORRSIG/NOISE;
    INTSNR = sum(SNRSWEEP);
    return
end

% Lifetime fitting
function [TINT,SIGINT,FITVECTOR,TFIT,FITCURVE,FITVAL,LIFETIME1,LIFETIME2,LT1LT2,FWHM,R2,RESID,SRR,FFTX,FFTY,CURRLTFITTYPE] = ...
    ltfitfn(T,SIGNAL,CURRLTFITTYPE,PEAKSIGN,SIGPEAK,SNR,CURRPULSEWIDTH,LTMIN,LTMAX,BASEFIT,CURRBASEFITTYPE,FLOATBASE,IRWNum,PROBE,TROUBLESHOOT)
%%% This function is performs lifetime fitting (convolution of Gaussian
% with vibrational decay) on an input signal vector for a given time
% vector. It can guess what type of fit to use based on frequency.

    if isempty(CURRLTFITTYPE) % auto-choose lifetime fitting
        if IRWNum > 1680 && IRWNum < 1800 % carbonyl
            CURRLTFITTYPE = 1; % monoexponential
        else % not carbonyl
            if IRWNum < 2400 && IRWNum > 2000 % triple bond
                if PROBE < 520
                    CURRLTFITTYPE = 2; % biexponential (Coumarin 337 nitrile)
                else
                    CURRLTFITTYPE = 1; % monoexponential (normal)
                end
            else % double bond or CH
                if (IRWNum == 1300 || IRWNum == 1500 || IRWNum == 1590) && PROBE < 770
                    CURRLTFITTYPE = 4; % stretched/compressed biexponential
                else
                    CURRLTFITTYPE = 2; % biexponential
                end
            end
        end
    end
    % Lifetime fitting
    if CURRLTFITTYPE > 0
        % Prepare for lifetime fitting
        deltat = zeros(length(T)-1,1); % calculate t spacing
        for i=1:length(T)-1
            deltat(i) = round(1e2*(T(i+1)-T(i)))/1e2; % round to nearest 0.01 ps (had issues without rounding)
        end
        tspacing = unique(deltat); % length 1 <-> evenly spaced t
        % Check t spacing (even) and length (odd)
        if isscalar(tspacing) && mod(length(T),2) == 1
            TINT = T; SIGINT = SIGNAL; % data are good for fitting
        else % interpolate for convolution
            TINT = linspace(T(1),T(end),length(T)*2-1).';
            SIGINT = spline(T, SIGNAL, TINT);
            if TROUBLESHOOT == 1
                figure; plot(TINT,SIGINT,'o',T,SIGNAL,'o'); title('Interpolation')
            end
        end
        if max(TINT) < 100 % auto-pad signal if tint doesn't go past 100 ps
            padsignal = 1;
        else
            padsignal = 0;
        end
        % Prepare for signal padding
        tintspace = TINT(2)-TINT(1);
        originallength = length(TINT); % should be an odd number
        if mod(originallength-1,4) == 0
            padlength = round(originallength/2)+1;
        else
            padlength = round(originallength/2);
        end
        tlong = linspace(TINT(1),TINT(end)+tintspace*padlength,length(TINT)+padlength).';
        baselong = BASEFIT(1)+BASEFIT(2)*exp(-BASEFIT(3)*(tlong-BASEFIT(4)))+...
            BASEFIT(5)*tlong;
        if CURRBASEFITTYPE > 0
            siglong = baselong;
            siglong(1:length(SIGINT)) = SIGINT;
        else
            if PEAKSIGN == 1
                siglong(length(SIGINT)+1:length(tlong)) = min(SIGINT(end-5:end));
            else
                siglong(length(SIGINT)+1:length(tlong)) = max(SIGINT(end-5:end));
            end
        end
        % Signal padding for short-tailed data
        if padsignal == 1
            TINT = tlong;
            SIGINT = siglong;
            if TROUBLESHOOT == 1
                figure; plot(TINT,SIGINT,'-o',T,SIGNAL,'o'); title('Padding');
            end
        end
        % Prepare convolution input (same spacing, half length, symmetrically around zero)
        tc = linspace(0.5*min(TINT),0.5*max(TINT),0.5*(length(TINT)+1)).';
        fun = @(r) r(1)+r(2)*exp(-r(3)*(TINT-r(4)))+r(5)*TINT+ ... % baseline
            conv(exp(-((tc-r(6))./(r(7)/(2*sqrt(log(2))))).^2), ...
            heaviside(tc).*(r(8)*exp(-(tc)./r(9))+r(10)*exp(-((tc)./(r(11))).^r(12)))) - SIGINT;
        % Initial guesses
        r0 = [BASEFIT,0,2.6,... % basecoefs, IRF center (ps), IRF width (ps)
            1,1,0.02,10,1]; % amp 1, τ1 (ps), amp 2, τ2 (ps), β (stretch)
        % Float baseline coeffs
        newlbbase = BASEFIT; newubbase = BASEFIT;
        for i=1:length(BASEFIT)
            if abs(BASEFIT(i)) < 1e-6 % set very small values to zero
                newlbbase(i) = 0;
                newubbase(i) = 0;
            else
                newlbbase(i) = min((1-FLOATBASE)*BASEFIT(i),(1+FLOATBASE)*BASEFIT(i));
                newubbase(i) = max((1-FLOATBASE)*BASEFIT(i),(1+FLOATBASE)*BASEFIT(i));
            end
        end
        % Lower and upper bounds
        lb = [newlbbase,-10,1*sqrt(2),... % basecoefs, IRF center (ps), IRF min (ps)
            0,0.1,... % amp 1, τ1min (ps)
            0,0.1,... % amp 2, τ2min (ps)
            1e-2]; % β min
        ub = [newubbase,20,3*sqrt(2),... % basecoefs, IRF center (ps), IRF max (ps)
            10*SIGPEAK,100,... % amp 1, τ1max (ps)
            10*SIGPEAK,100,... % amp 2, τ2max (ps)
            1e2]; % β max
        % Implement fit type options
        if CURRLTFITTYPE == 1 % Gauss*monoexp
            lb(10) = 0; ub(10) = 0; % first exp only
        end
        if CURRLTFITTYPE == 2 % Gauss*biexp
            lb(12) = 1; ub(12) = 1; % no stretching
        end
        if CURRLTFITTYPE == 3 % Gauss*stretchexp
            lb(8) = 0; ub(8) = 0; % second exp only
        end
        % Parameter constraints from config
        if ~isempty(CURRPULSEWIDTH) % pulse width was specified
            lb(7) = CURRPULSEWIDTH*sqrt(2); ub(7) = CURRPULSEWIDTH*sqrt(2);
        end
        if ~isempty(LTMIN) % minimum lifetime was specified
            lb(9) = LTMIN; lb(11) = LTMIN;
        end
        if ~isempty(LTMAX) % maximum lifetime was specified
            ub(9) = LTMAX; ub(11) = LTMAX;
        end
        % Run lsqnonlin - minimizes error function by adjusting r
        options=optimoptions(@lsqnonlin);
        if SNR > 10 % let run longer for high-SNR
            options.MaxFunctionEvaluations = 1e6; % let the fit run longer
            options.MaxIterations = 1e6; % let the fit run longer
            options.FunctionTolerance = 1e-12; % make the fit more accurate
            options.OptimalityTolerance = 1e-12; % make the fit more accurate
            options.StepTolerance = 1e-12; % make the fit more precise
        else
            options.MaxFunctionEvaluations = 1e5; % let the fit run longer
            options.MaxIterations = 1e5; % let the fit run longer
            options.FunctionTolerance = 1e-10; % make the fit more accurate
            options.OptimalityTolerance = 1e-10; % make the fit more accurate
            options.StepTolerance = 1e-10; % make the fit more precise
        end
        if TROUBLESHOOT == 0
            options.Display = 'off'; % silence console output
        end
        FITVAL = lsqnonlin(fun,r0,lb,ub,options); % fit coeffs
        FITVECTOR = FITVAL(1)+FITVAL(2)*exp(-FITVAL(3)*(TINT-FITVAL(4)))+FITVAL(5)*TINT+ ... % baseline
            conv(exp(-((tc-FITVAL(6))./(FITVAL(7)/(2*sqrt(log(2))))).^2), ...
            heaviside(tc).*(FITVAL(8)*exp(-(tc)./FITVAL(9))+FITVAL(10)*exp(-((tc)./(FITVAL(11))).^FITVAL(12))));
        if TROUBLESHOOT == 1
            figure; tiledlayout(2,1); nexttile([1 1]); 
            plot(TINT,SIGINT,'o',TINT,FITVECTOR,'-'); 
            nexttile([1 1]); plot(TINT,SIGINT-FITVECTOR); title('Fitting without crop')
        end
        % Trim padded signal to original length
        TINT = TINT(1:originallength);
        SIGINT = SIGINT(1:originallength);
        FITVECTOR = FITVECTOR(1:originallength);
        % Calculate pulse width from fit
        IRFwidth = FITVAL(7); % ps
        FWHM = IRFwidth/sqrt(2); % ps
        % Biexponential lifetime assignments
        LIFETIME1 = min(FITVAL(9),FITVAL(11)); % ps
        LIFETIME2 = max(FITVAL(9),FITVAL(11)); % ps
        % Lifetime amplitude ratio
        if abs(FITVAL(9)) < abs(FITVAL(11)) % if fitval(9) is shorter (τ1)
            LT1LT2 = FITVAL(8)/FITVAL(10); % A1/A2
        else % means fitval(9) is longer (τ2)
            LT1LT2 = FITVAL(10)/FITVAL(8); % still A1/A2
        end
        if CURRLTFITTYPE == 1 % Gauss*monoexp
            LIFETIME1 = FITVAL(9); LIFETIME2 = 0;
        end

        if CURRLTFITTYPE == 3 % Gauss*stretchexp
            LIFETIME1 = FITVAL(11); LIFETIME2 = 0;
            LT1LT2 = FITVAL(12); % stretching factor
        end
        if CURRLTFITTYPE == 4 % Gauss*strbiexp
            LIFETIME1 = FITVAL(9); LIFETIME2 = FITVAL(11);
            LT1LT2 = FITVAL(8)/FITVAL(10);
        end
        % Calculate residuals and ssresid
        RESID = SIGINT-FITVECTOR;
        SRR = SIGPEAK./std(RESID(2:end));
        ssresid = sum(RESID.^2); % sum of squares of residuals
        tss = sum((SIGINT-mean(SIGINT)).^2); % total sum of squares
        R2 = 1-(ssresid/tss); % coefficient of determination (R^2)
        if R2 < 0 % remove nonsense R^2 values
            R2 = 0;
        end
        % FFT of residuals
        ts = TINT*1e-12; % time, s
        dt = abs(ts(1)-ts(2)); % time step, s
        fs = 1/dt; % freq BW, Hz
        resfft = fft(RESID); % raw FFT
        resfftshift = fftshift(resfft); % centered FFT
        FFTY = abs(resfftshift).^2/length(RESID); % power spectrum
        freqHz = (-length(RESID)/2+1/2:length(RESID)/2-1/2)*(fs/length(RESID)); % FFT frequency vector in Hz
        FFTX = freqHz/29979245800; % FFT freq in cm-1
        if TROUBLESHOOT == 1
            % Check fit
            figure; tiledlayout(3,1); nexttile([1 1]); plot(TINT,RESID); nexttile([2 1]); plot(T,SIGNAL,'o',TINT,FITVECTOR); title('Fitting')
            % Check FFT
            figure; plot(FFTX,FFTY); title('FFT of residuals');
        end
    else % no lifetime fitting
        TINT = T; SIGINT = SIGNAL; FITVECTOR = SIGNAL;
        RESID = SIGNAL; FFTX = T; FFTY = SIGNAL;
        LIFETIME1 = 0; LIFETIME2 = 0; LT1LT2 = 0;
        FITVAL = zeros(12,1); R2 = 0; FWHM = 0;
        originallength = length(T); SRR = 0;
    end
    % Brute-forced smoother fit curve
    TFIT = linspace(TINT(1),TINT(originallength),originallength*10-9).'; % smoother t vector
    FITCURVE = spline(TINT,FITVECTOR,TFIT);
    return
end

% Plotting time-domain spectra
function makesubpanel(T,SIGNAL,SIGSDS,TBASE,SIGBASE,BASECURVE,TFIT,FITCURVE,...
    TINT,RESID,LABEL,LEGEND,CHANNEL)
%%% This function makes subplots for time-domain spectra given raw data and
% fits.

    % Make colors for figure
    plotcolor(1,:) = [0.5 0.5 0.5]; % Data, residuals
    plotcolor(2,:) = [0.9 0.1 0.1]; % CH1 baseline data
    plotcolor(3,:) = [0.7 0.4 0.1]; % CH1 baseline fit
    plotcolor(4,:) = [0.2 0.8 0.8]; % CH2 baseline data
    plotcolor(5,:) = [0.2 0.8 0.5]; % CH2 baseline fit
    plotcolor(6,:) = [0.2 0.2 0.8]; % Lifetime fit
    plotcolor(7,:) = [0 1 1];
    if CHANNEL == 2 % set colors by channel
        col1 = plotcolor(1,:);
        col2 = plotcolor(4,:);
        col3 = plotcolor(5,:);
        col4 = plotcolor(6,:);
    else
        col1 = plotcolor(1,:);
        col2 = plotcolor(2,:);
        col3 = plotcolor(3,:);
        col4 = plotcolor(6,:);
    end
    hold on; xlim([min(T) max(T)]); xlabel('Time delay (ps)'); ylabel(LABEL);
    errorbar(T,SIGNAL,SIGSDS,'o','Color',col1,'LineWidth',2);
    plot(TBASE,SIGBASE,'o','Color',col2,'LineWidth',2); 
    plot(T,BASECURVE,'-','Color',col3,'LineWidth',2);
    if SIGNAL ~= RESID % plot fit if truly fitted (resid = signal if no fit)
        plot(TFIT,FITCURVE,'-','Color',col4,'LineWidth',2);
        LEGEND(4) = {'Fit'}; yyaxis right; % plot residuals on secondary axis
        plot(TINT,RESID); LEGEND(5) = {'Residuals'};
        ylabel('Residuals (AU)');
    end
    legend(LEGEND); lg = legend; lg.EdgeColor = [1 1 1];
    ax = gca; ax.FontSize = 12; hold off;
    return
end

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

% Pairwise statistical comparisons
function [g1mean,g2mean,PVAL,COHEND,COHENLOWER,COHENUPPER] = statcompare(group1,group2)
% This function performs statistical testing to compare two groups
% (unpaired).
    
    % Ensure column vectors
    if length(group1) == width(group1)
        group1 = group1.';
    end
    if length(group2) == width(group2)
        group2 = group2.';
    end
    % Combine into a single array
    if length(group1) > length(group2)
        stack = NaN(length(group1),2);
    else
        stack = NaN(length(group2),2);
    end
    stack(1:length(group1),1) = group1; stack(1:length(group2),2) = group2;
    % One-way ANOVA
    [PVAL,~,stats] = anova1(stack,[],'off');
    g1mean = stats.means(1);
    g2mean = stats.means(2);
    % Effect size
    effsize = meanEffectSize(group1,group2,Effect="cohen",...
        ConfidenceIntervalType="bootstrap", BootstrapOptions=statset...
        (UseParallel=true),NumBootstraps=3000);
    effsizearray = table2array(effsize);
    COHEND = effsizearray(1);
    COHENLOWER = effsizearray(2);
    COHENUPPER = effsizearray(3);
    ann(1) = {"p = "+string(PVAL)};
    ann(2) = {"Cohen's d = "+string(COHEND)+" ["+string(COHENLOWER)+...
        ", "+string(COHENUPPER)+"]"};

    gardnerAltmanPlot(group1,group2,Effect="cohen",...
    ConfidenceIntervalType="bootstrap", ...
    BootstrapOptions=statset(UseParallel=true),NumBootstraps=3000);
    xticklabels({'Cytoplasm','Nucleoplasm','Effect size'})
    yyaxis left; ylabel('Lifetime (ps)'); 
    yyaxis right; ylabel("Cohen's d");
    annotation('textbox',[0.5 0.7 0.2 0.2],'String',ann,...
        'FitBoxToText','on','FontSize',12,'EdgeColor',[1 1 1]);
    title('')
    ax = gca; ax.FontSize = 12;
    return
end

% Write tif (from HW)
function Save_tif(name, A, B, C, D, E)
%%% This function writes up to five matrices into a hyperstack .tif file.
    if ~isempty(A)
        tiffObject = Tiff(name, 'w');
        tagstruct.ImageLength = size(A,1);
        tagstruct.ImageWidth = size(A,2);
        tagstruct.Compression = Tiff.Compression.None;
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagstruct.BitsPerSample = 64;
        tagstruct.SamplesPerPixel = size(A, 3);
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tiffObject.setTag(tagstruct);
        tiffObject.write(A);
        tiffObject.close;
    end
    if ~isempty(B)
        tiffObject = Tiff(name, 'a');
        tagstruct.ImageLength = size(B,1);
        tagstruct.ImageWidth = size(B,2);
        tagstruct.Compression = Tiff.Compression.None;
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagstruct.BitsPerSample = 64;
        tagstruct.SamplesPerPixel = size(B, 3);
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tiffObject.setTag(tagstruct);
        tiffObject.write(B);
        tiffObject.close;
    end
    if ~isempty(C)
        tiffObject = Tiff(name, 'a');
        tagstruct.ImageLength = size(C,1);
        tagstruct.ImageWidth = size(C,2);
        tagstruct.Compression = Tiff.Compression.None;
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagstruct.BitsPerSample = 64;
        tagstruct.SamplesPerPixel = size(C, 3);
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tiffObject.setTag(tagstruct);
        tiffObject.write(C);
        tiffObject.close;
    end
    if ~isempty(D)
        tiffObject = Tiff(name, 'a');
        tagstruct.ImageLength = size(D,1);
        tagstruct.ImageWidth = size(D,2);
        tagstruct.Compression = Tiff.Compression.None;
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagstruct.BitsPerSample = 64;
        tagstruct.SamplesPerPixel = size(D, 3);
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tiffObject.setTag(tagstruct);
        tiffObject.write(D);
        tiffObject.close;
    end
    if ~isempty(E)
        tiffObject = Tiff(name, 'a');
        tagstruct.ImageLength = size(E,1);
        tagstruct.ImageWidth = size(E,2);
        tagstruct.Compression = Tiff.Compression.None;
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagstruct.BitsPerSample = 64;
        tagstruct.SamplesPerPixel = size(E, 3);
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tiffObject.setTag(tagstruct);
        tiffObject.write(E);
        tiffObject.close;
    end
    return
end
