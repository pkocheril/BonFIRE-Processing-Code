%% Batch processing of 2DVF spectra and BonFIRE images 
%%% v2 - improved auto-lifetime fit selection
%%% v3 - bug fixes
%%% v4 - peak fitting in post-batch analysis
%%% v5 - lifetime-weighted IR spectra, bug fixes
%%% v6 - bug fixes, normalization update for improved fitting, tracktiming
%%% v7 - cleaned up plotting
%%% v8 - can partly reprocess summary.xml, improved importing, smarter 
% fitting, now interpolates residuals to original t, added basicltfit 
% function, added pkfit and Voigt functions, updated colors, added 
% monoexponential error estimation
%%% v9 - bug fixes
%%% v10 - added tcompare and normsweep functions to stack time delay sweeps
%%% v11 - added pkfitnt function for BF+NDR-TPA Fano fitting
%%% v12 - fixed SBR calculation and updated .tif saving function

% Initialize
% cd '/Users/pkocheril/Documents/Caltech/WeiLab/Data/2025_03_03_PK/'
clear; clc; close all;

% Main configuration options
loadprevious = []; % [] = auto-detect, 0 = new analysis, 
% 1 = load individual .dats, 2 = read structure.xml, 3 = re-process partial
runtype = 2; % 0 = process all, 1 = test run with a few files,
% 2 = examine a single file, 3 = examine a single image
targetfolders = []; % indices of folders to process, [] = dialog
t0pos = []; % specify t0 position(s) (mm), [] = autofind

% Additional configuration options
ltfittype = []; % [] = auto-choose, 0 = no fitting, 1 = Gaussian*monoexp,
% 2 = Gaussian*biexp, 3 = Gauss*stretchexp, 4 = Gauss*strbiexp
basefittype = []; % [] = auto-choose, 0 = no baseline fit, 1 = linear, 
% 2 = exponential, 3 = exponential+linear
fitchannels = 2; % specify data channels to fit, [] = dialog
writeprocyn = 0; % 1 = write batch processed files, 0 = not
writefigsyn = 0; % 1 = write figure files, 0 = not
powernormtype = 1; % 0 = no normalization, 1 = normalize by IR power,
% 2 = normalize by probe and IR powers, 3 = 2 + PMT gain correction
setpulsewidth = []; % define pulse width (ps) in fit, [] = float
floatbase = 0.1; % fraction to float baseline coeffs (0.1 -> +/- 10%)
cutlow = 0.03; % lower baseline fit cutoff, (default 3rd %ile)
cuthigh = 0.8; % upper baseline fit cutoff, (default 80th %ile)
verbose = 2; % 2 = all info on figure, 1 = regular figure annotation, 
% 0 = no figure annotation
showprogress = 1; % show progress as fitting proceeds over pixels in image
snrcutoff = -Inf; % SNR minimum for lifetime image fitting (-Inf = fit all)
tracktiming = 0; % 0 = default, 1 = track timing of calculations
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

if loadprevious == 3
    summary = readstruct('summarystructure.xml','FileType','xml');
    [indx,~] = listdlg('PromptString',{'Select folders to re-process.'},...
    'SelectionMode','multiple','ListString',subfolders);
    targetfolders = indx;
    loadprevious = 0;
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

if tracktiming == 1
    tic
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
                    fileinfo = split(filename,"_"); % split at "_"s
                    channel = char(fileinfo(end)); % 'CH2'
                    if numel(fileinfo) > 5 % old naming - lots of underscores
                        IRpowermeterstrfull = char(fileinfo(end-1)); IRpowermeterstr = IRpowermeterstrfull(2:end); IRpowermeter = sscanf(IRpowermeterstr,'%f'); % 'P0.4409' --> 0.4409
                        DFGWNstrfull = char(fileinfo(end-2)); DFGWNstr = DFGWNstrfull(4:end); DFGWN = sscanf(DFGWNstr,'%f'); % 'DFG3797.9' --> 3797.9
                        idlerWNstrfull = char(fileinfo(end-3)); idlerWNstr = idlerWNstrfull(6:end); idlerWN = sscanf(idlerWNstr,'%f'); % 'Idler2949.8' --> 2949.8
                        imgdim = char(fileinfo(end-4)); imgdimsplit = split(string(imgdim),"X"); ysteps = imgdimsplit(1); xsteps = imgdimsplit(end); xsteps = sscanf(xsteps,'%f'); ysteps = sscanf(ysteps,'%f'); % '5X5' --> 5 and 5 (Y by X)
                        prbWL = 1; % <-- unknown from filename
                        if IRpowermeter ~= 0
                            IRpower = IRpowermeter*300; % assuming 300 mW scale
                        end
                        % filebase = strjoin(fileinfo(1:end-5),'_');
                    else % new naming - fewer underscores
                        IRpower = 0; DFGWN = 0; idlerWN = 0;
                        % if numel(fileinfo) == 4
                        %     filebase = strjoin(fileinfo(1:end-(numel(fileinfo)-2)),'_');
                        % else
                        %     filebase = strjoin(fileinfo(1:end-(numel(fileinfo)-1)),'_');
                        % end
                    end
                    % Filename from new VI = prefix_other_stuff_mXn_CHx.raw
                    fileinfoX = split(filename,'X'); % split at "X" (for mXn)
                    fileprefix = split(strjoin(fileinfoX(1:end-1),"X"),'_'); % use 1:end-1 in case of multiple "X" in filename
                    filebase = char(strjoin(fileprefix(1:end-1),"_")); % use 1:end-1 in case of multiple "_" in filename
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
                    if length(infofile) < 12 % import failed
                    	infotable = readtable(infofilename); infofile = table2array(infotable(3:end,2));
                    end
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
                        if length(infofile) > 33
                            wheelfilter = infofile(34);
                        end
                        if length(infofile) > 34
                            concentration = infofile(35);
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
                % data = data.*1e3; % convert V to mV
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
                    channels(1) = 1; %ch1data = ch1data.*1e3; % convert V to mV
                end
                ch3file = string(currentfile(1:end-7))+'CH3'+string(currentfile(end-3:end));
                if isfile(ch3file)
                    [~,ch3data,~] = ...
                        loaddata(workingdirectory,subfolders,ii,ch3file,...
                        currfiletype,trimlastT,trimfirstT,xsteps,ysteps,guessTlist);
                    channels(3) = 3; %ch3data = ch3data.*1e3; % convert V to mV
                end
                ch4file = string(currentfile(1:end-7))+'CH4'+string(currentfile(end-3:end));
                if isfile(ch4file)
                    [~,ch4data,~] = ...
                        loaddata(workingdirectory,subfolders,ii,ch4file,...
                        currfiletype,trimlastT,trimfirstT,xsteps,ysteps,guessTlist);
                    channels(4) = 4; %ch4data = ch4data.*1e3; % convert V to mV
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
                    cb.Label.VerticalAlignment = "bottom"; cb.LineWidth = 2; datacursormode; 
                    xlabel('X'); ylabel('Y'); ax = gca; ax.FontSize = 12; ax.LineWidth = 2;
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
                                [~,ch1int,ch1fitvector,~,ch1fitcurve,ch1fitval,...
                                    ch1lifetime1,ch1lifetime2,ch1a1a2,ch1fwhm,ch1r2,ch1resid,ch1srr,ch1fftx,ch1ffty,currltfittype] = ...
                                    ltfitfn(t,ch1,currltfittype,ch1peaksign,ch1peak,ch1snr,...
                                    currpulsewidth,ltmin,ltmax,ch1basefit,currbasefittype,floatbase,IRWN,prbWL,troubleshoot);
                            else % run with currltfittype = 0
                                [~,ch1int,ch1fitvector,~,ch1fitcurve,ch1fitval,...
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
                                [~,ch3int,ch3fitvector,~,ch3fitcurve,ch3fitval,...
                                    ch3lifetime1,ch3lifetime2,ch3a1a2,ch3fwhm,ch3r2,ch3resid,ch3srr,ch3fftx,ch3ffty,currltfittype] = ...
                                    ltfitfn(t,ch3,currltfittype,ch3peaksign,ch3peak,ch3snr,...
                                    currpulsewidth,ltmin,ltmax,ch3basefit,currbasefittype,floatbase,IRWN,prbWL,troubleshoot);
                            else % run with currltfittype = 0
                                [~,ch3int,ch3fitvector,~,ch3fitcurve,ch3fitval,...
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
                                [~,ch4int,ch4fitvector,~,ch4fitcurve,ch4fitval,...
                                    ch4lifetime1,ch4lifetime2,ch4a1a2,ch4fwhm,ch4r2,ch4resid,ch4srr,ch4fftx,ch4ffty,currltfittype] = ...
                                    ltfitfn(t,ch4,currltfittype,ch4peaksign,ch4peak,ch4snr,...
                                    currpulsewidth,ltmin,ltmax,ch4basefit,currbasefittype,floatbase,IRWN,prbWL,troubleshoot);
                            else % run with currltfittype = 0
                                [~,ch4int,ch4fitvector,~,ch4fitcurve,ch4fitval,...
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
                                [~,spcmint,spcmfitvector,~,spcmfitcurve,spcmfitval,...
                                    spcmlifetime1,spcmlifetime2,spcma1a2,spcmfwhm,spcmr2,spcmresid,spcmsrr,spcmfftx,spcmffty,currltfittype] = ...
                                    ltfitfn(t,spcm,currltfittype,spcmpeaksign,spcmpeak,spcmsnr,...
                                    currpulsewidth,ltmin,ltmax,spcmbasefit,currbasefittype,floatbase,IRWN,prbWL,troubleshoot);
                            else % run with currltfittype = 0
                                [~,spcmint,spcmfitvector,~,spcmfitcurve,spcmfitval,...
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
                                    string(prbWL)+' nm, SNR = '+string(snr)+', SNR_{v2} = '+string(snrv2)};
                                if ~isempty(currltfittype)
                                    if currltfittype == 1
                                        annot(length(annot)+1) = {'SRR = '+string(srr)+', τ_{mono} = '+lt1+...
                                            '±'+string(sprintf('%0.1g',ltfiterror(lifetime1,snr)))+' ps, τ_{p} = '+pulsewidth+' ps, r^2 = '+sr2};
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
                                    pmtunitlabel = '(mV/mW_{IR})';
                                    spcmunitlabel = '(cpms/mW_{IR})';
                                else
                                    if powernormtype == 2
                                        pmtunitlabel = '(mV/mW^{2})';
                                    else % powernormtype == 3
                                        pmtunitlabel = '(mV_{gc}/mW^{2})';
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
                                        ch1resid,ch1label,ch1legend,1);
                                else % SPCM
                                    makesubpanel(t,spcm,spcmsds,tbase,spcmbase,spcmbasecurve,tfit,spcmfitcurve,...
                                        spcmresid,spcmlabel,spcmlegend,5);
                                end
                                nexttile([1 1]); xticks([]); yticks([]);
                                printoutdim = [0.67 0.7 0.3 0.25]; nexttile([2 3]);
                            else
                                % CH1-4 and SPCM
                                if length(nonzeros(channels)) > 2
                                    tiledlayout(3,5); nexttile([1 2]); % CH1
                                    makesubpanel(t,ch1,ch1sds,tbase,ch1base,ch1basecurve,tfit,ch1fitcurve,...
                                        ch1resid,ch1label,ch1legend,1);
                                    nexttile([1 1]); xticks([]); yticks([]);
                                    printoutdim = [0.4 0.7 0.2 0.25];
                                    nexttile([1 2]); % SPCM
                                    makesubpanel(t,spcm,spcmsds,tbase,spcmbase,spcmbasecurve,tfit,spcmfitcurve,...
                                        spcmresid,spcmlabel,spcmlegend,5);
                                    nexttile(9,[1 2]); % CH3
                                    makesubpanel(t,ch3,ch3sds,tbase,ch3base,ch3basecurve,tfit,ch3fitcurve,...
                                        ch3resid,ch3label,ch3legend,3);
                                    nexttile(14,[1 2]); % CH4
                                    makesubpanel(t,ch4,ch4sds,tbase,ch4base,ch4basecurve,tfit,ch4fitcurve,...
                                        ch4resid,ch4label,ch4legend,4);
                                    nexttile(6,[2 3]);
                                else % just CH2
                                    tiledlayout(3,2); nexttile([1 2]); xticks([]); yticks([]);
                                    printoutdim = [0.05 0.7 0.9 0.25]; nexttile([2 2]); 
                                end
                            end
                            % Plot CH2
                            makesubpanel(t,signal,sigsds,tbase,sigbase,basecurve,tfit,fitcurve,...
                                resid,ch2label,ch2legend,2);
                            % Annotate with experiment info
                            annotation('rectangle',printoutdim,'Color',[1 1 1],'FaceColor',[1 1 1]); % box to cover
                            annotation('textbox',printoutdim,'String',annot);
                            if writefigsyn == 1
                                saveas(fig,outname+'_proc.png')
                            end
                        end
                        if tracktiming == 1
                            toc; disp('End pixel');
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
                if tracktiming == 1
                    toc; disp('End file');
                end
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
analysistype = []; % 0 = nothing, 1 = load data, 2 = contours, 
% 3 = lifetimes, 4 = peak fitting, 5 = peak overlay, 6 = lifetime-weighted spectra,
% 7 = Tsweep comparison, [] = dialog
savecontours = 0; % 1 = save contours, 0 = not
xaxischoice = 0; % 0 = probeWL, 1 = sumfreq, 2 = detuning (specify),
% 3 = detuning (Rh800 max), 4 = detuning (Rh800 0-0)
logcontour = 0; % 0 = linear scale, 1 = log scale
subset = []; % [] = dialog, targetfolders = run all
plotfits = 0; % 1 = show individual peak fits, else don't
updatesummary = 0; % 1 = update summary structure, 0 = don't

if isempty(analysistype) % dialog to select analysis type
    [indx,~] = listdlg('PromptString',{'Select post-batch analysis type.'},...
    'SelectionMode','multiple','ListString',{'Load data','Contour maps',...
    'Lifetime comparisons','Peak fitting','Peak overlay',...
    'Lifetime-weighted spectra with overlay','Time delay sweep comparison'});
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
IRpowers = prb; prbpowers = prb; noise = prb; dcsbr = prb; modfr = prb;
pkht = NaN(maxfiles,length(subset));
csig = NaN(height(wIR),width(time),length(subset)); rawsig = csig;
csiga = NaN(height(prb),width(timealign),length(subset));
ch1sig = csig;
snr = prb; snrv2 = prb; srr = prb;

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
        modfr(j,i) = summary.(sf{i}).('file'+string(j)).parameters.modfreq;
        ch1sig(j,1:tlength,i) = summary.(sf{i}).('file'+string(j)).ch1data.corrsig;
        snr(j,i) = summary.(sf{i}).('file'+string(j)).summary.snr;
        snrv2(j,i) = summary.(sf{i}).('file'+string(j)).summary.snrv2;
        srr(j,i) = summary.(sf{i}).('file'+string(j)).fit.srr;
    end
end

pkht = abs(pkht);
ltavg = (ltr.*lt1+lt2)./(ltr+1);
photothermal = squeeze(rawsig(:,end,:));
defaultfont = 25;
tD = time(1,:);

if max(ismember(analysistype,2)) % plotting contours
    for i=subset
        nfiles = summary.(sf{i}).nfiles;
        if isscalar(unique(rmmissing(nonzeros(round(prb(:,i)./10))))) % no probe tuning --> IR sweep
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
        if length(rmmissing(w1)) > 1
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

if max(ismember(analysistype,3)) % lifetime comparisons
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
    % Default lifetime comparison
    ltlegend = strings(length(subset),1);
    figure('visible',figvis);
    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
    % τ1
    nexttile([1 2]); hold on;
    for j=1:length(subset)
        i = subset(j);
        [linecolor,marker,marksize] = colormarkset(j,length(subset),[]);
        plot(xdata,lt1(:,i),'.-','Marker',marker,'MarkerSize',marksize,'Color',linecolor,'LineWidth',2);
        ltlegend(i) = string(wIR(1,i))+' cm{-1}';
    end
    xticks([]); ylabel('τ_{1} (ps)'); % no xlabel
    hold off; xlim([min(xdata) max(xdata)]); ylim([0 5]); 
    legend(ltlegend);
    % τ2
    nexttile([1 2]); hold on;
    for j=1:length(subset)
        i = subset(j);
        [linecolor,marker,marksize] = colormarkset(i,length(subset),[]);
        plot(xdata,lt2(:,i),'.-','Marker',marker,'MarkerSize',marksize,'Color',linecolor,'LineWidth',2);
    end
    xlabel(label); ylabel('τ_{2} (ps)');
    hold off; xlim([min(xdata) max(xdata)]); ylim([0 10]); 
    legend(ltlegend);
    % A1/A2
    nexttile([2 1]); hold on;
    for j=1:length(subset)
        i = subset(j);
        [linecolor,marker,marksize] = colormarkset(i,length(subset),[]);
        plot(xdata,ltr(:,i),'.-','Marker',marker,'MarkerSize',marksize,'Color',linecolor,'LineWidth',2);
    end
    % set(gca,'Yscale','log');
    xlabel(label); ylabel('A_{1}/A_{2}');
    xlim([min(xdata) max(xdata)]); ylim([1e-2 1e2]); 
    legend(ltlegend);
end

peakfits.title = 'All fits'; % set up structure

if max(ismember(analysistype,4)) % peak fitting - Gauss2 for probe, Gauss+line for IR
    for i=subset
        nfiles = summary.(sf{i}).nfiles;
        sigmatrix = csig(1:nfiles,:,i);
        if isscalar(unique(rmmissing(nonzeros(round(prb(:,i)./10))))) % no probe tuning --> IR sweep
            xstring = 'ω_{IR} (cm^{-1})';
            w1 = wIR(1:nfiles,i);
            w2 = prb(1:nfiles,i);
            nonswept = string(round(w2(1)))+' nm';
            sumfreq = 1e7./w2+w1;
            % Sorting
            t1 = time(i,:);
            [t1,tind] = sort(t1);
            [w1,wind] = sort(w1);
            sigmatrix = sigmatrix(wind,tind);
            xval = w1; yval = max(sigmatrix,[],2); yval = yval./max(yval);
            % Separate setup
            if min(w1) > 2000 % idler sweep -- Gauss+Fano
                startcolor = [0 0.35 0.5]; % turquoise
                endcolor = [0.8 0.2 0.2]; % red
                [gx,gy,soln,fwhms] = pkfitnt(xval,yval);
                peakfits.('folder'+string(i)).center = soln(5);
                peakfits.('folder'+string(i)).amp = soln(3);
                peakfits.('folder'+string(i)).bftpa = soln(3)./soln(2);
                % % Check Fano fitting
                % figure; hold on; plot(xval,yval,'o'); plot(gx,soln(1)+soln(2)*gx,'--')
                % plot(gx,soln(1)+soln(2)*gx+soln(6).*((soln(7).*soln(8)+gx-soln(5)).^2)./(soln(8).^2+(gx-soln(5)).^2),'--'); plot(gx,gy,'--');
            else % DFG sweep
                startcolor = [0.5 0 0.5]; % purple
                endcolor = [0.2 0.8 0.2]; % green
                npeaks = 4;
                % Fitting
                [gx,gy,soln,fwhms] = pkfit(xval,yval,npeaks,4,'positive'); % DFG sweep -- 4Gauss+polynomial
                peakfits.('folder'+string(i)).center = soln(8);
                peakfits.('folder'+string(i)).amp = soln(6);
            end
            peakfits.('folder'+string(i)).fwhms = fwhms;
            peakfits.('folder'+string(i)).fitvector = soln;
        else % probe sweep
            xstring = 'ω_{probe} (cm^{-1})';
            w1 = 1e7./prb(1:nfiles,i);
            w2 = wIR(1:nfiles,i);
            nonswept = string(round(w2(1)))+' cm^{-1}';
            sumfreq = w1+w2;
            % Color setup
            startcolor = [0 0 1]; % blue
            endcolor = [1 0.5 0]; % orange
            % Sorting
            t1 = time(i,:);
            [t1,tind] = sort(t1);
            [w1,wind] = sort(w1);
            sigmatrix = sigmatrix(wind,tind);
            % Fitting
            xval = w1; yval = max(sigmatrix,[],2); yval = yval./max(yval);
            [gx,gy,soln,fwhms] = pkfit(xval,yval,2,1,'positive'); % 2 Gaussians, linear background
            peakfits.('folder'+string(i)).center1 = soln(8);
            peakfits.('folder'+string(i)).center2 = soln(11);
            peakfits.('folder'+string(i)).amp1 = soln(6);
            peakfits.('folder'+string(i)).amp2 = soln(9);
            peakfits.('folder'+string(i)).fwhm1 = fwhms(1);
            peakfits.('folder'+string(i)).fwhm2 = fwhms(2);
            peakfits.('folder'+string(i)).fitvector = soln;
        end
        if plotfits == 1
            figure; hold on; box on;
            plot(xval,yval,'.','MarkerSize',20,'Color',startcolor);
            plot(gx,gy,'LineWidth',3,'Color',endcolor);
            xlabel(xstring); ylabel('Relative peak (AU)');
            xlim([min(w1) max(w1)]);
            ax = gca; ax.FontSize = 12; ax.LineWidth = 2;
        end
        peakfits.('folder'+string(i)).foldername = string(subfolders(targetfolders(i)));
        peakfits.('folder'+string(i)).xval = xval;
        peakfits.('folder'+string(i)).yval = yval;
        peakfits.('folder'+string(i)).xfit = gx;
        peakfits.('folder'+string(i)).yfit = gy;
        peakfits.('folder'+string(i)).xstring = xstring;
    end
end

if max(ismember(analysistype,5)) % peak overlay
    [overset,~] = listdlg('PromptString',{'Select folders to overlay.'},...
    'SelectionMode','multiple','ListString',subfolders(targetfolders));
    overleg = subfolders(targetfolders(overset));
    figure; hold on; box on; % overlay figure

    % Lifetime-weighted spectra
    if max(ismember(analysistype,6))
        tiledlayout(3,1); nexttile([2 1]); hold on;
        % Make legend
        for j=1:length(overset)
            i = overset(j);
            % Decide color and markers
            [linecolor,marker,marksize] = colormarkset(j,length(overset),[]);
            plot([-1 -1],[-1 -1],'-','Marker',marker,'MarkerSize',marksize,'LineWidth',2,'Color',linecolor);
        end
        legend(overleg); lg = legend; lg.EdgeColor = [1 1 1]; lg.AutoUpdate = 'off';
        % Plot data
        for j=1:length(overset)
            i = overset(j);
            % Decide color and markers
            [linecolor,marker,marksize] = colormarkset(j,length(overset),[]);
            ltwt = lt1(:,i).*pkht(:,i);
            [maxwt,maxltind] = max(ltwt);
            peakfits.('folder'+string(i)).ltweighted.maxlt = lt1(maxltind,i); % lifetime at peak of lt-weighted sweep
            peakfits.('folder'+string(i)).ltweighted.maxsrr = srr(maxltind,i); % SRR at peak of lt-weighted sweep
            % Fitting
            xval = wIR(:,i); yval = ltwt./maxwt;
            [gx,gy,soln,fwhms] = pkfit(xval,yval,1,2,'positive');
            peakfits.('folder'+string(i)).ltweighted.xval = xval;
            peakfits.('folder'+string(i)).ltweighted.yval = yval;
            peakfits.('folder'+string(i)).ltweighted.center = soln(8);
            peakfits.('folder'+string(i)).ltweighted.amp = soln(6);
            peakfits.('folder'+string(i)).ltweighted.fwhms = fwhms;
            peakfits.('folder'+string(i)).ltweighted.fitvector = soln;
            peakfits.('folder'+string(i)).ltweighted.xfit = gx;
            peakfits.('folder'+string(i)).ltweighted.yfit = gy;
            % Plotting
            plot(xval,yval./max(gy),marker,'MarkerSize',marksize,'LineWidth',2,'Color',linecolor);
            plot(gx,gy./max(gy),'LineWidth',3,'Color',linecolor);
        end
        xlim([2100 2300]); ylim([0 1.2]);
        xlabel('ω_{IR} (cm^{-1})'); ylabel('Lifetime-weighted spectra (AU)');
        ax = gca; ax.FontSize = 15;
        nexttile([1 1]); hold on;
    end
    
    % Regular spectra
    for j=1:length(overset)
        % Make legend
        i = overset(j);
        % Decide color and markers
        [linecolor,marker,marksize] = colormarkset(j,length(overset),[]);
        plot([-1 -1],[-1 -1],'-','Marker',marker,'MarkerSize',marksize,'LineWidth',2,'Color',linecolor);
    end
    xlabel(peakfits.('folder'+string(i)).xstring); ylabel('Normalized peak (AU)');
    if max(ismember(analysistype,6)) == 0 % make legend if not made already
        legend(overleg); lg = legend; lg.EdgeColor = [1 1 1]; lg.AutoUpdate = 'off';
    end
    ax = gca; ax.FontSize = 15; xl1 = []; xl2 = [];
    for j=1:length(overset)
        i = overset(j);
        % Decide color and markers
        [linecolor,marker,marksize] = colormarkset(j,length(overset),[]);
        % Plot data
        xval = peakfits.('folder'+string(i)).xval;
        yval = peakfits.('folder'+string(i)).yval;
        gx = peakfits.('folder'+string(i)).xfit;
        gy = peakfits.('folder'+string(i)).yfit;
        plot(xval,yval./max(gy),marker,'MarkerSize',marksize,'LineWidth',2,'Color',linecolor);
        % Plot fits
        plot(gx,gy./max(gy),'LineWidth',3,'Color',linecolor);
        % Figure out bounds
        if isempty(xl1) || min(xval) < xl1
            xl1 = min(xval);
        end
        if isempty(xl2) || max(xval) > xl2
            xl2 = max(xval);
        end
    end
    xlim([xl1 xl2]); ylim([0 Inf]);
end

if max(ismember(analysistype,7)) % Tsweep comparison
    for i=1:length(subset)
        j = targetfolders(subset(i));
        nfiles = summary.('folder'+string(j)).nfiles;
        sigmat = squeeze(csig(1:nfiles,:,i));
        if i == 3 || i == 5 || i == 10 % mod freq sweep
            [sortmod,sortinds] = sort(modfr(1:nfiles,i));
            leg = string(round(sortmod)./1000)+' kHz';
            colormap = 'rainbow';
        else % power sweep
            [sortmod,sortinds] = sort(IRpowers(1:nfiles,i));
            leg = string(round(10*sortmod)./10)+' mW';
            colormap = [];
        end
        timelength = length(rmmissing(time(i,:)));
        tcompare(time(i,1:timelength),sigmat(sortinds,1:timelength),leg,colormap,'normalize');
        title(summary.('folder'+string(j)).foldername);
    end
end

if updatesummary == 1
    summary.peakfits = peakfits; writestruct(summary,'summarystructure.xml',FileType="xml");
end


%% Pairwise comparisons
statcomparison = 0;
if statcomparison == 1
    clear; close all; clc;
    cd '/Users/pkocheril/Documents/Caltech/WeiLab/Data/2024_03_12/'
    % Load lifetime images (pre-masked)
    cond1 = double(tiffreadVolume("03. 80 uM Rh800 LLPS/FOV1_condens_LT.tif"));
    buff1 = double(tiffreadVolume("03. 80 uM Rh800 LLPS/FOV1_buffer_LT.tif"));
    cond2 = double(tiffreadVolume("03. 80 uM Rh800 LLPS/FOV2_condens_LT.tif"));
    buff2 = double(tiffreadVolume("03. 80 uM Rh800 LLPS/FOV2_buffer_LT.tif"));
    
    condensate = [nonzeros(cond1); nonzeros(cond2)];
    buffer = [nonzeros(buff1); nonzeros(buff2)];
    name2 = 'Condensate';
    name1 = 'Buffer';
    testtype = 'mean'; % mean or median
    % Run comparison and make summary plot
    figure; box on;
    [buffermean,condmean,pval,cohend,cohenlower,cohenupper] = ...
        statcompare(buffer,condensate,name1,name2,testtype);
    bufferstd = std(buffer);
    condstd = std(condensate);
    ax = gca; ax.FontSize = 12; ax.LineWidth = 2;

    % Get y values from histograms
    binstep = 0.1;
    binedg = min(buffer):binstep:max(buffer);
    figure('visible','off'); histogram(buffer,'BinEdges',binedg);
    ax = gca; histhandle = ax.Children; bufferdata = histhandle.Values;
    figure('visible','off'); histogram(condensate);
    ax = gca; histhandle = ax.Children; conddata = histhandle.Values;

    t = binedg(1:end-1); signal = bufferdata;
    [tfit,fitcurve,~,~] = basicltfit(t,signal,2);
    figure; hold on; box on;
    histogram(buffer,'BinEdges',binedg,'FaceColor',[40 189 189]./255);
    histogram(condensate,'BinEdges',binedg,'FaceColor',[180 93 226]./255);
    % plot(t,signal,'ko','LineWidth',2);
    plot(tfit+binstep./2,fitcurve,'-','LineWidth',3,'Color',[40 189 189]./255);
    t = binedg(1:end-1); signal = conddata;
    [tfit,fitcurve,fitval,resid] = basicltfit(t,signal,2);
    plot(tfit+binstep./2,fitcurve,'-','LineWidth',3,'Color',[180 93 226]./255);
    xlim([min(binedg) max(binedg)]); ylim([0 Inf]);
    xlabel('Lifetime (ps)'); ylabel('Counts');
    ax = gca; ax.FontSize = 12; ax.LineWidth = 2;

    % Example Voigt plotting
    [voigty,voigtx] = Voigt(1,1,1,1,1:0.1:500);
end


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
        SIGNAL = 1e3*SIGNAL; SIGSDS = 1e3*SIGSDS; % convert V to mV
        % Power corrections
        if PINHOLEYN == 1
            PRBPOWER = 0.27*PRBPOWER; % pinhole transmission is 27%
        end
        PRBPOWER = PRBPOWER/(10^(NDF)); % power loss from ND filter
        PRBWL = round(PRBWL);
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
        IRPOWER = max([IRPOWER 1]); % don't divide by IR power less than 1 mW
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
    SBR = (max(SIGNAL(ilow:ihigh))-mean(SIGNAL(ihigh:end)))/mean(SIGNAL(ihigh:end));
    % Baseline fitting
    if CURRBASEFITTYPE == 0 % no baseline fitting, set basecurve to 0
        BASECURVE = zeros(height(T),width(T));
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
        baseopt.FunctionTolerance = abs(max(SIGNAL))*1e-12; % make the fit more accurate
        baseopt.OptimalityTolerance = abs(max(SIGNAL))*1e-12; % make the fit more accurate
        baseopt.StepTolerance = abs(max(SIGNAL))*1e-12; % make the fit more precise
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
        if IRWNum > 1702 && IRWNum < 1800 % carbonyl
            CURRLTFITTYPE = 1; % monoexponential
        else % not carbonyl
            if IRWNum < 2400 && IRWNum > 2000 % triple bond
                if PROBE < 520 && PROBE > 390 % 515.6 nm or 400 nm
                    CURRLTFITTYPE = 2; % biexponential (Coumarin 337 nitrile or small cyanocoumarin)
                else
                    CURRLTFITTYPE = 1; % monoexponential (normal probe)
                end
            else % double bond or CH
                if (IRWNum == 1300 || IRWNum == 1500 || IRWNum == 1590) && PROBE < 770 && PROBE > 690
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
            4.5,1,1.5,10,1]; % amp 1, τ1 (ps), amp 2, τ2 (ps), β (stretch)
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
            0.1]; % β min
        ub = [newubbase,20,3*sqrt(2),... % basecoefs, IRF center (ps), IRF max (ps)
            10*SIGPEAK,100,... % amp 1, τ1max (ps)
            10*SIGPEAK,100,... % amp 2, τ2max (ps)
            1.8]; % β max
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
            options.FunctionTolerance = abs(max(SIGINT)-min(SIGINT))*1e-12; % make the fit more accurate
            options.OptimalityTolerance = abs(max(SIGINT)-min(SIGINT))*1e-12; % make the fit more accurate
            options.StepTolerance = abs(max(SIGINT)-min(SIGINT))*1e-12; % make the fit more precise
        else
            options.MaxFunctionEvaluations = 1e5; % let the fit run longer
            options.MaxIterations = 1e5; % let the fit run longer
            options.FunctionTolerance = abs(SIGPEAK)*1e-10; % make the fit more accurate
            options.OptimalityTolerance = abs(SIGPEAK)*1e-10; % make the fit more accurate
            options.StepTolerance = abs(SIGPEAK)*1e-10; % make the fit more precise
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
        RESID = spline(TINT,RESID,T);
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

% Basic lifetime fitting
function [TFIT,FITCURVE,FITVAL,RESID] = basicltfit(T,SIGNAL,CURRLTFITTYPE)
%%% This function performs lifetime fitting but without the extra "smarts"
% of the previous function.
    
    [sigmax,maxind] = max(SIGNAL);
    % Ensure column vectors
    if length(T) == width(T)
        T = T.';
    end
    if length(SIGNAL) == width(SIGNAL)
        SIGNAL = SIGNAL.';
    end
    % Make sure lengths match
    if length(T) ~= length(SIGNAL)
        if length(T) < length(SIGNAL) % SIGNAL is longer
            SIGNAL = SIGNAL(1:length(T));
        else % T is longer
            T = T(1:length(SIGNAL));
        end
    end
    
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
    siglong = SIGINT(end).*ones(length(tlong),1);
    siglong(1:length(SIGINT)) = SIGINT;
    
    % Signal padding for short-tailed data
    TINT = tlong;
    SIGINT = siglong;

    % Prepare convolution input (same spacing, half length, symmetrically around zero)
    tc = linspace(0.5*min(TINT),0.5*max(TINT),0.5*(length(TINT)+1)).';
    fun = @(r) r(1)+r(2)*exp(-r(3)*(TINT-r(4)))+r(5)*TINT+ ... % baseline
        conv(exp(-((tc-r(6))./(r(7)/(2*sqrt(log(2))))).^2), ...
        heaviside(tc).*(r(8)*exp(-(tc)./r(9))+r(10)*exp(-((tc)./(r(11))).^r(12)))) - SIGINT;
    % Initial guesses
    r0 = [0,0,0,0,0,... % basecoefs
        T(maxind),2.6,... % IRF center (ps), IRF width (ps)
        1,1,0.02,10,1]; % amp 1, τ1 (ps), amp 2, τ2 (ps), β (stretch)
    % Lower and upper bounds
    lb = [-Inf,-Inf,0,-Inf,-Inf,... % basecoefs
        min(T),0.1,... % IRF center (ps), IRF min (ps)
        0,0.1,... % amp 1, τ1min (ps)
        0,0.1,... % amp 2, τ2min (ps)
        1e-2]; % β min
    ub = [Inf,Inf,0,Inf,Inf,... % basecoefs, 
        max(T),3*sqrt(2),... % IRF center (ps), IRF max (ps)
        10*sigmax,100,... % amp 1, τ1max (ps)
        10*sigmax,100,... % amp 2, τ2max (ps)
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
    
    % Run lsqnonlin - minimizes error function by adjusting r
    options=optimoptions(@lsqnonlin);
    options.MaxFunctionEvaluations = 1e6; % let the fit run longer
    options.MaxIterations = 1e6; % let the fit run longer
    options.FunctionTolerance = abs(max(SIGINT)-min(SIGINT))*1e-12; % make the fit more accurate
    options.OptimalityTolerance = abs(max(SIGINT)-min(SIGINT))*1e-12; % make the fit more accurate
    options.StepTolerance = abs(max(SIGINT)-min(SIGINT))*1e-12; % make the fit more precise
    options.Display = 'off'; % silence console output
    FITVAL = lsqnonlin(fun,r0,lb,ub,options); % fit coeffs
    FITVECTOR = FITVAL(1)+FITVAL(2)*exp(-FITVAL(3)*(TINT-FITVAL(4)))+FITVAL(5)*TINT+ ... % baseline
        conv(exp(-((tc-FITVAL(6))./(FITVAL(7)/(2*sqrt(log(2))))).^2), ...
        heaviside(tc).*(FITVAL(8)*exp(-(tc)./FITVAL(9))+FITVAL(10)*exp(-((tc)./(FITVAL(11))).^FITVAL(12))));

    % Trim padded vectors to original length
    TINT = TINT(1:originallength);
    SIGINT = SIGINT(1:originallength);
    FITVECTOR = FITVECTOR(1:originallength);

    RESID = SIGINT-FITVECTOR;
    RESID = spline(TINT,RESID,T);

    % Brute-forced smoother fit curve
    TFIT = linspace(TINT(1),TINT(end),length(TINT)*10-9).'; % smoother t vector
    FITCURVE = spline(TINT,FITVECTOR,TFIT);
    return
end

% Fitting peaks to Gaussians
function [GX,GY,SOLN,FWHMS] = pkfit(XVAL,YVAL,nGauss,nPoly,peakorient)
%%% This function fits multi-peak data to a configurable number of
% Gaussians (1-4) with a variable-order polynomial baseline (1-4).

    % Ensure column vectors
    if length(XVAL) == width(XVAL)
        XVAL = XVAL.';
    end
    if length(YVAL) == width(YVAL)
        YVAL = YVAL.';
    end
    % Make sure lengths match
    if length(XVAL) ~= length(YVAL)
        if length(XVAL) < length(YVAL) % YVALL is longer
            YVAL = YVAL(1:length(XVAL));
        else % T is longer
            XVAL = XVAL(1:length(YVAL));
        end
    end

    if isempty(nGauss)
        nGauss = 1; % guess single peak if not specified
    end
    if isempty(nPoly)
        nPoly = 1; % guess linear baseline if not specified
    end

    % Multi-Gaussian fit with polynomial baseline
    GX = linspace(min(XVAL),max(XVAL),length(XVAL).*10); GX = GX.';
    fn = @(x) x(1)+x(2)*XVAL+x(3).*XVAL.^2+x(4).*XVAL.^3+x(5).*XVAL.^4+... % polynomial baseline
        x(6).*exp(-x(7).*(XVAL-x(8)).^2)+... % Gauss1
        x(9).*exp(-x(10).*(XVAL-x(11)).^2)+... % Gauss2
        x(12).*exp(-x(13).*(XVAL-x(14)).^2)+... % Gauss3
        x(15).*exp(-x(16).*(XVAL-x(17)).^2)+... % Gauss4
        -1.*YVAL;
    opt=optimoptions(@lsqnonlin);
    opt.MaxFunctionEvaluations = 1e6; % let the fit run longer
    opt.MaxIterations = 1e6; % let the fit run longer
    opt.FunctionTolerance = (max(YVAL)-min(YVAL)).*1e-12; % make the fit more accurate
    opt.OptimalityTolerance = (max(YVAL)-min(YVAL)).*1e-12; % make the fit more accurate
    opt.StepTolerance = (max(YVAL)-min(YVAL)).*1e-12; % make the fit more precise
    opt.Display = 'off'; % silence console output

    if max(XVAL) < 5000 && min(XVAL) > 825 % IR sweep
        gwidth = 1.5e-2; % guess ~13 cm-1
        gmax = 1e-1;
        gmin = 1e-3;
    else % IR sweep
        gwidth = 1e-5; % guess ~500 cm-1
        gmax = Inf; gmin = -Inf;
    end

    % Figuring out initial guesses
    [YMAX,maxind] = max(YVAL);
    TEMPRES1 = YVAL-exp(-gwidth.*(XVAL-XVAL(maxind)).^2);
    [YMAX1,maxind1] = max(TEMPRES1);
    TEMPRES2 = TEMPRES1-exp(-gwidth.*(XVAL-XVAL(maxind1)).^2);
    [YMAX2,maxind2] = max(TEMPRES2);
    TEMPRES3 = TEMPRES2-exp(-gwidth.*(XVAL-XVAL(maxind2)).^2);
    [YMAX3,maxind3] = max(TEMPRES3);
    % Ylength = round(length(YVAL)./3);
    % [YMAX1,maxind1] = max(YVAL(1:Ylength));
    % [YMAX2,maxind2] = max(YVAL(Ylength:2*Ylength));
    % [YMAX3,maxind3] = max(YVAL(2*Ylength:end));

    x0 = [min(YVAL) (YVAL(end)-YVAL(1))./(XVAL(end)-XVAL(1)) 0 0 0 ...
        YMAX-min(YVAL) gwidth XVAL(maxind) ...
        YMAX1 gwidth XVAL(maxind1)...
        YMAX2 gwidth XVAL(maxind2)...
        YMAX3 gwidth XVAL(maxind3)];
    lb = [-Inf -Inf -Inf -Inf -Inf ...
        -Inf gmin min(XVAL) ...
        -Inf gmin min(XVAL) ...
        -Inf gmin min(XVAL) ...
        -Inf gmin min(XVAL)];
    ub = [Inf Inf Inf Inf Inf ...
        10*YMAX gmax max(XVAL) ...
        10*YMAX gmax max(XVAL) ...
        10*YMAX gmax max(XVAL) ...
        10*YMAX gmax max(XVAL)];

    % Set Gaussian orientation
    if strcmp(peakorient,'positive')
        lb(6:3:15) = 0; % no negative peaks
    end
    if strcmp(peakorient,'negative')
        ub(6:3:15) = 0; % no negative peaks
    end

    % Set number of Gaussians (can use 0 for polynomial fit)
    nGauss = round(nGauss);
    if nGauss < 4 && nGauss >= 0
        % 3 -> 15, 2 -> 12, 1 -> 9, 0 -> 6
        lb(((nGauss*3)+6):3:15) = 0; ub(((nGauss*3)+6):3:15) = 0;
    end

    % Set degree of polynomial (can use -1 for bare Gaussians)
    nPoly = round(nPoly);
    if nPoly < 4 && nPoly >= -1
        lb((nPoly+2):5) = 0; ub((nPoly+2):5) = 0; % 1 -> 3:5 zeroed, etc.
    end

    SOLN = lsqnonlin(fn,x0,lb,ub,opt);
    GY = SOLN(1)+SOLN(2)*GX+SOLN(3).*GX.^2+SOLN(4).*GX.^3+SOLN(5).*GX.^4+... % polynomial baseline
        SOLN(6).*exp(-SOLN(7).*(GX-SOLN(8)).^2)+... % Gauss1
        SOLN(9).*exp(-SOLN(10).*(GX-SOLN(11)).^2)+... % Gauss2
        SOLN(12).*exp(-SOLN(13).*(GX-SOLN(14)).^2)+... % Gauss3
        SOLN(15).*exp(-SOLN(16).*(GX-SOLN(17)).^2); % Gauss4

    stdev = sqrt(1./(2*SOLN(7))); FWHMS(1) = 2*sqrt(2*log(2))*stdev;
    stdev = sqrt(1./(2*SOLN(10))); FWHMS(2) = 2*sqrt(2*log(2))*stdev;
    stdev = sqrt(1./(2*SOLN(13))); FWHMS(3) = 2*sqrt(2*log(2))*stdev;
    stdev = sqrt(1./(2*SOLN(16))); FWHMS(4) = 2*sqrt(2*log(2))*stdev;
    return
end

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
function [g1mean,g2mean,PVAL,COHEND,COHENLOWER,COHENUPPER] = statcompare(group1,group2,NAME1,NAME2,TESTTYPE)
%%% This function performs statistical testing to compare two groups (unpaired).

    if contains(TESTTYPE,"median")
        TEST = "mediandiff";
        YSTRING = "Median difference";
    else
        TEST = "cohen";
        YSTRING = "Cohen's d";
    end
    if isempty(NAME1)
        NAME1 = 'Group 1';
    end
    if isempty(NAME2)
        NAME2 = 'Group 2';
    end
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
    effsize = meanEffectSize(group1,group2,Effect=TEST,...
        ConfidenceIntervalType="bootstrap", BootstrapOptions=statset...
        (UseParallel=true),NumBootstraps=3000);
    effsizearray = table2array(effsize);
    COHEND = effsizearray(1);
    COHENLOWER = effsizearray(2);
    COHENUPPER = effsizearray(3);
    ann(1) = {"p = "+string(PVAL)};
    ann(2) = {YSTRING+" = "+string(COHEND)+" ["+string(COHENLOWER)+...
        ", "+string(COHENUPPER)+"]"};

    gardnerAltmanPlot(group1,group2,Effect=TEST,...
    ConfidenceIntervalType="bootstrap", ...
    BootstrapOptions=statset(UseParallel=true),NumBootstraps=3000);
    xticklabels({NAME1,NAME2,'Effect size'})
    yyaxis left; ylabel('Lifetime (ps)'); 
    yyaxis right; ylabel(YSTRING);
    annotation('textbox',[0.4 0.7 0.2 0.2],'String',ann,...
        'FitBoxToText','on','FontSize',12,'EdgeColor',[1 1 1]);
    title('')
    ax = gca; ax.FontSize = 12;
    return
end

% Save_tif v2 (from HW)
function Save_tif(fname, imageStack, sliceLabels)
%%% This function writes a 3D array into a hyperstack .tif file.

    TIF = Tiff(fname, 'w');
    for k = 1:size(imageStack, 3)
        % Set TIFF tags for floating-point storage
        tagStruct.ImageLength = size(imageStack, 1);
        tagStruct.ImageWidth = size(imageStack, 2);
        tagStruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagStruct.BitsPerSample = 64; % 64-bit double precision
        tagStruct.SamplesPerPixel = 1;
        tagStruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagStruct.SampleFormat = Tiff.SampleFormat.IEEEFP; % Floating-point format
        tagStruct.Software = 'MATLAB';
        % Set tags and write the slice
        if ~isempty(sliceLabels)
          tagStruct.PageName = sliceLabels{k};
        end
        TIF.setTag(tagStruct);
        TIF.write(imageStack(:, :, k));
        if k < size(imageStack, 3)
           TIF.writeDirectory(); % Create a new directory for the next slice
        end
    end
    TIF.close();
    return
end

% Create a Voigt lineshape for plotting
function [VOIGTY,VOIGTX] = Voigt(HEIGHT,CENTER,GFWHM,LFWHM,X)
%%% This function creates a Voigt lineshape (x and y vectors) for an
% input height, center, Gaussian FWHM, Lorentzian FWHM, and x vector.
% X must be odd in length for VOIGTX to equal X.

    % Make sure X is a column vector
    if length(X) == width(X)
        X = X.';
    end
    % Create X vector for convolution input (center around 0)
    if mod(length(X),2) == 0 % if length(X) is even
        X(end+1) = X(end)+(X(end)-X(end-1)); % extend by 1 point
    end
    INPUTX = linspace(0.5*(min(X)-CENTER),0.5*(max(X)-CENTER),0.5*(length(X)+1)).';
    % Make Gaussian
    Gstdev = GFWHM./(2*sqrt(2*log(2)));
    Gdecay = -1./(2*Gstdev.^2);
    GAUSSIAN = exp(Gdecay.*(INPUTX).^2); % centered at 0
    % Make Lorentzian
    gamma = LFWHM./2;
    LORENTZIAN = (gamma./pi).*((INPUTX).^2+gamma.^2).^-1; % centered at 0
    % Make Voigt
    VOIGTY = conv(GAUSSIAN,LORENTZIAN);
    VOIGTY = HEIGHT.*VOIGTY./max(VOIGTY);
    VOIGTX = linspace(min(X)-CENTER,max(X)-CENTER,length(X)).';
    VOIGTX = VOIGTX+CENTER;
    return
end

% Automatic color mapping and marker selection
function [LINECOLOR,MARKER,MARKSIZE] = colormarkset(INDEX,LOOPLENGTH,COLORMAP)
%%% This function automatically chooses colors and marker sets for plotting
% data in a loop based on sets of pre-defined colors and markers.
    
    % Default colors
    colorset(1,:) = [0 0.35 0.5]; % turquoise
    colorset(2,:) = [0.8 0.2 0.2]; % red
    colorset(3,:) = [0.1 0.7 0.5]; % seafoam green
    colorset(4,:) = [0.8 0.1 0.6]; % pink
    colorset(5,:) = [1 0.5 0.2]; % orange
    colorset(6,:) = [0.2 0.8 0.8]; % light blue
    colorset(7,:) = [0.8 0.8 0.2]; % yellow
    colorset(8,:) = [0.15 0.6 1]; % moderate blue
    colorset(9,:) = [0.2 0.8 0.2]; % green
    colorset(10,:) = [0.5 0.5 1]; % lavender
    colorset(11,:) = [1 0.35 0.35]; % salmon
    colorset(12,:) = [0.8 0.2 0.8]; % magenta
    colorset(13,:) = [0.5 0.8 0.5]; % pale green
    colorset(14,:) = [0.1 0.4 0.8]; % blue
    colorset(15,:) = [0.8 0.6 0.1]; % yellow-orange
    colorset(16,:) = [0.5 0 0.5]; % violet
    colorset(17,:) = [0.6 0.8 0.1]; % yellow-green
    colorset(18,:) = [0.6 0.1 0.8]; % purple
    colorset(19,:) = [0 0 1]; % dark blue
    colorset(20,:) = [0.2 0.2 0.8]; % indigo
    colorset(21,:) = [0.3 0 0.3]; % dark purple
    colorset(22,:) = [0.5 0.2 0.5]; % misc
    colorset(23,:) = [0 0.5 0.5]; % misc
    colorset(24,:) = [0.2 0.6 0.2]; % misc
    colorset(25,:) = [0.6 0.2 0.6]; % misc
    colorset(26,:) = [0.6 0.6 0.2]; % misc
    colorset(27,:) = [0.2 0.8 0.8]; % misc
    colorset(28,:) = [0.8 0 0.8]; % misc
    colorset(29,:) = [0 160 255]./255; % ATTO740
    colorset(30,:) = [0 80 255]./255; % MARS2233
    colorset(31,:) = [80 0 255]./255; % MARS2238
    colorset(32,:) = [160 0 255]./255; % O-NHS
    colorset(33,:) = [0 128 0]./255; % ATTO740
    colorset(34,:) = [0 255 0]./255; % MARS2233
    colorset(35,:) = [0 255 160]./255; % MARS2238
    colorset(36,:) = [0 128 128]./255; % O-NHS
    colorset(37,:) = [255 80 0]./255; % ATTO740
    colorset(38,:) = [255 160 0]./255; % MARS2233
    colorset(39,:) = [255 228 0]./255; % MARS2238
    colorset(40,:) = [160 255 0]./255; % O-NHS
    colorset(41,:) = [128 0 80]./255; % ATTO740
    colorset(42,:) = [255 0 160]./255; % MARS2233
    colorset(43,:) = [255 0 80]./255; % MARS2238
    colorset(44,:) = [255 128 128]./255; % O-NHS
    colorset(45,:) = [0.15 0.6 1]; % ATTO725 nitrile
    colorset(46,:) = [0.2 0.8 0.8]; % Rh800 nitrile
    colorset(47,:) = [0.8 0.2 0.2]; % 13C15N Rh800 nitrile
    
    % Default markers
    mark = strings(2,1);
    mark(1) = "*"; mark(2) = "<"; mark(3) = "square";
    mark(4) = ">"; mark(5) = "x"; mark(6) = "diamond";
    mark(7) = "^"; mark(8) = "+"; mark(9) = "pentagram";
    mark(10) = "v"; mark(11) = "o"; mark(12) = "hexagram"; 
    mark(13) = "_"; mark(14) = "|"; mark(15) = "x";
    
    % Set marker type
    if INDEX <= length(mark) % if less than 15
        MARKER = mark(INDEX); MARKSIZE = 8; % use unique marker and markersize 8
    else
        MARKER = '.'; MARKSIZE = 20; % otherwise, uses '.', size 20
    end

    % Set color scheme
    if (isempty(COLORMAP) || strcmp(COLORMAP,'default')) && LOOPLENGTH <= height(colorset) % default colorset
        LINECOLOR = colorset(INDEX,:); % if total set can be described with specified number of colors
    else
        if (isempty(COLORMAP) || strcmp(COLORMAP,'default')) % not enough colors in default list
            COLORMAP = 'rainbow'; % default to rainbow
        end
        if strcmp(COLORMAP,'rainbow') % red-green-blue rainbow (optimized)
            redamt = exp((INDEX-LOOPLENGTH)*2/LOOPLENGTH);
            greenamt = ((-4/LOOPLENGTH^2)*(INDEX-0.5*LOOPLENGTH)^2+1)^2;
            blueamt = exp(-2*INDEX/LOOPLENGTH);
            LINECOLOR = [redamt greenamt blueamt];
        else % not rainbow
            startcolor = [1 0 0]; peakcolor = [0 1 0]; endcolor = [0 0 1]; % if unrecognized colormap, default to rainbow (unoptimized)
            if strcmp(COLORMAP,'fire') % orange-red-purple gradient
                startcolor = [255 134 29]./255; peakcolor = [200 60 60]./255; endcolor = [100 0 100]./255;
            end
            if strcmp(COLORMAP,'fireinv') % inverted orange-red-purple gradient
                startcolor = [100 0 100]./255; peakcolor = [200 60 60]./255; endcolor = [255 134 29]./255;
            end
            if strcmp(COLORMAP,'rwb1') || strcmp(COLORMAP,'unionjack')  % red-white-blue gradient 1
                startcolor = [231 48 41]./255; peakcolor = [235 242 250]./255; endcolor = [48 124 188]./255;
            end
            if strcmp(COLORMAP,'rwb2') % red-white-blue gradient 2
                startcolor = [103 0 30]./255; peakcolor = [240 240 240]./255; endcolor = [6 48 97]./255;
            end
            if strcmp(COLORMAP,'rwb1inv')  % inverted red-white-blue gradient 1
                startcolor = [48 124 188]./255; peakcolor = [235 242 250]./255; endcolor = [231 48 41]./255;
            end
            if strcmp(COLORMAP,'parula') % blue-green-yellow gradient
                startcolor = [68 1 84]./255; peakcolor = [38 138 141]./255; endcolor = [253 231 36]./255;
            end
            if strcmp(COLORMAP,'parulainv') % inverted blue-green-yellow gradient
                startcolor = [253 231 36]./255; peakcolor = [38 138 141]./255; endcolor = [68 1 84]./255;
            end
            if strcmp(COLORMAP,'parula2') % alternate yellow-green-blue gradient
                startcolor = [240 249 185]./255; peakcolor = [65 182 196]./255; endcolor = [6 29 88]./255;
            end
            if strcmp(COLORMAP,'DFG') % purple-white-green gradient (DFG contour)
                startcolor = [0.5 0 0.5]; peakcolor = [1 1 1]; endcolor = [0.2 0.8 0.2];
            end
            if strcmp(COLORMAP,'idler') % turquoise-white-red gradient (idler contour)
                startcolor = [0 0.35 0.5]; peakcolor = [1 1 1]; endcolor = [0.8 0.2 0.2];
            end
            if strcmp(COLORMAP,'probe') % blue-white-orange gradient (probe contour)
                startcolor = [0 0 1]; peakcolor = [1 1 1]; endcolor = [1 0.5 0];
            end
            if mod(LOOPLENGTH,2) == 1 % maplength is odd
                segmentlength = (LOOPLENGTH+1)./2;
                map1 = [linspace(startcolor(1),peakcolor(1),segmentlength) linspace(peakcolor(1),endcolor(1),segmentlength)]; map1(segmentlength) = [];
                map2 = [linspace(startcolor(2),peakcolor(2),segmentlength) linspace(peakcolor(2),endcolor(2),segmentlength)]; map2(segmentlength) = [];
                map3 = [linspace(startcolor(3),peakcolor(3),segmentlength) linspace(peakcolor(3),endcolor(3),segmentlength)]; map3(segmentlength) = [];
            else % maplength is even
                segmentlength = (LOOPLENGTH./2)+1;
                map1 = [linspace(startcolor(1),peakcolor(1),segmentlength) linspace(peakcolor(1),endcolor(1),segmentlength)]; map1(segmentlength) = [];
                map2 = [linspace(startcolor(2),peakcolor(2),segmentlength) linspace(peakcolor(2),endcolor(2),segmentlength)]; map2(segmentlength) = [];
                map3 = [linspace(startcolor(3),peakcolor(3),segmentlength) linspace(peakcolor(3),endcolor(3),segmentlength)]; map3(segmentlength) = [];
            end
            map = [map1.' map2.' map3.'];
            LINECOLOR = map(INDEX,:);
        end
    end
    return
end

% Lifetime error estimation
function [ERROR] = ltfiterror(LIFETIME,SNR)
%%% This function estimates lifetime fitting error (for monoexponentials)
% as a function of lifetime and SNR.
    
    % Recreate fitting accuracy contour
    ltsim = logspace(log10(0.03),log10(30),101); % simulated lifetimes
    snrsim = [1 2 3 4 5 10:2:20 25:5:500]; % simulated SNRs
    fitacc = [0.483804593905656	0	0	0.383845664248049	0	0.333244007826899	0.310578948292675	0.297259852729952	0	0.257538549862602	0	0.233785640620373	0.208919965614068	0.198342924006017	0.181546871287507	0	0	0.151809938667558	0.147336989449974	0	0.119735000963353	0.112888740498700	0	0	0.0895806454888083	0.0838245739924131	0	0	0	0.0646063027635810	0.0592683218788609	0	0	0.787874237016465	0.0456961870573680	0	0.791794202819557	0.215874235475063	0.0391074163494739	0.0330494878529820	0	0.0290125349161729	0.0264266917803289	0.459719763011264	0.680856622670870	0.0214532833217901	0.159741684380995	0.946315815131875	0.811809123601886	0.849088137615829	0.777533146437771	0.408781913626687	0.0139225417696381	0.940755841816458	0.818414386546064	0.907767080914942	0.972867274642639	0.999368327848378	0.787360575015030	0.961718187311911	0.948452748824119	0.813971788554422	0.752609553564011	0.965623664156581	0.939749102379463	0.943853289799975	0.942875465633406	0.801629351414017	0.896185270252043	0.846768807444589	0.990376133601570	0.841931991916348	0.958841909528662	0.977707535293894	0.924272829438390	0.876685851631593	0.916031101408771	0.960893597779441	0.964578587502319	0.989617874886143	0.960099346245160	0.910927733877796	0.968450403761591	0.958157497569471	0.927331324744680	0.834007006904153	0.999374068186111	0.985370687975871	0.956417725568522	0.936712508555749	0.944634981019878	0.953238962295002	0.997290625217334	0.977228951745186	0.955791238200467	0.966362777059741	0.986260174319666	0.998988600256319	0.934281487958203	0.958868481422997	0.980406379394751	0.479418919233227	0.456992354531961	0	0.393095928082889	0.0407215370713845	0	0	0	0	0.271289724923100	0	0.227688065612867	0.212849835751995	0.225220148992722	0	0.169866550371871	0	0.156756987810724	0.124877099615333	0.131866417684911	0	0	0.593923402885305	0.0966168740241538	0	0.0849084228987380	0	0	0.0144307764447851	0.0648758547062831	0	0.272757136183310	0.0550171609579878	0.909140756560937	0	0.524607108735818	0.629061975264068	0.735807539546798	0.0350133433085230	0.638714214376652	0.983796680430198	0.951687824984386	0.920047158683768	0.755076920604296	0.909493415925039	0.676980689083545	0.680189136336649	0.593264076799814	0.991450481720743	0.0170434107385328	0.902214990280624	0.973398565542857	0.997952304601189	0.716047407644647	0.929229598449449	0.943181103002606	0.853476330997019	0.840872934528153	0.929504225984439	0.885402259730013	0.960796790505905	0.925485936642382	0.879789571642725	0.930618371301443	0.841211897923514	0.936710663511060	0.979684470137656	0.994312470931473	0.996457340857636	0.998884712553845	0.983751326182135	0.942903918859370	0.985081676133304	0.936338835337006	0.954755101957253	0.998670047941392	0.916052410315706	0.982297524202897	0.990823871376003	0.976768680549128	0.999368789867406	0.995431174089029	0.973741239279880	0.983377259324432	0.964620656010953	0.965780112329196	0.992654937464456	0.987911102316113	0.989786191174257	0.964098923247668	0.974181140514607	0.998453986030628	0.985849644733883	0.998174882091025	0.998699972664124	0.997590581128239	0.989638213684343	0.994449026395530	0.991113517171518	0.988260100488151	0.985756437426795	0.491474026202997	0	0.450769873123041	0.400361150686201	0.369987248893965	0	0.318399059580232	0	0.271744859449566	0	0	0	0.232653810974085	0.208463154786867	0.192670001724040	0	0.158579273996870	0.148609101455112	0.140437038498592	0.134736808116276	0.129238400097827	0.112576157914861	0	0.910790824902753	0.462163713198732	0.0914194963626742	0.0951727762288217	0	0.0725444115257543	0.0710801996929477	0.461163627138662	0.0557800656755822	0.208327664837695	0.388719952915238	0.0494473891541083	0.946715070694136	0.220256225637135	0.837563508796182	0.927160024155850	0.711720464664979	0.452551411194609	0.884113926351603	0.936364004763225	0.837985739905324	0.0241967595647158	0.777931723411744	0.905513136451455	0.912545980193007	0.827903357334696	0.730018529422417	0.896576039599152	0.826460376996309	0.889514045084862	0.845207521566369	0.949029040196566	0.953211467866643	0.991896511631529	0.954352782595231	0.986931640390240	0.983634222008743	0.910716382279238	0.975257391710372	0.904646439902731	0.969247219709067	0.977121026596896	0.940043555557227	0.971189727649855	0.980030964081337	0.959939225156804	0.929009124953635	0.987765380729508	0.994908600562106	0.990026653397389	0.974358232182150	0.989910101021043	0.987629071309257	0.981335688037574	0.981571673731825	0.992576145268744	0.978505258908763	0.976728485160642	0.998033731215838	0.986386659220207	0.968086019632575	0.996152375097574	0.988068262522718	0.985945213241438	0.991060696906316	0.974957994498252	0.978382178480385	0.990381506793634	0.996184798153526	0.992806511799651	0.985329889238577	0.985131606284235	0.975939924719003	0.999049920644361	0.987661696573925	0.995720506777789	0.998723060285199	0.992202812245465	0	0.323794278070037	0.424627327585365	0.394142570991861	0	0	0	0.291965523945798	0	0.261514108406003	0.241732939400914	0	0	0.193107945125380	0	0.179829591256471	0	0	0.146366935436136	0	0	0	0.107326168342701	0.0992640173085619	0	0.0892805492419690	0.0784846990037247	0	0.638649242167795	0.0649974079560182	0.0627955247044049	0.181293694998341	0.0534617323534294	0.0550448100087962	0.777884195232579	0.123473193811874	0.750140062857243	0.612423352031591	0.657272414773510	0.845191594132474	0.711761883771201	0.671761820811986	0.290325457607812	0.881258530440765	0.988654306432600	0.819187627341751	0.911350193392795	0.980197249222944	0.937674444201368	0.915318711143892	0.957696833886996	0.919879962348043	0.996949202223431	0.886611977211294	0.927348270063558	0.953669761331421	0.906901762271086	0.950678695478949	0.962344332099902	0.980030928281674	0.975938872584838	0.956363522313733	0.925806251448931	0.976417407546747	0.916197737140399	0.933905100623831	0.959848170493511	0.998615866689079	0.991898280054972	0.985129841409954	0.993015500960884	0.998835573695210	0.988859367994624	0.973301903711986	0.988545328820557	0.996185925245943	0.983362184368755	0.992474803958358	0.980235722417861	0.983024230613264	0.978207910049515	0.987856996619887	0.997762357947963	0.997978770459735	0.998519331080939	0.994081359651222	0.990471353459753	0.993211522983858	0.997421620879557	0.997421652245615	0.983941702710649	0.994446470372536	0.986169904425952	0.997043976258050	0.996356898145014	0.995767965965580	0.995186484323578	0.994358319917764	0.996460420767382	0.981068764243605	0.987787214712447	0	0.512118676715350	0	0	0	0.371761895054805	0	0.321772006206410	0.768649410131455	0	0.245092458597625	0	0	0.752849487809040	0	0.183115708354983	0	0	0.143645799102745	0.128580327217924	0.120272548281823	0	0.113088564717419	0.0993498734613603	0.0992431647474308	0.0860880717664205	0.0861823755553002	0.0745433039235344	0	0.0681783342713612	0.204543469942892	0.0608374007626303	0.0530256608824013	0.909160992086525	0.297906531353296	0.726719188642898	0.612491811480897	0.0408125257261953	0.0360469422825457	0.595075798763713	0.972933251686917	0.913831766748078	0.813103886271594	0.879184805046323	0.955492121801004	0.901870905699382	0.863196245156554	0.970660342728749	0.973317520392339	0.939282601578171	0.881566700504722	0.967027348099531	0.875889520384667	0.944518239623801	0.986837656798351	0.990871675584651	0.996252887013056	0.992312326510456	0.968209060526424	0.996711549242824	0.979308254236556	0.969584036501483	0.973683963499212	0.974502793582309	0.991507963198611	0.981480420102649	0.965960994228950	0.965669347431078	0.978320652450697	0.987318892192878	0.998387398972102	0.989236653705538	0.992677618148294	0.981517738068396	0.998824391867695	0.994343317360703	0.998516344625457	0.986736039849326	0.968321425570420	0.989268902218699	0.984425995748382	0.997480245341991	0.994695068996761	0.981901854412806	0.987981291100391	0.990776400876452	0.985793322239421	0.999682441632192	0.997492513972920	0.998550930172252	0.998147087359131	0.987838659718852	0.999965914751833	0.999012526267511	0.975448874427752	0.999535592458747	0.995420431063137	0.998789494173679	0.999831450303543	0.992355737084761	0.991316686159788	0	0.438890386688960	0.413769222221821	0	0	0.355900786236510	0.345266253815035	0.326546739542648	0.294320053495180	0.263369646809253	0.819338904849779	0	0	0.512829730589965	0.478142710860782	0	0	0	0	0	0.133257151973372	0	0.480260022538965	0.951243760545679	0.0980968442007302	0.0896496751578465	0.271919543928291	0.0730628370006298	0.801810819531021	0.0708295218845721	0.114293085405104	0.932645759264442	0.838156479348658	0.559983282989959	0.726316607050377	0.806208158929391	0.820076045847076	0.788531522758336	0.888965804634614	0.863525367869279	0.988409862217030	0.839660010778615	0.827201341955919	0.839079431801120	0.965428764225049	0.983744117760627	0.906899764967317	0.918855138771645	0.950005383730044	0.988318931200984	0.966152138310081	0.968085901515996	0.955625657844840	0.989354888033897	0.999992643878315	0.997702929223649	0.996749004802709	0.983760789755935	0.950780999679005	0.994814406140909	0.982125838633077	0.993855362902739	0.980353922026477	0.971616594735875	0.969091563689424	0.977691934314230	0.985628724928129	0.998827131649826	0.992579770790195	0.992575919849433	0.974604894928926	0.996991154576333	0.995734873046936	0.990515733020639	0.985629802624649	0.993219045833993	0.998490615133748	0.995045422868998	0.994135759596510	0.994886816531936	0.998994709097521	0.997528573368556	0.987418074205846	0.994514499292007	0.996728391827360	0.998636329149723	0.995420300346707	0.990838871625546	0.994402458983942	0.990056237674765	0.996304374318564	0.996615261523937	0.998329310891816	0.996795877905619	0.998579440797599	0.999194437541842	0.993064890701737	0.998918186047467	0.996829207886058	0.998469617692460	0.999352828138965	0.528537610573521	0.470569915740068	0	0.411000629717117	0	0.360053638273717	0.334133389352411	0	0	0.266375693573791	0	0	0.230502198558747	0	0.194617093260541	0	0	0.361456785576495	0	0	0	0.120623728330396	0.311340425863032	0.176696379285106	0.854426971212756	0.494778810540203	0.551855428043464	0.0410899952988328	0.0775076646641580	0.0680135590810124	0.792714065070207	0.668383283211685	0.758752784242940	0.983178600780930	0.918474066899474	0.901636898418992	0.843090354202964	0.991079733489695	0.813940724054747	0.795256047364294	0.863975725895102	0.964210141191370	0.987544230983275	0.987251884318481	0.840074900627466	0.953090495746646	0.976175879493088	0.915297350525823	0.958751539234855	0.969124834054295	0.997090029338396	0.959273949154124	0.979749543755907	0.984707522903702	0.990960117745872	0.974276634047220	0.979774648262992	0.995034047578908	0.983651107936539	0.993831107019199	0.995434899386723	0.995975277714126	0.986654632313668	0.995337331500201	0.991081257163970	0.985032427169988	0.978504311310622	0.996037662387339	0.991250777869240	0.981969276609821	0.996977665840082	0.995749527391339	0.997985240972697	0.990364800907328	0.998771103125314	0.998072959487438	0.999183200689078	0.992231333524479	0.999139619497717	0.985961672251310	0.994061935024110	0.997796463011943	0.990181649333840	0.994101042609261	0.998551957401317	0.995745313541002	0.998263021005924	0.987986007928664	0.996708439152982	0.997884288056622	0.993287632451762	0.998851801815266	0.999900989224725	0.997669229292682	0.999141910554786	0.995721730248725	0.994969621104601	0.999139696670403	0.999853055691150	0.996591061230140	0.999305646599389	0	0	0.443276840175042	0.444001672185223	0.412426515082873	0.371197389372819	0.363598811204297	0.318511206095761	0	0.256326387944873	0.255277045347245	0.231820762101044	0	0.209099658333749	0.201433963327901	0	0.528462662052673	0	0.147584291765158	0.141744143482141	0.883447923738907	0.125638830076288	0.137799731231055	0.980277382622961	0.101197681698608	0.537852763717686	0.651104011223050	0.799260388214781	0.100404497494929	0.663964670190783	0.749584494647988	0.0606065816401997	0.555545033085016	0.689412924845404	0.901219818700315	0.908367211941609	0.744854337292168	0.933117596752976	0.996757277203216	0.826333341921513	0.877322670949743	0.994390222931019	0.973253631994014	0.945004351637536	0.987658386910235	0.882159052783567	0.958058937245081	0.998634692974542	0.982394750320584	0.986584152004135	0.997671293223388	0.984489333029461	0.979407694297106	0.984887943187103	0.980553247845607	0.999798243024618	0.991961906161781	0.994410127902199	0.994032510795986	0.973380291934263	0.996168469885509	0.997157998492819	0.991356872741463	0.996725690658512	0.991064602264542	0.990111738692575	0.994459006533450	0.999913295803139	0.986837427261382	0.995119229789121	0.998466751628845	0.999230670636787	0.993938734887230	0.989059989263884	0.996054059217016	0.998987309239351	0.987097850270155	0.995123168703027	0.999211349258565	0.994008660155556	0.996877319584794	0.999194586303539	0.998610708282913	0.993223651581290	0.997311698142577	0.999294197613703	0.999246197921600	0.996772663085137	0.995295689203181	0.999837877677816	0.994467756893766	0.997492252770499	0.999470558012676	0.999133797533862	0.996389754714806	0.997462791963882	0.998582175883753	0.997814190056281	0.996296525839375	0.999538212080916	0.997995080071554	0.551304637907631	0	0.434807408551596	0.458444197756896	0.405736644556577	0.417288551808994	0	0	0.293901665137122	0.268971970220442	0	0	0	0.222587984422348	0.183826795620626	0	0.759589890495511	0.152052600796833	0.170501566282481	0	0	0	0.117582075179869	0.668222804624721	0.506052070574255	0.0933155412126944	0.0824223716855014	0.0791177416604826	0.665340192327773	0.0720094096188471	0.325312265561470	0.896509101268069	0.827788352542602	0.934811027414662	0.876082692230714	0.934596926224982	0.951241816076973	0.878370046373128	0.987759350874360	0.985168880798452	0.969546724649012	0.956753667706404	0.956933431343069	0.986964112401363	0.940586324169431	0.958014371242313	0.960793822374147	0.961769065832132	0.989343502239275	0.998323794811049	0.982995189333416	0.965606083114607	0.972721524051958	0.983907111103768	0.981064036040205	0.990021530463819	0.975024661032180	0.994574352653236	0.999633119772141	0.997525848955397	0.989905641964685	0.998271484385836	0.994609428514051	0.985282156416879	0.992704057600756	0.995998844417924	0.993845029314021	0.987643106847795	0.997669705774597	0.997639566472926	0.999135620484608	0.999723303962448	0.993828155967504	0.996327811007256	0.999798224748091	0.997640910742039	0.995217951189165	0.998171724170761	0.991823836857447	0.995078767878201	0.993501820431882	0.997327474938696	0.999104305953197	0.999913240755734	0.994797287197228	0.995933608588195	0.998401466362230	0.999615999721510	0.995833383410295	0.999295363434954	0.999929699936200	0.999229593179102	0.999250817448565	0.993417010956373	0.999302899087962	0.997673737850687	0.998090147128780	0.999787709965484	0.997723277788170	0.997244342733722	0.998738036691053	0.476885376732394	0.494930208538301	0	0	0.395113383700112	0	0.340247953226352	0	0.311889343387970	0.266059667776421	0.257767198217458	0	0.222016910132780	0.222570784349213	0.220288611583999	0.189248658700086	0	0.710201343290548	0.594513658819917	0.463557593306953	0.134594020193008	0.116342911573788	0.740283056679958	0.332048105696521	0.885007559525010	0.969287370414502	0.0894600710238455	0.0810459171184839	0.993579834798676	0.791892821295892	0.845905052529073	0.834158334143094	0.853532511714498	0.858286772360590	0.940454404186811	0.761083409399893	0.959151750439477	0.945715604683104	0.991567751880178	0.922824812381428	0.997038199354276	0.982363882229491	0.934451557226276	0.990830076529385	0.997360528719442	0.983334330065263	0.949025368745853	0.922734939182067	0.977176972667237	0.977048958714204	0.993959063464493	0.993620227972526	0.990606300648299	0.965305525976398	0.982223822738624	0.983240858484219	0.977173483120701	0.990678805491778	0.999226147782112	0.980956383445434	0.994133255744664	0.989257445407376	0.999317367607527	0.996833509360995	0.994588420948724	0.999694820443323	0.999213400567515	0.995445126976678	0.991903794300450	0.992815962709641	0.999952209626888	0.997675630078126	0.999767220650617	0.992747041450391	0.997882832211928	0.992953016000068	0.996649875388412	0.996628008425018	0.997409172080017	0.995491410754611	0.998574316844343	0.998973825501175	0.994034699038019	0.995359038245981	0.999606211593377	0.999326793110343	0.996164637727863	0.998549454630283	0.994539443430192	0.999587258212035	0.999671950995142	0.996550812268148	0.998086943082044	0.997060205141827	0.997045641201730	0.998604395545190	0.999740737349646	0.998619310500283	0.997021405460038	0.999931401080043	0.998562281830441	0.527148400314927	0.470048251125490	0	0.425177750817676	0	0	0	0.329763628188808	0.282153855150140	0	0.252205616302082	0.244021084568096	0	0.239170363738958	0	0.187116056837951	0.180826876188316	0.168611884549966	0.158460025123815	0.142441277125011	0	0.123529254351283	0.181867379468971	0.107844796822514	0.102156613459082	0.705296830769668	0.0847901835417663	0.0842122031604453	0.757558485369804	0.715902825953765	0.531564030377580	0.941787639172776	0.869028645847927	0.815327949245480	0.804099414785190	0.896578197932178	0.992215104498690	0.882717412860151	0.979025169810673	0.966240807924583	0.973294176056235	0.930090185866877	0.924867694876119	0.986943778690644	0.954948409478194	0.985599147262702	0.994186151479885	0.995241136699797	0.983849869572398	0.975277433037973	0.979472433196688	0.983040695250177	0.963079569365620	0.979112968910886	0.998945968377402	0.996078566927636	0.998058499390931	0.987485441907760	0.999341470548359	0.990300005537099	0.979847772325529	0.999280332166177	0.993680907020726	0.997200590689896	0.994937693816011	0.992639611031062	0.995260153453388	0.992205957683011	0.999432555211269	0.996754588144721	0.992143480839728	0.991245054636910	0.995762779521338	0.991700683833174	0.998199937512951	0.994000971064206	0.999994008360116	0.992455324089814	0.998145872547735	0.997725255825491	0.998408039117171	0.998862913608298	0.997757930394808	0.997234983522229	0.999608014883959	0.998468295595962	0.998922788678952	0.999556891637127	0.997859022983989	0.997041867296067	0.997259626951775	0.999634493167819	0.998481693132660	0.998520607391156	0.999536747032514	0.997379756126273	0.997561452569855	0.999598287570149	0.999876920228302	0.999024195496597	0.997866932868722	0.513545068442741	0.477414747106794	0	0.416646500793490	0	0	0.350651843583116	0	0	0.288349632737273	0	0.262277093556172	0	0.221136045786213	0.201202710591416	0	0	0.171137478230446	0.754033072478394	0	0.133368478492590	0.0213040552132250	0.414517845509844	0.498798906039404	0.886395480034327	0.470609267432933	0.933715495926412	0.962967355972667	0.541438325162419	0.575239016915420	0.768762363087832	0.621450507995814	0.948477142551927	0.660206291523828	0.927170657489271	0.828547526315278	0.981750641241195	0.969889625696621	0.890618768478262	0.918458664884214	0.985495331946771	0.947250199912231	0.957361047364284	0.990592956682979	0.985713255668052	0.995490480829716	0.981472047307204	0.979243482187927	0.987339447083071	0.997436763230991	0.977289135646932	0.989010124317846	0.967365607749094	0.996285907407844	0.990747464375870	0.999137847683302	0.996338160446189	0.996610805617199	0.999806510369135	0.995834702463210	0.993582812325353	0.995624041633729	0.988506189555121	0.993030731426992	0.996453454717737	0.998881200391689	0.999674106914262	0.998384546713889	0.995866872283332	0.998491047715052	0.996397354609654	0.992287426897088	0.998581892196077	0.997928116664510	0.996508603458756	0.997492036385676	0.996718390473298	0.996065273832972	0.999515979537304	0.998560996864468	0.997821594822782	0.999301745186294	0.999067341912649	0.997050617818834	0.999731503341251	0.998290328988176	0.997809652590261	0.998220637000523	0.997433689772669	0.997731295497744	0.998843012020436	0.997903636014863	0.999777337584295	0.997805761634227	0.999720461648060	0.998955044923343	0.997592473880220	0.997969596673456	0.998251562964014	0.999354227385710	0.999555223726682	0	0	0.465819768942781	0	0.425292025971796	0	0	0	0	0.279924549485154	0	0	0	0.213264073652372	0.216836703900207	0.190598831436671	0.181162132811878	0.165790062089770	0.195254615061086	0.445123869761679	0.137135579994915	0.117921406644826	0.415849115844406	0.383276689509045	0.768969915770346	0.0953299585053081	0.897133837209029	0.867558842953185	0.762067198889410	0.866003030180358	0.747800684776138	0.957721539736283	0.947364037547722	0.975537298822730	0.940704181632257	0.919648193269213	0.857693033325971	0.921535044464315	0.932821934405482	0.957698472827910	0.972034828477124	0.976848321170910	0.975743643234744	0.992878625889660	0.986060829675611	0.997114245963007	0.968769439628643	0.975885949669794	0.988763284780003	0.989858289856152	0.981504812902769	0.976126457384927	0.978305923931247	0.999151843949427	0.997278981352415	0.992820814775046	0.990056556578435	0.984373683874993	0.995168118464390	0.999546805020740	0.997447413279540	0.994704042800689	0.990764664386716	0.998053567089956	0.984704218404859	0.995770520898283	0.999256821191718	0.999756455320049	0.998318208167301	0.997744860215412	0.999309248451838	0.998333515452735	0.995406725471841	0.997878193580637	0.997927779128491	0.997043826307174	0.998523381726228	0.999483936151679	0.998245666152784	0.998186626143053	0.999064532040460	0.997928076588539	0.998865483136338	0.997577675393725	0.998882721967508	0.999332429809357	0.999590230641161	0.999565624369756	0.998606266239687	0.999560981320431	0.999435663193578	0.999034345965189	0.998829550910837	0.999398127877340	0.998909970337026	0.997257963403250	0.999312549591984	0.998651030517809	0.997908225571359	0.998771631643921	0.999754593004001	0	0	0.478177365817623	0	0	0.369480475665517	0.346512208406786	0	0	0.289573093083792	0	0.259876553295433	0.0127132586360529	0.416296396903796	0.197723318441785	0	0.177013718874110	0	0.306914999731788	0.726985796706282	0.131659261689226	0.125623042720839	0.111846582731652	0.988109817090939	0.105134859215770	0.990143977848566	0.438633749773488	0.679326142901666	0.812890102424682	0.642475936443848	0.734871334514987	0.973990139038234	0.932200575841049	0.748559532890938	0.906570288649080	0.892421848737352	0.954291231454002	0.958568543760002	0.969980819842113	0.986737164126433	0.990993747596532	0.997050264677480	0.984235572874953	0.974586430199253	0.979610061213415	0.988466643110325	0.991367220063348	0.992159760526382	0.997217035783724	0.996311257160849	0.998126936704716	0.997682613581148	0.995900765509264	0.994418791679130	0.997718843966241	0.998684119185911	0.998033100832511	0.998656797508790	0.996304904203597	0.990589635935823	0.998424895591333	0.994291381260513	0.998915135078861	0.996280238958064	0.998656559264415	0.994909621284961	0.992743060528723	0.997602184277516	0.995003785116501	0.994232099644041	0.993890734077853	0.999091464853819	0.997017269406379	0.999694117539043	0.997012139943934	0.996246263893193	0.999140456586280	0.999956349936246	0.998558924371280	0.998291901303115	0.998395448952522	0.998267781747587	0.999948480906079	0.997998510762190	0.999557407285090	0.998917512028811	0.998696720779167	0.999269233693683	0.999117806305159	0.999412215064014	0.997272874614065	0.999564898006518	0.997004934451860	0.999770701582669	0.998509168349666	0.998132133410987	0.999985715504409	0.998541771570124	0.999414861418593	0.999637913285722	0.999105036212568	0	0.518177863357010	0.462009757563161	0	0	0	0.343186483897742	0.323239382903331	0	0.288593352658588	0.262072155543209	0	0.238957148647120	0.217328166334293	0.978788475675686	0.188998221152600	0.176049374626139	0.295077883235065	0.148083223558755	0.606726675433549	0.958158738466327	0.673196056508194	0.665971020913369	0.564984599805887	0.104832502100230	0.756604625177224	0.472115185329367	0.802469863145556	0.484934188137051	0.930829628923330	0.902506356994438	0.862308838422854	0.905850518410919	0.864531384812856	0.903882993765360	0.919553678932014	0.978605429025965	0.962717325561028	0.991703174807981	0.993699402427309	0.967302813827938	0.933287441833870	0.991250321444005	0.976362277842396	0.974579150662771	0.988402668565792	0.979818864789516	0.976657318514364	0.979494144213076	0.996592589858569	0.996918613612873	0.999991728161292	0.990695923800339	0.998202324527657	0.987550449928503	0.987998066855600	0.997064662707091	0.993639850238529	0.998375660739785	0.999567328287971	0.998842323986858	0.991043589621111	0.992768572596433	0.999643726919729	0.998310238825951	0.998088520721115	0.998309766736627	0.995924149704275	0.998724205517021	0.999429930784482	0.997222191182136	0.998989682908736	0.999181344965240	0.999215933920794	0.998930720898724	0.998794424500087	0.999724172166411	0.999543935941917	0.999331856493358	0.999965183027278	0.999577804998515	0.998726897249697	0.997967396937152	0.999264732394488	0.999807853847763	0.999686646619601	0.999563010351012	0.998170967242764	0.998988301437537	0.999698470624472	0.999852669235216	0.999816663128457	0.999635517590317	0.999919807154937	0.999120591598531	0.999748393239799	0.999067931211525	0.999306340905916	0.999703034203793	0.999555227013363	0.997805788043245	0.529921344783827	0.484554548672660	0.463013713775694	0.469226160200553	0.398980302210974	0	0.351555618933308	0	0.307905610230869	0.280979646796584	0	0.532735071388568	0.242098921288925	0	0.202520903859502	0.414357690745574	0	0.162927779766803	0.156262981098665	0.618263715641989	0.131087464808808	0.424607038005592	0.947023293988280	0.398529374173739	0.864405148761526	0.858730695012336	0.825327926899916	0.0825719558047009	0.831164105901263	0.986022251926880	0.997421402458311	0.974868255545997	0.862374957198132	0.984008591934453	0.980693458961600	0.876416026082271	0.942735764716441	0.928516689554720	0.964432869733207	0.964322679272054	0.941673311615853	0.998032541904616	0.988582366330039	0.974488039382788	0.976433453770883	0.999859043077512	0.988185827046120	0.994347990377426	0.989181876452872	0.988247634462576	0.986626868516770	0.992672930578609	0.997297242719930	0.996819523157912	0.992733219252625	0.997277138888064	0.994785723406567	0.999451916260885	0.999738871677076	0.998687200896898	0.998347049205676	0.996253701224179	0.997352837814388	0.995610389143435	0.997562727717042	0.996224255520409	0.996550173694346	0.997486777511028	0.997222112358019	0.999079705427224	0.998380364420649	0.999310424164097	0.999037251157544	0.999466621386343	0.995735274566832	0.997679668200560	0.998452915146959	0.998159753180667	0.998048909787925	0.999317772158605	0.999267531699102	0.999907497505376	0.998074830494781	0.999222460813818	0.997290309102174	0.999905301280430	0.998248842647199	0.999868173207822	0.998231730844526	0.999815989387440	0.998696301140988	0.997585987047524	0.999001798479298	0.999441889691718	0.999834688561059	0.999179846693989	0.998820851090110	0.998345888343366	0.999739763778451	0.999785687782820	0.999903194192908	0.524697271661206	0.526262405773560	0	0	0.401744680886217	0.375460414282156	0.364905780965165	0	0	0.820513036090232	0.277359250139488	0.244160541264652	0	0.655552769575160	0.201649905295927	0.193168377277682	0.709765495356812	0.343332819488516	0.154481570508123	0	0.551275309309048	0.456342244781295	0.500166893094201	0.615492186734346	0.643324139595595	0.688592030814057	0.656922627335646	0.492822400176368	0.848855429557276	0.982699116756034	0.947624074551424	0.960657334674775	0.935973895292025	0.951991670203166	0.990795925858043	0.894656226768469	0.983187643795322	0.980041237114720	0.970441406018103	0.966958812158765	0.979301144254467	0.996178612952186	0.966781561974979	0.984684613444289	0.977689362999406	0.971639177611821	0.999546069188549	0.990107557971058	0.990098901800539	0.995628634177837	0.997846747682900	0.985195908943167	0.992895189343512	0.994630600207877	0.993417093084362	0.998104920272716	0.994835090364213	0.996971346325539	0.998151710419849	0.998038674254678	0.997193551108701	0.997672748547779	0.997653182387544	0.997967980651044	0.999904521104056	0.995203348461138	0.995271712231858	0.999484208029665	0.999052516359670	0.998319822313682	0.997403197904428	0.999401962994911	0.998231719368647	0.999354382815722	0.999897732891253	0.999435849819602	0.999879697435202	0.999157876042002	0.999163275995137	0.998671725250663	0.998864833521707	0.998070576844269	0.998533887189119	0.998790635773601	0.999484891137770	0.998892549731596	0.999316322842188	0.998714083151456	0.999821630917832	0.999717653911040	0.998442365837750	0.999850625332151	0.999930217212168	0.999731570819789	0.999982227225862	0.999815299709811	0.999182541549706	0.999659623260325	0.999582901996144	0.999866934852000	0.999910435187514	0	0	0.477118261624651	0.419319310514279	0.389122099001151	0	0	0.369457440239969	0.101025332915757	0	0.257363646901887	0.271661109583967	0.226997692761103	0.917579968897258	0.200942837770071	0.195393668314492	0.174895186175178	0.0107083098884136	0.166277557543498	0.223907544038667	0.502537287704458	0.944463580585261	0.987734066586456	0.853610651372931	0.869365970999522	0.715461612532405	0.727107081859398	0.969556153694338	0.928501458593234	0.887534957348441	0.942375995209891	0.870642745230547	0.966206592189105	0.818704731351060	0.947058134340672	0.989038505503800	0.940478986724559	0.954021087799997	0.986892492209098	0.967692324137763	0.974036521937793	0.989846216200117	0.991184906934355	0.984165091772801	0.997124534790874	0.996255819778109	0.983622301589600	0.992618385425000	0.989013546186581	0.993729482272927	0.993086274242128	0.999306613602293	0.992104670128155	0.999198208806600	0.999626783612527	0.997433949709155	0.994493691049512	0.998237156413585	0.999676420712854	0.996474534118680	0.998755347380226	0.992709792170558	0.999609307423444	0.995777980736606	0.998310796740445	0.998992931705997	0.999871089161515	0.999257300670883	0.998822551268086	0.995483550260666	0.999463897103389	0.999152127925999	0.996569996470975	0.999295473791080	0.999823512637832	0.999881037460390	0.999089249281722	0.997940455837032	0.999055812177016	0.997327372961474	0.999362860947938	0.999197609819273	0.999401657100433	0.998999826139164	0.999967447142099	0.997757951955110	0.998818177239611	0.998868765001113	0.999679115625297	0.998982856179638	0.998814693764046	0.999435288670454	0.999878879155432	0.999928120568074	0.999910319193835	0.999494047523327	0.999234544299469	0.999927081249777	0.999759225841652	0.999054427881045	0.999712064182559	0	0.532584271556807	0.471996378483008	0.416099330702449	0	0	0.337463492247674	0.351720105767488	0.308754550497318	0	0.260930004766544	0	0.227487404916671	0.807354622208257	0.202736591890536	0.184514461004642	0.178951214848295	0.763817602770846	0.274990135357787	0.754351550658328	0.873565524321723	0.954308663035425	0.827448931156128	0.752177241461403	0.0998774065034993	0.712134505354175	0.999299316310675	0.900494112424350	0.911087808340520	0.971018493110739	0.881404165846261	0.901979749159270	0.940916192798860	0.987066283156800	0.906940291850494	0.898240571097624	0.895658129277927	0.993769227092776	0.970156433977744	0.953571611656338	0.971817693680747	0.968182991110014	0.989503736885940	0.996014835385646	0.989082909851738	0.979517428685996	0.989853137805105	0.998157629401762	0.996621444302755	0.997246172812005	0.999214423329985	0.995905465147479	0.994067646026822	0.998442079366698	0.994367324646063	0.994283250841111	0.998771857877161	0.998165642630610	0.997427464573547	0.997209039337105	0.997109238573671	0.996081226458236	0.998677535532176	0.995460879793027	0.998065472976534	0.998780466049913	0.999525181746763	0.999656028285789	0.997631851157487	0.999458858127849	0.999038980256736	0.999569059939397	0.999548614580507	0.999013204282427	0.998710055276976	0.998999834166809	0.999875434060300	0.999762732598453	0.999112458714922	0.999370756913267	0.999829614247800	0.999298289168214	0.999805100704747	0.999307962249354	0.999975779639416	0.999129610151977	0.999915045520115	0.999336590201803	0.999301691103611	0.998972854744368	0.999849646165971	0.999970359373583	0.999120518697252	0.999608607856589	0.999506721889364	0.998876970516267	0.998793021978554	0.999707527918467	0.998699596647030	0.999518797612466	0.999879550278487	0	0.543383962024723	0	0.904076845494681	0.397911468966096	0.384969826861781	0.0608552649490794	0.331257575720644	0	0.287809859964155	0.266744713512139	0.408665283301229	0	0.213287452442191	0.202441651563247	0.932968415438061	0.183888760573908	0.830664175250989	0.163639607775793	0.140899720761947	0.443658346980162	0.780770091119924	0.625539023294328	0.910419399328169	0.998990778158813	0.704569224725931	0.991075250743843	0.908249882222875	0.984994224846922	0.908762064956579	0.957495709257935	0.860866047488719	0.955516559323947	0.988333007942453	0.965468005950529	0.956770500389673	0.955909641204190	0.972805684290584	0.995885503415516	0.993898788467604	0.975140473405908	0.973912821667342	0.997531619290111	0.991288117910733	0.972570175668972	0.996431436258677	0.990543264799189	0.981718096524422	0.997275601381477	0.990840618670798	0.997400336272377	0.999992784884674	0.997198080565952	0.997233973280303	0.995730357950556	0.999722819022394	0.997625090157939	0.997610764646933	0.996712894020766	0.999345269254181	0.997762253258294	0.998743708435909	0.997298074218492	0.998824192676737	0.999413364270726	0.999996762421051	0.999469247318682	0.998073926027025	0.998791632760835	0.999247086543629	0.999208503180172	0.999824332108530	0.997971451303391	0.998941489902990	0.999154539764551	0.999565287972756	0.999436648597466	0.999444397597621	0.999613501663376	0.999675562912604	0.998650038288559	0.999773541377304	0.999714731818913	0.999906141829189	0.998473641577351	0.999858837215728	0.999439895270311	0.999522951886755	0.999450032124569	0.999067606509006	0.998402421110008	0.999544866496765	0.998782229479245	0.999666098062222	0.999794923475136	0.999724758306677	0.999414123979358	0.999618926464446	0.998898422560239	0.999250401588602	0.999737886262515	0.546546512139959	0.511708251414920	0	0	0.423890893920433	0.377471084204057	0.342211903706770	0.327010014425578	0.318779930685388	0	0	0.249545920206424	0.239829033874520	0.234639047002764	0.206481907711347	0	0.567024565300739	0.162080970343369	0.162457106541426	0.940330818044006	0.986849672763758	0.851244968295446	0.119687932670746	0.943090138252758	0.677874302424344	0.943193528419920	0.935513810481635	0.682566302018577	0.800170567320155	0.924426591898317	0.926968623423492	0.990951299529504	0.993564185201496	0.985690757666489	0.957779629899828	0.949941348466609	0.996345731346541	0.973787797890564	0.984683448805633	0.962632400809132	0.993693807787382	0.995956535732242	0.981597021511842	0.982091882661558	0.985080983475084	0.979987587994654	0.996811816421144	0.999390415249220	0.990861778437297	0.992080614002013	0.999588041236667	0.999306383492331	0.994462847589722	0.994322659163982	0.996057801760811	0.997486972697485	0.997789422525152	0.998982560554530	0.999431920114926	0.995056954999542	0.996453823159430	0.997913232441111	0.998175346547312	0.998757872948062	0.998710326969938	0.999731268552730	0.996512994861414	0.999353606741260	0.998777831274271	0.998550611844313	0.997549022902413	0.999168154812992	0.998861123552938	0.998845455961115	0.997629043878799	0.999118004725194	0.999844210750342	0.998558639731027	0.999897797065208	0.998249407692044	0.999886162547900	0.998729512761209	0.998700174890348	0.998583544211721	0.997784685519413	0.999920894028725	0.999809019775503	0.999900901364108	0.999199492628991	0.999989720587665	0.999301590955044	0.999868848432623	0.999524135905946	0.999829051940893	0.999813970815566	0.999818299601327	0.999150722737277	0.999640473505400	0.999880096751322	0.999486739796339	0.999925448706119	0.521716372592684	0.496595962038105	0.464108755623874	0.447732265639139	0.420656139197998	0.375839431703917	0.352221586791098	0.330075717110759	0.308098799319879	0.300508297558399	0	0.265058661608357	0.858446156500379	0.794006098230032	0.213421815554235	0.287606561024117	0.172950328409778	0.166675533415829	0.751541637915416	0.560776368302984	0.857802236740378	0.121963744098198	0.867523294394983	0.913794255132532	0.928542775995046	0.854137742049302	0.874787295651541	0.964695281169749	0.996725153974994	0.925291937712877	0.901143435078469	0.933332408445187	0.963968757858193	0.819388255529315	0.973582697626891	0.997059765655834	0.950737397520584	0.957016581535475	0.977421120615583	0.995114333375680	0.976774969427099	0.991218394637530	0.987386306729291	0.986856831095378	0.984812518886817	0.998283054153092	0.986092778747678	0.983264286180456	0.990878841253594	0.998950790492810	0.991649307760871	0.998751848214816	0.994648297732322	0.997618591910986	0.998221767459831	0.999378144587176	0.992497702362122	0.999875183123491	0.998797480966698	0.997061595751105	0.995753354813184	0.999559592890735	0.999311963895843	0.995546737238199	0.997860901868114	0.999413306308347	0.999886273165440	0.999497330964717	0.999521116683438	0.998217737248115	0.998431765759112	0.998459278852685	0.999862558818911	0.998510622746148	0.998043151143526	0.999199029174987	0.999099891499183	0.997858980403275	0.999574625410023	0.999541969688759	0.999750831036862	0.997517372893599	0.998816333553015	0.999872021939561	0.998770439878299	0.999176904466431	0.999632263071625	0.999649732363725	0.999142637685628	0.999481448115953	0.998774143754075	0.999069910579622	0.999837600808375	0.999944991688166	0.999694666361577	0.999572431455923	0.999775886063892	0.999923948358537	0.999571196938753	0.999918091843811	0.999661432256836	0.523714349371240	0	0	0.432090282195847	0.427095448453948	0.378289769305540	0.0285548935315604	0	0	0.283527314830860	0.272322622361213	0.270096586087310	0.230622602779271	0	0.202703202391722	0.927473610450113	0.179898204400110	0.168138216975258	0.637873651057563	0.995286377958725	0.762013450346114	0.512289833743615	0.117491812641762	0.967902780411013	0.740641667822822	0.903413350746573	0.906582006624045	0.803103566855724	0.998620386843029	0.990768967618687	0.971215682421368	0.902808886570561	0.930155628333567	0.976711125938979	0.974526939129251	0.994558632783062	0.962830126814554	0.993413573929002	0.938934900363412	0.973523187848108	0.986374183764456	0.987879360240627	0.991819763500964	0.974789421681223	0.997385294292023	0.991198299406971	0.985342187163005	0.999720250232627	0.994888343292287	0.993554831618387	0.996124549946825	0.999218440105567	0.999356476743901	0.998055647720703	0.998669225255941	0.997713991131644	0.993657477882757	0.997169556832099	0.998412998159437	0.998173403105316	0.998103160582875	0.998239886726774	0.996943482913764	0.998943022949291	0.997797352192480	0.997977261836784	0.998893637479753	0.999220629949916	0.997510384906785	0.999030576776814	0.997808387693598	0.998079493432850	0.998784103610578	0.999244281220575	0.998539526828666	0.999418094592571	0.998937189945869	0.998847900575364	0.998511052948987	0.998912904585777	0.999794090354085	0.999994620878687	0.999881920316176	0.999292089730928	0.998511192243424	0.999345391374785	0.999825772572681	0.999694268561169	0.999992797437373	0.999363783835825	0.999317516953678	0.999979848997553	0.999658783913201	0.999905126797869	0.999099520579382	0.999780376568372	0.999966643487470	0.999938472074426	0.999974272620325	0.999943627191739	0.999316600231347	0.531293549634593	0.487165909573038	0.529429289695668	0.470801655126800	0.478926971239831	0.395425866719336	0.354806194310231	0.353086910114082	0.303144828174926	0.302370304824181	0.275000067234826	0.986209817781003	0.236764605530915	0.216584659104640	0.300846051146539	0.188668651114023	0.695342174293631	0.633629793525098	0.163773616263422	0.147792226416366	0.773446333208225	0.555294961747423	0.953454352337105	0.851421859094848	0.801070518637849	0.433063942600624	0.946157683740372	0.945017349735072	0.778542112781396	0.920402718468260	0.991246908272017	0.969549427945297	0.929825669569151	0.948492491827812	0.996509916845771	0.968784855734546	0.968117408656700	0.972840873192459	0.970064729714058	0.963039836964540	0.987318706282848	0.993079870447325	0.993473089156945	0.993343085563771	0.990586155670891	0.994636533224811	0.999547494976783	0.997994861563426	0.993017062667184	0.999612134232680	0.995118561558745	0.994036991280087	0.994774954518010	0.998344122269561	0.998463934702327	0.996499907462193	0.999992122209550	0.997398505649070	0.999962115359303	0.997593522449181	0.996570507397216	0.998091531532155	0.999969223329760	0.999762610792325	0.997592041484870	0.997814455823159	0.999142985439038	0.996842940643952	0.999492024818310	0.999129830747476	0.999849904087039	0.999678427163909	0.999513869515120	0.998263289203592	0.998941307707276	0.998892218305999	0.999379405779326	0.999319718289374	0.999387502466996	0.999691581176925	0.999939184039046	0.998559698198723	0.999973399918044	0.999996147079823	0.999826734567788	0.999802571840673	0.999529393104507	0.999787705003872	0.999313603014604	0.999786297125421	0.999751894596535	0.999705321613059	0.999631431808622	0.999325958165955	0.999851026631626	0.999925271937830	0.999976955748870	0.999702115077677	0.999272254371207	0.999718125665123	0.999831078065734	0	0	0.485566908810306	0.465648134306411	0	0	0.373245805113794	0.326196752099845	0.967719854114726	0.293063216450944	0.126263335378152	0.252012202515072	0	0	0.989966484813464	0.345743334445934	0	0.760560840376017	0.152171531767616	0.502622043528292	0.662423570351533	0.121279242452486	0.671477155327866	0.847414851350272	0.725877610516361	0.876758451192613	0.687225827345310	0.910680080737853	0.885883227678670	0.899808538957052	0.890924318671103	0.975173873029239	0.983929828323414	0.995640602289032	0.979868846035031	0.974886226596991	0.999906210566412	0.988630499236267	0.997271933552632	0.995206886726634	0.979039224024660	0.988123696051824	0.993042083539128	0.993276142107773	0.991127644640976	0.996607405803399	0.987681431214241	0.990672247386578	0.996467051221785	0.997373666340672	0.999115475963668	0.995468358194270	0.999869125203938	0.998936397291795	0.995575759091271	0.996846082673810	0.995455883845876	0.997937761999587	0.999863335858546	0.999634513509362	0.999500395737913	0.997939837614932	0.998292320727154	0.998720548046315	0.998089198314787	0.999285789957293	0.999356970032057	0.998663275168066	0.999826628791682	0.999433756949099	0.998206466023248	0.999087042045116	0.999567174617637	0.998985900847899	0.999042981989587	0.999091508076467	0.999910364299689	0.999842617628110	0.999616004564166	0.999644101890329	0.998915365948428	0.999160081107774	0.999574999660815	0.999646205788688	0.999708250838760	0.999967793862952	0.999959903413906	0.999892165544223	0.999058367697512	0.999838738670195	0.999167420841848	0.999311292428054	0.999761598129824	0.999879411448944	0.999718776229844	0.999679301972552	0.999853570791093	0.999947224363938	0.999911221283941	0.999855473877654	0.999782751477611	0.524065211063766	0.303939475360127	0.473232590592308	0	0	0	0.368491241864407	0.329884338610291	0.305203175628594	0.301452706317348	0	0.266008289691019	0.270189714016235	0.229241499978257	0.201536526927221	0.416119760732492	0.258060360029406	0.723542854264876	0.967967434552612	0.970323951856200	0.922170296949754	0.870173985449941	0.876740249504956	0.932624225287941	0.708231811580192	0.959252062129151	0.985784242370124	0.912242893787048	0.924484776772659	0.957859883150392	0.918083366841774	0.946101014781849	0.928172435368247	0.995206225567658	0.956805489361040	0.946491835171898	0.989200911683814	0.990894471366607	0.974910834047545	0.992228469371787	0.979136412688385	0.980073392029319	0.998737322935415	0.996957028032992	0.997556676454770	0.997065751789258	0.993174553370889	0.998194315725798	0.996621325991683	0.995768564439731	0.999507593260293	0.998626810117778	0.994958418897938	0.999242485661167	0.999563564463106	0.998279598566908	0.996037151606097	0.998716931495750	0.996980597843702	0.997806719975323	0.999194116061489	0.996570145175163	0.999869304240607	0.997145223202139	0.999259326271402	0.999187318653230	0.999986197719159	0.999910760197669	0.998344064193485	0.999953961248677	0.998980962555719	0.999029617296489	0.999206054775149	0.999162650964425	0.999493181949501	0.999655840181326	0.999224243439462	0.999812483054357	0.999557098493484	0.999958197578251	0.999988517877181	0.999490922831904	0.999960705618675	0.999155251945393	0.999430559330081	0.999464088242771	0.999446771002876	0.999523120960965	0.999481295194893	0.999702342093155	0.999933435549702	0.999864601197136	0.999954137363464	0.999855531900082	0.998779753753740	0.999885837509356	0.999861068113642	0.999910748138749	0.999851445435691	0.999593785245189	0.999281352733855	0.553907381219184	0.510423387466307	0	0	0.402245813750298	0	0	0	0.311651122262329	0.289409510389865	0.0965589223772938	0.268543772469202	0.0585720417751423	0	0.210498490902030	0.191114711628317	0.180779824813570	0.165620587429043	0.156761678117453	0.942528976452633	0.317991136180068	0.577367366724403	0.931673573215570	0.950140820732621	0.752235926035535	0.989608367866636	0.901724608918317	0.980093016228696	0.893901092541908	0.864013196929403	0.999378517083722	0.890840828927700	0.930227460627881	0.999754371897482	0.963807900291623	0.994087197573307	0.992449123079780	0.992937967300113	0.961596561838173	0.999635420088547	0.990868631448303	0.989694726590706	0.991821561939547	0.983129830349442	0.996488426833886	0.992475358888874	0.992182561955006	0.999780715417206	0.998047435519080	0.999489135085043	0.996317630844971	0.999343803458695	0.998086053701591	0.996688789995536	0.996884374152334	0.996964055374399	0.998585706782415	0.999211220366675	0.999946192344197	0.997111678329901	0.999228251002268	0.996657286137756	0.998964459240589	0.998879577366747	0.999932021825192	0.999706530906453	0.998649282771595	0.999386587387940	0.998854618431599	0.999957769600098	0.999159509615901	0.999618979848121	0.999370069282132	0.999441291848589	0.999746812074896	0.999844436170097	0.999454821128798	0.999351671952139	0.999782062454427	0.999639173088296	0.998791348246242	0.999214686117935	0.999684252499722	0.999496008911615	0.999216571170539	0.999907227835518	0.999249047280185	0.999860084457028	0.999806921729813	0.999392914488656	0.999522534617483	0.999403751840206	0.999463795260288	0.999351080895070	0.999683514407145	0.999833080597546	0.999935256863114	0.999684605012187	0.999587038942740	0.999985228961783	0.999434115926859	0.641259096116458	0	0	0	0.396336171792982	0.0890788921363171	0.346817421940249	0	0	0.616335333115588	0.269215016488220	0.301184313963078	0.251799777365205	0.763240974765302	0.216601477951751	0.00668074963238996	0	0.824849805772324	0.978021754206644	0.957868001847289	0.949160104029173	0.124850439481062	0.942926811101821	0.770726688477757	0.728640817615106	0.985549695278144	0.857154412606136	0.927143826114379	0.937285530446532	0.983396926372169	0.990714325684800	0.967846815320327	0.897423528420133	0.954953143453925	0.982451610127076	0.983476108143199	0.989756004021319	0.998581870722979	0.994479511060820	0.985224972178391	0.976134389865864	0.997053888101584	0.994066431726622	0.992015638362588	0.992383732419833	0.998173432734113	0.995663770644854	0.989016335025087	0.989516038774491	0.995476030041247	0.996115989297925	0.996503022744668	0.999196143599760	0.998000411273351	0.997931890070268	0.995629652441526	0.997805859356458	0.996229181254981	0.998672632637355	0.999987318715631	0.999494400270736	0.998776975828988	0.998646092849034	0.999035591305373	0.999471639248670	0.998383866476762	0.999395496047175	0.998178535766332	0.999729830024696	0.998420933424638	0.998038151632654	0.999499276048326	0.999943652576426	0.999614531922079	0.999706084505288	0.999305081717365	0.998420211761422	0.998949303361438	0.999891671113951	0.999956441799422	0.999367248368011	0.999813941119235	0.999758304038642	0.999592865618735	0.999126327160904	0.999255020367820	0.999768219961171	0.998724023941950	0.999742127624833	0.999496837247955	0.999536096865787	0.999909323272616	0.999891775421752	0.999920287397212	0.999951360286523	0.999292427712909	0.999972323810503	0.999619790231589	0.999795132088774	0.999255352958243	0.999750790922258	0	0.492632199075193	0	0	0	0.376385221331683	0	0	0	0	0	0.0364630725079821	0.250604969362020	0.288641895721663	0.484230958407535	0.190117682009560	0.185014234360920	0.711176924163851	0.781989830805800	0.578220012607810	0.892711406738964	0.421808511470182	0.535808780176545	0.829656074210865	0.975804224511452	0.827513528980249	0.911502295627982	0.833471724972531	0.859386845041739	0.939581602605962	0.968964078213973	0.984354483291005	0.996948043608622	0.987645015156114	0.977171394240643	0.971086380573372	0.990393709155313	0.961202001869240	0.992161760408205	0.996840138703693	0.977744440944898	0.989270912094994	0.992697998364309	0.995796974922396	0.995503898092799	0.997186502400374	0.995226134825086	0.997110962983916	0.990892836627235	0.998249283576783	0.997399543807217	0.996060428397889	0.994114592894654	0.999948343483697	0.997505602974852	0.995780853024379	0.998825677850874	0.996852192565358	0.997398222153759	0.999450735029707	0.999915981175254	0.998740537225063	0.997180093230482	0.999307150432750	0.999258947794397	0.997853389917350	0.999227392522627	0.999279249843153	0.999865137496611	0.998381665733167	0.998859254879265	0.999516878145589	0.999653991601160	0.999970637915502	0.999748350486753	0.999506080095664	0.999931527317282	0.999987284944234	0.999531128041756	0.999663015795475	0.999108148631310	0.999609987010585	0.998854796315944	0.999801031569340	0.999888469879704	0.999789693639810	0.998719002605748	0.999732930795390	0.999637019173056	0.999529237307254	0.999184994772787	0.999659465074508	0.999856285855896	0.999948532386581	0.999736610583154	0.999878706556073	0.999989788279723	0.999943686677570	0.999974339153714	0.999944587365878	0.999650137655861	0	0	0	0.436639234751020	0.417978804590886	0.168746519377048	0	0.352536162224083	0	0.294575773428311	0.846684388439119	0.257100130547655	0.251530020993324	0	0.218278240984720	0.758135858958969	0.969870105669200	0.177516474097796	0.355897485370975	0.906895973118916	0.707736418827464	0.727060494789808	0.861267326315255	0.979174838416547	0.904915211108796	0.880864593038212	0.972976500276099	0.989361646893972	0.996000138602600	0.974064119307303	0.975254743650591	0.960942339484935	0.970895124000193	0.996876692379616	0.976143111655807	0.993564543551637	0.992007220016879	0.986008287486752	0.995506018445321	0.992474970620283	0.990173880040526	0.999494992138353	0.990616537516762	0.996918534008781	0.995988861250260	0.993789174141751	0.996353448454242	0.997020394306380	0.990922053858849	0.993229983729925	0.997900051104780	0.999586566725877	0.998706574999201	0.999465840999007	0.995338793305637	0.999266623706542	0.998814665748892	0.999664645303040	0.999361469306107	0.998973108420124	0.999127154519865	0.997621193823656	0.998173994280157	0.998913680507454	0.998726831173551	0.999188630390482	0.999032667378923	0.999445026376658	0.998579444914268	0.999920568914620	0.998150171967422	0.999822860800294	0.999033859131470	0.999841285851859	0.999695464856143	0.999831074545267	0.999358964202501	0.999208932384958	0.999941735238825	0.999652969903924	0.998994354275545	0.998892931493311	0.999903483944756	0.998751821288108	0.999994575848613	0.999620350276065	0.999741416619830	0.999865125806450	0.999041576360917	0.999995037123431	0.999822589889786	0.999890383450553	0.999335115276327	0.999544453776855	0.999727627868151	0.999842675262301	0.999930657178809	0.999861253518655	0.999713822290339	0.999851837359921	0.999500421561027	0.575972111827561	0.548748631418245	0	0.570097243302292	0.408320618127939	0.387009997798428	0.0189677786143099	0.328096036566629	0.308669151715736	0.285140348824253	0.273572999452097	0.341476100595821	0	0.217400910520028	0.849793913782314	0.448397088660422	0.181050582694570	0.169593278490395	0.154266811171456	0.533165539621687	0.494810477112703	0.808367038787400	0.818380910520039	0.737417630734690	0.994539224202948	0.917448200451719	0.902093726394266	0.970772689639140	0.986040024019646	0.922526658104022	0.900980004649862	0.999553245324590	0.941535686119643	0.974561858198225	0.981367262503064	0.993725996265188	0.991232525007371	0.991546489986651	0.983512032870597	0.990772116994236	0.992484175711869	0.991644632126367	0.991007779271104	0.984129458145096	0.997494428313352	0.993234139579730	0.996378079516465	0.997168312554497	0.996429219014161	0.998644910286132	0.995503107376502	0.999471214227808	0.999624527519596	0.999439447561967	0.998965234079856	0.999364376333550	0.996937924377240	0.999041607393198	0.999289504785515	0.998337275779774	0.999113757658221	0.998446864792289	0.998762465788484	0.999369740504745	0.996726335737412	0.999038329348416	0.998688146989606	0.999696574975614	0.999535729753105	0.999708616935722	0.999468995034728	0.998226651802031	0.999334620919419	0.999466253807753	0.998558823193373	0.999868451070427	0.999477899206733	0.999631317101378	0.999882112614351	0.999780185631472	0.998732481260751	0.999619101988858	0.999599373174641	0.999689874974210	0.999311965538372	0.999979741462777	0.999403048815794	0.999918893063073	0.999682819922022	0.999978272178084	0.999741896364916	0.999916744136384	0.999845650058678	0.999338315262919	0.999764204363373	0.999988540034575	0.999882067063227	0.999625816584126	0.999669878198761	0.999644484427605	0.999473210444547	0.531161577108217	0.520943637459422	0	0	0	0.399462759699842	0	0.378133048747411	0.308317037079352	0.287528551458876	0.267964386396618	0.292889376379157	0.388973320101809	0.809332794930342	0.354322402675927	0.931651449039130	0.744154715711415	0.802350040131034	0.166508843161006	0.182564422046207	0.966812759788123	0.128176257369595	0.774592712735356	0.566761688766145	0.756056916015667	0.801831373707611	0.935202528799681	0.996679507786962	0.999123339940115	0.973861539539894	0.961328123368490	0.986085371435265	0.970654780080363	0.980121560773439	0.996238868638124	0.990851314604058	0.991716781266527	0.996278218360426	0.989649518770873	0.996918517498604	0.988887941293918	0.986161281936907	0.995885248751467	0.991187217116941	0.999529073526566	0.991138760785317	0.999692955657190	0.997311132136234	0.996038504428710	0.998278972812202	0.998737503720956	0.995299749003137	0.999717302844212	0.999103770304030	0.998907637189407	0.998884125042535	0.997966056952485	0.999891554307008	0.999799364798751	0.999583157910517	0.997929939078971	0.999730640701207	0.998779645798277	0.999740416249947	0.998437725692921	0.999332038730251	0.998520155568357	0.999719084351638	0.999661729597814	0.999695359334714	0.998578285266962	0.999690135494426	0.998468294884984	0.999660727344609	0.999655814366613	0.998260557717868	0.999177104720072	0.999450906198786	0.999513109422848	0.999390732238696	0.999818002648838	0.999783037339540	0.999866820121595	0.999698135264231	0.999811900220477	0.999650809597076	0.999677720806141	0.999578568496158	0.999771537233929	0.999640007574147	0.999751511436644	0.999966269195553	0.999463668048704	0.999773339455160	0.999544743561976	0.999993797315481	0.999660467003644	0.999932809514976	0.999845576093061	0.999496124899994	0.999879255584026	0	0	0.489790598036245	0	0.404692751258075	0	0	0.418330754107265	0	0	0.756979730175089	0	0.252418307411873	0.655940945754751	0.757524595345463	0.199236312065652	0.434080397825294	0.993995153811613	0.411125850764533	0.821139711010917	0.668010494160262	0.932463461773138	0.805538304981147	0.893773805561486	0.816074540679019	0.888416992886679	0.778903658243116	0.926829238744091	0.923551269202753	0.955700950657212	0.969566793136470	0.965891882941178	0.953126641450843	0.964273822621088	0.975914948911695	0.974043228659529	0.983685883386028	0.985881219312669	0.976691696037000	0.985168778520505	0.988189150043138	0.994363164358259	0.996650126694513	0.988308251813764	0.999285993094479	0.998065365001148	0.992730615001552	0.999200861372807	0.998861295494854	0.997866676158746	0.997680827367972	0.996385167998843	0.998681174986814	0.997170669303925	0.999551017270188	0.996908333053107	0.999374268770324	0.999945188720120	0.999817554155667	0.997935859541495	0.998799857964253	0.998959069257815	0.999228139073837	0.999901775468630	0.998793987153194	0.999868691404716	0.999249913716287	0.998249025148210	0.999808597079267	0.999793382933222	0.999781387275543	0.999351672519101	0.999087992566986	0.999679240993292	0.999567806966269	0.999776713378673	0.999690968111785	0.999856344687306	0.999721441961982	0.999815609096704	0.999997581232181	0.999595653731481	0.999834689496006	0.999872981775708	0.999349050172947	0.999850402364208	0.999935475658423	0.999396775982977	0.999696077420664	0.999267891203960	0.999270488177674	0.999805740626477	0.999922672035319	0.999172453879596	0.999864030404278	0.999943841740496	0.999754760612649	0.999525557171423	0.999598523776145	0.999984413429837	0.999961486561304	0.953074842898814	0	0.125193953619012	0.452466383723016	0	0	0.503038019353263	0.332028397698871	0.327249851819329	0.465386776696336	0.266442228615318	0.256175491953410	0.466808630211513	0.687791057955475	0.465685876072207	0.0979043322841376	0.180791548203625	0.295216316086113	0.156009230331253	0.799665539596172	0.693954714870502	0.975259236393564	0.985476896409255	0.760660422150151	0.852676523405524	0.866285974453001	0.977000719687552	0.971407249955523	0.780232619693438	0.945543645713734	0.996151136241472	0.986693737753960	0.968882042061795	0.984631378855354	0.933101291549032	0.969126868298929	0.993334459023954	0.992195874397526	0.989544468994658	0.996505954731254	0.999996759475400	0.993914923276632	0.995447495672782	0.997199124996120	0.994863453610100	0.994959531011419	0.998266161774448	0.995705705816665	0.999556234488116	0.996798987361884	0.995919595987362	0.997693432955560	0.999910457506895	0.998309982420013	0.999566501328117	0.997330393310753	0.999922659010362	0.999740074989376	0.998135039329383	0.998399871855310	0.998378702239592	0.999134691506420	0.999233338185506	0.999429323769180	0.999470930158390	0.999455726545914	0.999740722182121	0.999700699873295	0.999304930407327	0.999340664351312	0.999854167374053	0.999464627713913	0.999763590605300	0.999006565288191	0.998573603084985	0.999941838082549	0.999429624007228	0.999084640086140	0.999470391947777	0.998398463401168	0.999892848370021	0.999946387508272	0.999753038413965	0.999874956768600	0.999753068918323	0.999731382837561	0.999741430477665	0.999864161305953	0.999681387699725	0.999788678553118	0.999499589631415	0.999665876096958	0.999702214593310	0.999937695647154	0.999920887594874	0.999668193817397	0.999801868442599	0.999527403239690	0.999892969460350	0.999969237196008	0.999554520865829	0.525865013170822	0	0	0.516033787776704	0.0547440666419569	0.417880068984234	0.353393177911025	0	0.758987900935569	0.628655732434171	0.804688376615942	0	0.633160730100616	0.860651770408415	0.824403729741166	0.709377879098656	0.591132000978000	0.395543929460580	0.737780356544465	0.669257219505185	0.135179519955335	0.948650419677634	0.890684090122716	0.931948360116169	0.886330924417163	0.758132850727821	0.965950837171004	0.920495511752368	0.980715328344324	0.993413752339946	0.889782597309610	0.920204768374702	0.966282031407518	0.980377111687052	0.990551962700782	0.995379747641097	0.976779162526344	0.999734225999964	0.986978587834986	0.997271587375698	0.996741288247222	0.997819906754115	0.998840290147066	0.998059667994299	0.998532708411379	0.996506081580295	0.997762841929242	0.995487957024881	0.995111352173173	0.997114310209194	0.999361745345937	0.992597179854096	0.995241397423835	0.999276584940491	0.998100012733191	0.999538698282540	0.998477285563704	0.999034064608266	0.998199173248270	0.997633008599460	0.999000278611972	0.997526047733658	0.999552959630260	0.999310868577484	0.999362536603402	0.999461520961088	0.999619836798987	0.999209161904495	0.999915562578440	0.999417689149298	0.999287901245485	0.999317958115618	0.999523581162441	0.999809272540547	0.999571241336555	0.999543835527218	0.999717014633357	0.999321402666170	0.999986616251644	0.999848367018906	0.999322648074381	0.999888138440999	0.999589862866424	0.999568620267882	0.999731541153807	0.999852270271164	0.999992265391704	0.999282592643856	0.999603394225986	0.999457260489645	0.999584313601844	0.999936786005413	0.999850992994255	0.999884551599911	0.999784834842844	0.999601133520240	0.999920876833203	0.999736577526930	0.999685314102602	0.999884215591628	0.999889277840281	0.565766699688002	0	0.504694737021879	0	0	0.843488674400304	0.144349565896774	0.0206621716468401	0.447748875113517	0	0.412928212516417	0.258653103293146	0.238599489560177	0.132939773124137	0.591357049853193	0.860183879052875	0.900820701937747	0.169811466709041	0.596457252243281	0.152532223060358	0.586823583393621	0.864329507930341	0.684356671752067	0.759498539665300	0.910053285855936	0.773981751044591	0.812197936761930	0.941755354544708	0.917397939093150	0.972768494854754	0.968993631034086	0.960561501383746	0.980342877202015	0.946373458833713	0.986822564122655	0.995251000904218	0.990719067529640	0.991279796257290	0.989514697397682	0.986086395658557	0.991602419790459	0.997460034234978	0.996140319241497	0.984771624899661	0.988391087358387	0.994582733405817	0.996515003878986	0.998056203737560	0.999511121641501	0.999887357129917	0.998579075553286	0.999881526049078	0.999312174351577	0.997094633020145	0.997915363733394	0.999980064999394	0.996014041945151	0.997643209811173	0.999360003421908	0.998544288974708	0.998354772869977	0.998775722096705	0.998453169390691	0.998993688370027	0.999836561924967	0.998840730908299	0.999430301431382	0.999707274318685	0.999219758996876	0.999518504632651	0.999091083286376	0.999988659545861	0.999294691312589	0.999156127943421	0.999937616284593	0.999734962715921	0.999186178716887	0.999991339821661	0.999633008396397	0.999523605447027	0.999870614372039	0.999841421812729	0.999716233188650	0.999959100434493	0.999578298498413	0.999748634768441	0.999385251051324	0.999771467952839	0.999642637354070	0.999889017678049	0.999914969517332	0.999915255850878	0.999721180284869	0.999978067573004	0.999874015822052	0.999930382205441	0.999920246768929	0.999934094388968	0.999988431644448	0.999756230353913	0.999899614304503	0.926411491239806	0.521447049868577	0	0.440917778977924	0.422648393033654	0.798483797323893	0.356256901915707	0.0750860562911266	0.793710935660989	0.290759956682626	0	0.254404887772884	0.232213883849515	0.757397105441716	0.776265637913238	0.933491553165532	0.517317835100472	0.780170380687572	0.892851707545435	0.560239846395402	0.140841283162986	0.932917390236991	0.801089278128298	0.918010877850757	0.769963683727814	0.996688653322389	0.992114062036525	0.923577520099020	0.978494922998127	0.938802234765020	0.986166279179659	0.996523031850302	0.997562475223353	0.977738212926816	0.988123124794210	0.965312984379963	0.990617471267650	0.998310250490266	0.993259011464834	0.993976602655351	0.999650279640958	0.996073495959120	0.993485272402775	0.998794952118547	0.994480228326025	0.995109601320695	0.997368462099620	0.999882030773662	0.994275864039331	0.998140867303391	0.998859800392293	0.999180190140966	0.999060841083681	0.994471316937888	0.997406726417624	0.999092332777052	0.999256660316013	0.999877424253722	0.999529104621284	0.999816142534100	0.998383872417055	0.998863564216652	0.999127315250887	0.999011646377873	0.999922371540793	0.998202681266112	0.999475601526146	0.999863398724760	0.999583361272864	0.998800818458544	0.998742293208832	0.999126153694431	0.999308069953520	0.999297065906223	0.999948085819014	0.999080788948772	0.999880719602728	0.999826004748839	0.999684800251877	0.999719603323347	0.999874917133601	0.999791806849945	0.999947796326396	0.999529174133607	0.999995060889122	0.999948412924862	0.999797969636953	0.999725016790134	0.999644078421543	0.999493752996473	0.999950647114028	0.999511104545196	0.999790783008759	0.999612533965658	0.999771872462512	0.999863424943705	0.999937724503490	0.999728477833386	0.999854319049389	0.999936769998631	0.999717350731823	0.544761176636301	0.565038821640313	0.466785852950041	0.470239689986042	0.425672055293028	0.287732087863292	0	0.990681124717498	0	0.279820107474022	0.275576874449297	0.264679640570871	0.237341620532593	0.292454149926570	0.202071002435282	0.566790210217359	0.536944802168177	0.845774813867638	0.756187047309043	0.810962390578528	0.901727704571418	0.899444126051767	0.879245346981439	0.963413373152443	0.997868383851497	0.872585273698297	0.929007611361688	0.945216607518266	0.945337400248664	0.981714942825308	0.958419925940133	0.975098046655370	0.975859864448268	0.968899459649138	0.993746699829687	0.995501620655392	0.995962433807957	0.999980345883767	0.988131806367476	0.999426908064409	0.990078266339998	0.990254233641725	0.994935107170543	0.994098612308892	0.995896660900380	0.998416918872137	0.999153956819624	0.996623799718235	0.999316773698089	0.999490048772032	0.996867853574192	0.996748007879604	0.997070784672964	0.998162740723511	0.998626226093370	0.997734592695517	0.998247997992998	0.998039096450702	0.998460134796124	0.999065201314521	0.999921171859455	0.998521456556952	0.998971655355045	0.999275984834936	0.999169619134785	0.999919823252792	0.998677626726604	0.999699375232367	0.999317585555684	0.999570575781189	0.999757786626321	0.999755269986567	0.999674290315373	0.999458374133035	0.999981874404000	0.999883338396202	0.999743915377745	0.999423547416159	0.999483699138419	0.999767172306946	0.999781821747598	0.999985874988367	0.999638093625359	0.999824580570151	0.999995831775340	0.999767918110482	0.999955155917198	0.999803851594866	0.999634566079914	0.999924687865741	0.999772437649282	0.999998671770900	0.999742201123036	0.999969179288724	0.999777580782443	0.999899575839358	0.999983844590712	0.999735604479776	0.999730809639887	0.999963780735636	0.999948698350116	0	0.510539106235301	0.493438188905705	0.891932299313859	0	0	0.366531732753881	0.613098952600916	0.785482464588635	0.306069476292733	0	0	0.702611780532345	0.237942745989933	0.843288547213303	0.204932326124496	0.908536102962397	0.754952803714362	0.807644423064316	0.892685963964846	0.965404879204751	0.967645789052921	0.682248684329594	0.997450060050935	0.972840960585020	0.931603576382291	0.977535384108423	0.912709779850627	0.989919990528778	0.985513761117257	0.964658235898100	0.989631118058055	0.971186673685502	0.955368750461407	0.964152785801696	0.987031889441442	0.974120974301383	0.991910757926288	0.974117086416081	0.995622700626920	0.999171989562373	0.992379513659977	0.999940928063788	0.995727351017804	0.995326154926812	0.998504325382064	0.999840713412577	0.997347987888091	0.997660697494319	0.995869757058766	0.998358340701845	0.996281750305422	0.996509657514502	0.999051810603279	0.999222982563175	0.998752538347109	0.999503963131198	0.999568708028222	0.998858538756986	0.999503543977942	0.998691416672594	0.999409361797952	0.998748743430831	0.999037397830721	0.999240239762443	0.999364188469713	0.999463765141635	0.999283289309239	0.999914567753146	0.999713552402132	0.999552182034539	0.999841956700001	0.999859374563442	0.999979072518230	0.999693222995526	0.999994738255016	0.999659791558794	0.999640660919579	0.999764937366177	0.999629292860124	0.999874698288067	0.999732633713644	0.999603287507171	0.999909562659495	0.999861562412035	0.999689168082274	0.999975545710913	0.999923434546516	0.999811639304023	0.999735476884991	0.999983629266511	0.999369815871791	0.999883997243352	0.999901484668373	0.999820625676288	0.999919907491763	0.999757400186478	0.999978065175961	0.999858204497680	0.999898015349797	0.999971942480733	0.183707693328377	0.524001894033505	0.488617510165318	0	0.402840844772891	0.385385175919422	0.361496706670402	0.329982826498861	0	0.305360671295386	0.172064087975021	0.114806513042085	0.346198429521672	0.217022728726647	0.872836024980936	0.480637199661142	0.834910225020893	0.607876656788060	0.831360645377269	0.932918568809338	0.797885763111908	0.898893846511338	0.823698459135763	0.898652886196903	0.896235184962873	0.958143860641605	0.979438789645498	0.987731742305551	0.998988568157834	0.960638953030275	0.923591180812634	0.952341660128238	0.993857618400819	0.977765590426056	0.993703416363577	0.987670549300949	0.971738734096551	0.982616429542017	0.999094600187485	0.986746532890344	0.990904847586705	0.999365421548621	0.994328915065636	0.999964254356593	0.992265272581393	0.997757761693193	0.999317484530230	0.998765093818506	0.998318260644786	0.996947618909806	0.997639689312142	0.998966180978745	0.997428622041613	0.998894867581353	0.999049146792054	0.999587981273849	0.998821695045118	0.999176221768525	0.999298758047637	0.998710053900523	0.998663509815814	0.998363926254931	0.999097131016280	0.999257611109377	0.999671953639017	0.999655608166627	0.999069901725494	0.999760510411346	0.999835666301031	0.999602367797495	0.999845886066292	0.999661748727712	0.999735894160520	0.999257985951936	0.999753571723822	0.999896504206094	0.999495360408491	0.999989306986857	0.999805243666753	0.999512395291677	0.999636787911140	0.999994322192524	0.999787056325992	0.999615347048908	0.999466587597033	0.999481849220651	0.999657738739798	0.999962281068365	0.999999638855984	0.999713827364317	0.999755758674894	0.999979087064605	0.999684450710158	0.999980793636318	0.999617588746285	0.999800442809691	0.999990056836329	0.999691592743267	0.999925170246719	0.999961244506291	0.999720817788388	0.824892288075957	0	0.496117201401302	0	0.133983749302233	0.668693857473662	0	0.349864388602699	0.308740083934049	0.306365846066428	0.278305985121772	0.759360683101166	0.829045430485782	0.493773293840665	0.983983033925138	0.271043511919926	0.862245230174963	0.193277326606384	0.816868476108178	0.953907685483110	0.655407460873846	0.874034906457011	0.881880217294960	0.867141174158856	0.947028328325262	0.990096364335499	0.908179159503827	0.958640441794845	0.997939404112146	0.978455930585644	0.978324068417572	0.989820581234782	0.964087812999347	0.953143633318366	0.988143061320793	0.993094739119297	0.993787229742050	0.996281771746410	0.996495357786267	0.994554022198045	0.998329468358308	0.994711510230670	0.996525418068506	0.992672499230249	0.999581267238883	0.999395894350767	0.999561366699684	0.998660778546106	0.998495741620049	0.999183246101696	0.996593367463731	0.999860820555711	0.999564120976872	0.999208167607492	0.997147582658295	0.998082563900557	0.998176602457088	0.998589000104116	0.999856653747468	0.998322688827471	0.999689031440031	0.999245513591606	0.999185106558737	0.998232900307256	0.999745580001551	0.999401656963818	0.998258615085063	0.999149445583687	0.999423778432582	0.999919616850389	0.999960998743067	0.999662241409939	0.999948004218011	0.999643366936142	0.998771283395716	0.999668105467893	0.999955509070260	0.999679651285171	0.999560697114913	0.999791331546144	0.999635694692847	0.999886594178804	0.999995829950394	0.999961198598998	0.999666764615528	0.999594524428520	0.999970004905799	0.999739421979574	0.999588169261525	0.999484449038749	0.999994709216521	0.999724166336685	0.999779021580043	0.999836453610889	0.999898640154013	0.999660830091528	0.999818433802574	0.999815212748010	0.999812995198029	0.999900637052908	0.999967159550405	0	0.498192902140033	0.460768806627449	0.432151785352130	0.122476530504628	0.394083778339102	0.370260208168153	0	0.330690414433641	0.713350869013238	0.267806462654394	0.529464072535158	0.795820112921944	0.699614602950140	0.766622552658018	0.191250771577215	0.683027411918039	0.636930012148067	0.926534724079530	0.750906388920449	0.267873919831134	0.811300088174518	0.931652217885261	0.884831337093096	0.904288924069927	0.895493932556883	0.865978806658533	0.916879754863486	0.903149414122301	0.918464136699063	0.934227456972477	0.996352437755373	0.980588855469249	0.988202695011127	0.990332416488430	0.997244727378473	0.996598292478102	0.996190100859626	0.996987712693876	0.999247208531731	0.999064622967224	0.993679650435524	0.995235120314831	0.996407770189941	0.994840931946953	0.997779866827292	0.994873424621112	0.993225628579552	0.999894992879525	0.995500648049986	0.997104942146508	0.999457202193385	0.998313210985012	0.998460802554073	0.995542934899848	0.998930262850213	0.999095012831562	0.999842467597935	0.999313274782346	0.999106735888522	0.999378594040762	0.999565789720101	0.999318050143053	0.999515127122022	0.999844847451837	0.999622444163065	0.998778211718894	0.999717025099000	0.999205702275972	0.999484241024952	0.999100011179228	0.999215810535168	0.999933893502909	0.999949260736696	0.999948706929994	0.999904170955053	0.999471616063303	0.999661141527770	0.999850876425661	0.999692762173266	0.999696084055063	0.999963212728129	0.999442901758962	0.999990867194324	0.999871239488081	0.999992514882063	0.999906682934793	0.999921730260915	0.999853349658550	0.999779871273989	0.999807644720272	0.999730045344378	0.999700310349424	0.999710324395764	0.999661621575569	0.999928156127235	0.999709222585595	0.999714838614439	0.999865456133700	0.999934082208377	0.999690508922377	0	0	0.289102052208040	0	0	0.0499988809312928	0.342842099737717	0.353576221921666	0.493003057473538	0.305481430770811	0.589618784855726	0.244323069642036	0.719478771675753	0.826593632213942	0.275054803838622	0.892541409183751	0.186837220162713	0.166697820123617	0.643196898439505	0.149389493018290	0.828598456209506	0.123942485686641	0.871925310228478	0.867500604270056	0.943013858187329	0.941262698584994	0.970926142641664	0.983781821508997	0.973548573984852	0.953697776923178	0.929055163007602	0.966048698336538	0.969310297540589	0.974426610308716	0.935850825110980	0.996432036750603	0.993674012059787	0.992343731444382	0.998460726935257	0.989166888529303	0.997428638659018	0.998401739734338	0.999145703239164	0.997463957048826	0.995422021869504	0.994873907692810	0.999209622204250	0.996637282537389	0.997020327938236	0.999474248301929	0.997986640186910	0.999997868617563	0.999456573871030	0.997797781315164	0.999949987436126	0.999214413257813	0.998507296598752	0.998849771875872	0.998556576634418	0.999146219219819	0.998751412111533	0.999245425718982	0.999934896592588	0.999849355944970	0.999981023578609	0.999576519232756	0.999788624486921	0.999988124453157	0.999464036416991	0.999517143277033	0.999711276061124	0.999887703340598	0.998902589788730	0.999615768172343	0.999474616843464	0.999747710851069	0.999723593126122	0.999879665842958	0.999825268029804	0.999973765145987	0.999687922535300	0.999711012043577	0.999547119808144	0.999472204653984	0.999968252923431	0.999954519672891	0.999558228897692	0.999898387536740	0.999737619635602	0.999962670838499	0.999869383981529	0.999716276024006	0.999987876622711	0.999835396571873	0.999815494420112	0.999744532990584	0.999962321059818	0.999881503964015	0.999889994036735	0.999959405501446	0.999910829241243	0.559600592223275	0.491626470139772	0.467797518333573	0.466214258600193	0.403371954500457	0.393201011986981	0	0.354770615477217	0.320401218462630	0.287651890210212	0.852564684074502	0.253239836038409	0.254406584425189	0.233384337014081	0.779624031329647	0.832326649405393	0.415734992122594	0.634021539135151	0.950414644485221	0.150740429596176	0.674967544297504	0.958063988205299	0.922212273363112	0.784959015070599	0.815084819120321	0.960383265932942	0.999742291544030	0.958712955218460	0.961344480445970	0.975794339880236	0.966915327297390	0.973359371185597	0.997119749861930	0.970283381993529	0.999731614914657	0.991554912640119	0.981095779432131	0.993869403209596	0.996566104204999	0.999393394677398	0.990045660391651	0.985019592152399	0.985576112167125	0.995466013684497	0.998893942930892	0.995811954161507	0.997204669012266	0.997574497973561	0.995968420651876	0.999955105491856	0.997850165930418	0.999280738091025	0.997076675189080	0.999890380520231	0.998903170688161	0.999873252680318	0.998896519216582	0.999686255857302	0.999848752244520	0.999539948812651	0.998764263130931	0.999356143925175	0.999449027416118	0.999355082244855	0.999553030745740	0.999656296517149	0.999179785335511	0.999901424511273	0.998844379791612	0.999862266036407	0.999527115088629	0.999771496395142	0.999103048836900	0.999037692950385	0.999634972741200	0.999464779889076	0.999079627667133	0.999360379370558	0.999837916036891	0.999807686863679	0.999766813271515	0.999858582871681	0.999558615892792	0.999670359454519	0.999822249783107	0.999485292982504	0.999681243964070	0.999501432846294	0.999773420983669	0.999970096673063	0.999769383667903	0.999976575773321	0.999727716973947	0.999856623055584	0.999978683559816	0.999989852192268	0.999919773190896	0.999717208838258	0.999828637286341	0.999766625625509	0.999874152952900	0.533624494006699	0.489867529710857	0.474650602051169	0.718082418744218	0.434851538425713	0.385872740261520	0.364700419976072	0.714044311433992	0.613551856478989	0.352416673831650	0.580932967414616	0.811704809553654	0.983468659792252	0.216603097569370	0.213776423426207	0.765276015196823	0.640833679972463	0.177403480316630	0.548223579596655	0.574229754519654	0.577838200635914	0.872450994292409	0.962092089886543	0.754446954603487	0.893369658511839	0.965631298256672	0.923516780561982	0.997376806643597	0.989579666444778	0.926693979173917	0.921901714948390	0.965382059026637	0.970640588924566	0.985379406105625	0.988539143833902	0.993983105475389	0.999964806281215	0.986781912545364	0.993311182997669	0.991956455385932	0.995948340969754	0.995485795745899	0.997175690967831	0.998435326550023	0.999690004184599	0.993639157134689	0.996813981240571	0.997699360963979	0.997286859028092	0.997601903175888	0.997223179875367	0.999832812241658	0.999302138726093	0.998671250838770	0.998821008538757	0.998311897249366	0.999195305501262	0.999763105725685	0.999836062343979	0.999342774681759	0.999226524941289	0.999656044124156	0.999542592000582	0.999221700566693	0.998857207055308	0.999441630160346	0.999560344892332	0.999858114140724	0.999728718115538	0.999538008696753	0.999562867939852	0.999692507773843	0.999852973661656	0.999684115668580	0.999525806791756	0.999585543365746	0.999108660322541	0.999387358240804	0.999012579912651	0.999269979225770	0.999730183168140	0.999996994055313	0.999575272093216	0.999417262021804	0.999457056533320	0.999378381680839	0.999738896883603	0.999877992522104	0.999944790928675	0.999994517784992	0.999709623861112	0.999787271669205	0.999802484368833	0.999559511787040	0.999985900098179	0.999853014250764	0.999859842171895	0.999880162636247	0.999782445548569	0.999895657403105	0.999993428364612	0	0	0.459099523831853	0.462033407442761	0.404724846724969	0.210099111994456	0	0.364511467471513	0.939561209476079	0.287716120603665	0.269427090323627	0.758209049754565	0.231254416069128	0.707433338207845	0.621741153523205	0.950739972647000	0.741811665723489	0.737070861848803	0.951093454818737	0.991209716574242	0.821082800376578	0.906115108952765	0.823599786946470	0.899886470135704	0.901151528858444	0.919834977113577	0.877924038199379	0.955957178581795	0.991842897593654	0.993431979642764	0.989988833086309	0.962910810598845	0.975209863766692	0.991486476425724	0.989711875664144	0.996999092002936	0.969742808029364	0.988672764395195	0.989972341056407	0.997755980068454	0.995052076056983	0.997381137540416	0.991493134068198	0.997674773633004	0.999342535292454	0.998491476228904	0.998675878414278	0.998958730061482	0.998336379853665	0.997206250106190	0.997607374580762	0.997598974052197	0.999136819333235	0.998587523103109	0.998373081929650	0.997597195233050	0.998427022543358	0.998637370999667	0.999220480936429	0.998574917775426	0.999809204436857	0.999883008169015	0.999495108153986	0.999878560465908	0.999954088021971	0.998899884354513	0.999539213325725	0.999646940703734	0.999536893590234	0.999937983681458	0.999872404248196	0.998963709995251	0.999689013108922	0.999453074794548	0.999949609385105	0.999767998178212	0.999928338301275	0.999976791200254	0.999393800040615	0.999777097094658	0.999832249973475	0.999898538268784	0.999432272312363	0.999772116531618	0.999902055433232	0.999897037993957	0.999865308815931	0.999957188100285	0.999469632316730	0.999990374980478	0.999892863051021	0.999781034218355	0.999707744009418	0.999975551558631	0.999871178990830	0.999799387753438	0.999983797028493	0.999967784357260	0.999583896827242	0.999909254029840	0.999899067240092	0	0.496727609109391	0.544735654592465	0.463031506428134	0	0.481897186187395	0	0.325347298296326	0.779379034666070	0.580738539091887	0.295515266219159	0.267243586381559	0.248986405067161	0.234872442277215	0.948917116780414	0.436367524708524	0.299246459647514	0.836151359906774	0.516792685663608	0.919379467198407	0.894482233894901	0.907083150302086	0.913601137067875	0.973061696015948	0.964874225413765	0.949881984454073	0.981259939814495	0.977721697991572	0.982642757350839	0.953546525074272	0.959486583415606	0.964904137546583	0.974398276163681	0.990341921683791	0.982003713497049	0.972792329234103	0.990483360081201	0.995296600417304	0.998520031913482	0.992396735279767	0.993631254781078	0.998289063097107	0.999498586308473	0.997843505185046	0.997120363825620	0.999554257136705	0.996466231871212	0.996028872721657	0.996921364436232	0.999142016255810	0.997440349036146	0.999104948924472	0.999354745594372	0.999708027366019	0.998761182932093	0.999280211719006	0.999047528155659	0.998643450770562	0.998547162084616	0.997681376500876	0.999538640714572	0.999884750506045	0.999466846342297	0.999800090866516	0.999789067505933	0.999859727483002	0.999938220417445	0.999467334579077	0.999857189548369	0.999444783714947	0.999729555380424	0.999788336149270	0.999492066735033	0.999844023522178	0.999993325933725	0.999535407256379	0.999772685053933	0.999671379497121	0.999995743290267	0.999483028546118	0.999471594384348	0.999494909266037	0.999988414861196	0.999848306380280	0.999741250473086	0.999971576257893	0.999747391435004	0.999886287433476	0.999894480475771	0.999544429072666	0.999935300548359	0.999975167356029	0.999800290852210	0.999964586521341	0.999699047714244	0.999609681179245	0.999955747802446	0.999893076006452	0.999972086473006	0.999860086658893	0.999916930080999	0	0	0.461192323593968	0.453489597478649	0.432199653203229	0	0.371812776001324	0.0785186598046004	0.234766870392466	0.309911967047419	0.391455609070208	0.813639374015217	0.249797765319663	0.130852419175473	0.201702459640508	0.196288407186544	0.723197037726205	0.958612156656229	0.650908492096948	0.929294529648184	0.903448076023208	0.839943006675688	0.974329881309495	0.875929672710720	0.993614028131898	0.879655352950608	0.964876924585776	0.974483172976691	0.946349592007100	0.962666614721027	0.956448447300434	0.991612954980612	0.978917645652253	0.998323088040948	0.988396103412790	0.990760084554858	0.972418329696263	0.999278918319075	0.996904428007202	0.999408675109397	0.999709719442590	0.997533746318206	0.996187970751792	0.998299233531069	0.999397201506996	0.997591472877565	0.999364139895220	0.995894651483199	0.997363453423843	0.999916451810433	0.997470422530613	0.998969453831704	0.997603657348907	0.998362090614106	0.999512294825015	0.998489708279279	0.998537067230764	0.999835487508536	0.998950477839912	0.999113565785494	0.999736827326781	0.999517978853880	0.999981085051108	0.999095033192328	0.998944653202918	0.999800027572368	0.998492500966030	0.999943640168276	0.999993403022115	0.999358604490260	0.999800642859538	0.999546855816011	0.999277502816450	0.999953920489706	0.999502330587484	0.999481405220968	0.999003712178587	0.999838471083797	0.999850130680436	0.999629808850032	0.999759920971540	0.999607784434559	0.999779723484887	0.999955492740917	0.999792377994734	0.999901260872748	0.999851065993449	0.999749107261489	0.999935428751678	0.999728270035910	0.999968209255257	0.999903757918358	0.999760404318200	0.999603217525474	0.999648807847806	0.999960368337220	0.999936795689721	0.999969917746712	0.999917212395313	0.999968914475298	0.999987040165001	0.923208050704300	0.360802113949494	0.962886994900872	0	0	0.875430118900273	0.368561662357670	0.352554447984423	0.00890056424217522	0	0.276879629614272	0.585008181816723	0.373663766785873	0.322098209868759	0.986804525759496	0.578989697237755	0.828740888792687	0.453963606859757	0.744988896587465	0.807036918765596	0.948799396703630	0.993809841198551	0.883769839330390	0.815094459406884	0.947492847899347	0.994195734059796	0.939255225611440	0.969369861518949	0.988128110380800	0.985441133900593	0.983948364955120	0.997116337598865	0.997837933450460	0.983419252352812	0.983631044338394	0.963688662034070	0.999826226255982	0.995889569824782	0.999072119535936	0.998595720271388	0.990301197888476	0.992633131378024	0.994348979305052	0.996743101458922	0.996843694863475	0.998486216823969	0.999804815081289	0.999991711953363	0.999575558127718	0.999695940384411	0.997124163190680	0.999054393643317	0.999742465827857	0.998130452444826	0.999144475187012	0.998881973936458	0.998323237673430	0.998841132294924	0.999098647822877	0.999927332873287	0.999831957711278	0.999366775644043	0.999169925265328	0.999345661352889	0.999947567698838	0.999528957077952	0.999300547990190	0.998961458002451	0.999734984039329	0.999669536438051	0.999668976806470	0.999902899361565	0.999657220468795	0.999559018187187	0.999555310375211	0.999966533521124	0.999766532884669	0.999411343792761	0.999555414431108	0.999475248807999	0.999638208593802	0.999767026037293	0.999446370263699	0.999876917578392	0.999974176532612	0.999952253404753	0.999951069871821	0.999923065647547	0.999921559634095	0.999993410620204	0.999933169326683	0.999899758196904	0.999843932334058	0.999825684176921	0.999902348957548	0.999813696403507	0.999917717493125	0.999809526348809	0.999890124912385	0.999934652493201	0.999891352011365	0	0.536668634392717	0.396557589000062	0.450243589835759	0.519970120929300	0	0.379743837820516	0.356646426183170	0.332459965923795	0.425226801834222	0.599352715793706	0.684464940330049	0.0918416278960045	0.654223533034507	0.205766012498878	0.704113838042983	0.744912090049396	0.594291033566181	0.878832955249593	0.950717988873373	0.612565008329168	0.931775860883486	0.965357443063193	0.892024985813675	0.962410240097464	0.982732184095319	0.990297210878175	0.926947674235588	0.956674214914993	0.983543984227162	0.992699265158120	0.923662196295612	0.992404005378301	0.995205441177656	0.979902978341551	0.983696687388316	0.986866200305814	0.995818306062614	0.998236878032016	0.998240360109282	0.998293436889279	0.999799769663867	0.996612627237756	0.997453750952257	0.994902113999870	0.997754917394923	0.998802571805034	0.998616298572558	0.999258151200034	0.998284232550243	0.999476221933328	0.996688790821553	0.998952272552957	0.998484031623209	0.998648253219136	0.999216313626092	0.999996051504356	0.999177584827838	0.999046729922100	0.999794261073367	0.999850271423365	0.999139231535851	0.999490863691880	0.999669235546445	0.999940848356813	0.999658548925408	0.999824592271771	0.998720137090189	0.999679095681107	0.999822834952044	0.999894628821044	0.999551740379229	0.999918006863139	0.999900466046764	0.999914600923398	0.999357644855750	0.999950323666487	0.999782229438453	0.999919962465228	0.999880504189633	0.999668397362373	0.999876996395996	0.999895496463546	0.999989807608335	0.999565221025482	0.999845842568681	0.999876274279623	0.999906697949363	0.999958997135611	0.999552554601755	0.999386727195239	0.999957562102167	0.999959878186984	0.999853184266356	0.999915968385692	0.999872785794462	0.999885640465404	0.999745251587530	0.999910491110935	0.999938768989251	0.999985521027876	0	0	0.507913743144600	0	0.709631184552752	0.604908004182683	0	0.528300870750969	0.0874269237322266	0.309398457508619	0.268658516565908	0.883189287222472	0.634013950498779	0.219274687950594	0.962114491262409	0.202836512719837	0.850207304447309	0.171298112710500	0.742770909003032	0.832986702487924	0.998254176162889	0.893874620316126	0.966360283134700	0.953956659807979	0.973212454359677	0.988626883955197	0.914636264942122	0.931279325979373	0.990631660614210	0.935031922622135	0.991402353611296	0.993293618395188	0.984076142728687	0.968744826057110	0.985993607268401	0.985977395912899	0.988566957856207	0.991801060316530	0.997619509341272	0.995913677619471	0.993484797276589	0.995574863514104	0.992594857663478	0.997941182687069	0.999138462091139	0.994296662706452	0.997264953882002	0.996819107859266	0.996914152861225	0.996221533746350	0.999169335546799	0.999842088122710	0.999273614736793	0.997900151570105	0.998419710927569	0.997272559704439	0.998669347889553	0.999312025085311	0.999494106522570	0.999647768419736	0.999455959333210	0.999357342289211	0.999622094730585	0.999961369157079	0.999952946264360	0.999057225498938	0.999885160483702	0.999273569796730	0.999591580573279	0.999937338926543	0.999784532966306	0.999718979630521	0.999655484607177	0.999662893192706	0.999518981138630	0.999983120237033	0.999553136439811	0.999780295565563	0.999633578536451	0.999970450588264	0.999660902245949	0.999990635082141	0.999789542133819	0.999756419113881	0.999450703247864	0.999885396591152	0.999969004989801	0.999873243035600	0.999737240102563	0.999857283211560	0.999795643461105	0.999872936129068	0.999916070292207	0.999893405389301	0.999909610599778	0.999835234990779	0.999983200348132	0.999652896656473	0.999875703466124	0.999585538952646	0.999906560243711	0.503135665456350	0.0700621085846840	0.466907332574756	0	0.402345600941542	0.395242600005157	0.369977905662818	0	0.149075989316048	0.121164189126986	0.914267220252135	0.478193385405960	0.250966546631734	0.254472997288480	0.669107362106643	0.443263234048434	0.543439583962528	0.995273473290893	0.961800834990469	0.804752040088401	0.866369118162654	0.998218211767571	0.830278791661544	0.948448589421389	0.905161531724147	0.970128549049374	0.902904314973587	0.976881932610023	0.983680692290582	0.968945865665572	0.964333197423532	0.971270754349358	0.956977611340777	0.993966511938493	0.996155005006315	0.996674794097476	0.989046164012622	0.987257595543195	0.986183066036024	0.999755462272772	0.997290830214699	0.999482463276039	0.994967573936098	0.990876530156786	0.999678509679638	0.998088862613762	0.995782234562622	0.995857938976119	0.999462944091203	0.998420302225459	0.998841061201289	0.999905607419707	0.997733022248934	0.999266482544951	0.999404605171719	0.998753009288558	0.997708664528186	0.999491853244772	0.999975417155170	0.999868497070427	0.998950868850471	0.999769288961458	0.999286276151221	0.998376093507026	0.999978865691337	0.999801700669472	0.999775083109447	0.998962000383038	0.999719717725897	0.999775021489043	0.999607908831821	0.999483353317415	0.999845740788460	0.999845983969926	0.999606553299439	0.999931640691915	0.999798429497108	0.999828615842665	0.999994010060752	0.999719896996960	0.999483539073763	0.999793322692035	0.999866173103220	0.999658552371536	0.999926789801545	0.999837178340513	0.999839670939431	0.999847982124677	0.999361086831786	0.999853318595214	0.999907721359261	0.999894993856798	0.999820930137966	0.999750535524850	0.999997695078985	0.999765798707846	0.999899636148944	0.999864599326593	0.999969476232410	0.999885174690208	0.999916731023657	0.557841190805778	0	0.497719266813089	0.433960524088008	0.419927731063103	0.970513349535979	0.354253008840170	0.559430570168368	0.976139662952910	0	0.549174481832400	0.786857507513785	0.878503359141188	0.876778923033119	0.910769735587707	0.842699708771818	0.866060917724243	0.935462282451346	0.659669806807974	0.967742808628994	0.969319317643171	0.849425914321284	0.765118423560246	0.979974692129238	0.985570348732968	0.945166733468917	0.937858737305247	0.972560818750802	0.987182366150862	0.958964466989611	0.952411713998150	0.957638690701262	0.974915966064450	0.994001014001741	0.985204530760363	0.997943010197004	0.997556425105943	0.989139068582247	0.988474379398142	0.996214029844017	0.997418606385489	0.997486942916271	0.997147242277576	0.997727899342158	0.997795781010074	0.999672435480189	0.995024065808048	0.999788774228906	0.996978369968209	0.998998976749541	0.999659856252740	0.998877310041885	0.999666023764189	0.999361651219709	0.998822789821646	0.999386050044876	0.998575614245396	0.999944781539227	0.999203949986583	0.999367609685007	0.999699740078670	0.999331761880172	0.999366489829573	0.999714638680002	0.999433278992632	0.999897239406682	0.999553902157657	0.999596364322371	0.999686619448778	0.999762272870268	0.999681154351267	0.999889196158029	0.999796684042712	0.999549064391113	0.999560099214214	0.999753853584220	0.999933537455690	0.999652565406202	0.999729842905770	0.999527996886261	0.999428729807839	0.999572730628117	0.999892376352250	0.999994275321416	0.999771866267868	0.999997687040130	0.999849027579440	0.999869260878860	0.999910841042122	0.999871627059006	0.999940093870794	0.999985549613645	0.999844246661025	0.999802412897551	0.999750703963863	0.999766390846196	0.999868755045359	0.999673632484958	0.999919064366901	0.999995596756108	0.999924612021049	0.559464862781518	0	0.482705661329982	0.00472294491434289	0.434620491218856	0.355297879149165	0.228176409386655	0.340039464047970	0.783102422253471	0.805057127164771	0.386360096789277	0.679209128366422	0.965700296453952	0.379624820688418	0.736957127411631	0.887341168177125	0.820154759944715	0.491675512773731	0.957611695403711	0.996679517771258	0.383812175120013	0.819324137403349	0.793769610199079	0.925636831069204	0.913550301564768	0.984521804809730	0.996585845844783	0.975820969666727	0.982736337078123	0.965248128950287	0.983533478537884	0.984963675207860	0.997182871020666	0.992428538598518	0.997522531453683	0.993173677370826	0.999584180489744	0.996770099356638	0.983533572837390	0.992566871768525	0.999441689955477	0.997821544277427	0.997315141153498	0.995612866701865	0.997715514776659	0.999258595795007	0.999254217230205	0.996995451233102	0.998611860570592	0.999574267465860	0.998406697689170	0.998486882430986	0.998883279233168	0.999498230548367	0.998795973708645	0.999532988202107	0.998289998437124	0.999319578032736	0.999603682716863	0.999694918291401	0.999360590339817	0.999054820788083	0.998491184213188	0.999547480496279	0.999094767856865	0.999720378470302	0.999876113443870	0.999945539023616	0.999623711505349	0.999947139230791	0.999814240653454	0.999799573120302	0.999993684935597	0.999651450534334	0.999415247792342	0.999859415681866	0.999574901803506	0.999769046804235	0.999844971099222	0.999809384800925	0.999997576145223	0.999848131403506	0.999952106433377	0.999789926704508	0.999895382743466	0.999622547672570	0.999763863605115	0.999992671897095	0.999970798682564	0.999837262239434	0.999746645064127	0.999863020848575	0.999991693942671	0.999885009813201	0.999912025769760	0.999737096848137	0.999937946361077	0.999997474797881	0.999906915049600	0.999955924738442	0.999805133055916	0.556890479574035	0.499207831384560	0.486629838149055	0.622121002608539	0.477695985607262	0.518092097028079	0.353205251107110	0.169614819185591	0.465982328885837	0.434533542413259	0.266444957235070	0.249402407539249	0.664832161009396	0.759675330796971	0.213976511983973	0.888136595791366	0.180842856219790	0.971356264476561	0.681074696034422	0.867335355772632	0.943575044912022	0.727233391225893	0.808567117830577	0.873184389069433	0.984918674410362	0.976605368398913	0.994875284716862	0.996891960297077	0.987530249820256	0.989644929760854	0.985068307395202	0.997732552698904	0.967355971143105	0.997802036430356	0.994209663346604	0.999112281247115	0.988992693779153	0.998187838982004	0.995463695671328	0.989772599275437	0.991352901925685	0.996581951833213	0.993981344429025	0.998984661705964	0.998607853383752	0.996811111595510	0.996855410777445	0.998003258270583	0.999548277890545	0.998917736530958	0.999090772629657	0.999142037511362	0.999764731624963	0.999443481245420	0.999375715633800	0.999626154080402	0.999756605425167	0.998567655483371	0.999261914565348	0.999805460466672	0.999857968032266	0.999519685586421	0.999406520103193	0.999425182880294	0.999621573236509	0.998528676272842	0.999886358482801	0.999602381565613	0.999182134479625	0.999835055548927	0.999978318427330	0.999802663787271	0.999955559738096	0.999818136466669	0.999607549460801	0.999815639872856	0.999910240319352	0.999676056776584	0.999950009364796	0.999582212009346	0.999746607823349	0.999543605922856	0.999954867284532	0.999998042749288	0.999988821743712	0.999757926386129	0.999719825655496	0.999709081193441	0.999890091241928	0.999805677464637	0.999927889888452	0.999484053900688	0.999992732777337	0.999767972245779	0.999832273636790	0.999895377922327	0.999786047661303	0.999831998299950	0.999822896986840	0.999918030848765	0.999958592832890	0.531242280551160	0.504121969946036	0.567767470023920	0	0.398669461206442	0	0.994053512204489	0.356182153699707	0.313857390972694	0.175297607432709	0.278678922029326	0.275343731678515	0.239689921038541	0.609482576050355	0.911425134556636	0.783567241677443	0.973072741743272	0.871121886080385	0.756465989236754	0.917158097048886	0.615398814198294	0.958384209034539	0.863335710051490	0.953492689733497	0.943338186786994	0.950742756961611	0.990872030909317	0.990302653081977	0.995797744587239	0.996383726774611	0.995651542777260	0.994683562478660	0.981955388935651	0.982659482885142	0.995582676109457	0.997723995262779	0.998832172398229	0.995291481185968	0.994824281775442	0.993748301911564	0.999289332225586	0.994611540784913	0.992745731718877	0.994480832678600	0.998584578241229	0.997366591942098	0.998242085005780	0.995816003857229	0.997947276314646	0.998524019947988	0.999869448372178	0.998700624474873	0.999966056322618	0.998957161227716	0.999766251197543	0.998530625848788	0.999852822465464	0.999552464278104	0.999900343574098	0.999856623989220	0.999749810509046	0.999517520721317	0.999924365331160	0.999547649630963	0.999834610553642	0.999542892538597	0.999622910688793	0.999911665891668	0.999833764039799	0.999703684849091	0.999494006812341	0.999809838089918	0.999920607029760	0.999881447865256	0.999162363927167	0.999718275111446	0.999661349382000	0.999987380326596	0.999845749041802	0.999846288635749	0.999610607921066	0.999740157303804	0.999753445710247	0.999946386251867	0.999932089506227	0.999735993374266	0.999868222214804	0.999717464370424	0.999652487083664	0.999859664709172	0.999862989286036	0.999768188264346	0.999718695864946	0.999868828461788	0.999989024555227	0.999846903323301	0.999762404444564	0.999630979658245	0.999723846765230	0.999813043276606	0.999742977272088	0.537057283091749	0.284666585341796	0.698873968169064	0	0.519523652823560	0.396813382128130	0.353011106642376	0.345251542197184	0.306552428342628	0.803126100604617	0.272074535056518	0.264654974991380	0.369215181460720	0.901235989573067	0.948492053739517	0.707622197853258	0.735795334596236	0.603879013823745	0.865007458742190	0.939244465053601	0.941017808416377	0.820818087023315	0.841074330172779	0.954896972581561	0.805159593160855	0.927439286186433	0.966990116229863	0.982382025538379	0.990554696463766	0.970966868948798	0.979103373851107	0.972224996502326	0.980708630387640	0.973080301091084	0.979603967435559	0.999725330436276	0.995775740568187	0.998025142521379	0.996453661220114	0.999753052270318	0.994314661591999	0.992273917822387	0.994602418074999	0.998935928013916	0.997635474411575	0.999095644941695	0.997083889556459	0.999331075041407	0.996289303051089	0.999264479383934	0.998595591722162	0.998764796346537	0.999406027144625	0.999633862288437	0.999999328113264	0.999694357763938	0.999977029384471	0.999567761927391	0.999515979668666	0.999205945908265	0.999996670017851	0.999213215713107	0.999900125788040	0.999761882476691	0.999568487084171	0.999894348762404	0.999706075863794	0.999916117582618	0.999655244435907	0.999893827280489	0.999499829577504	0.999316726387992	0.999891410874461	0.999957907218459	0.999650893598351	0.999983282393160	0.999946196380631	0.999824607919323	0.999708446097946	0.999980259578320	0.999878510900349	0.999915545024372	0.999655015920166	0.999953193717081	0.999768386854105	0.999722656140113	0.999959462573350	0.999836112494287	0.999876727343639	0.999557357145122	0.999935472693366	0.999731527655688	0.999734242031950	0.999981396214620	0.999996699768822	0.999854285726269	0.999836704377083	0.999724875791271	0.999839440197553	0.999946818637785	0.999979471074612	0.578371012030497	0	0.485047318379789	0.434568068907259	0.428506734062182	0.549815447210192	0.358262919913977	0	0.468081350900511	0.391031531875265	0.880464185673831	0.350596692748533	0.236186259344550	0.471562936330218	0.792726737909236	0.201643489192888	0.647284680192212	0.853178414399768	0.997428135294912	0.983060743806171	0.919030830639944	0.789927563588307	0.859895101288447	0.942136999099328	0.900520298761723	0.956354068016375	0.968770026922864	0.988714815143885	0.950390036572652	0.988946386972417	0.984610818418757	0.981284488352144	0.998487970741034	0.999753234305009	0.978063584881547	0.994486045954968	0.991305328084203	0.994157320353045	0.996632645292714	0.999474801407398	0.997567425624609	0.995792320613603	0.993839036541933	0.999360717139540	0.999276528529381	0.998313664695861	0.993979822534144	0.999718108608519	0.997824672825668	0.999554697094273	0.997722945609856	0.999986648111200	0.998782728999788	0.999610928255226	0.998263655765814	0.999265541162248	0.998954092519183	0.999797208741587	0.998527554691598	0.999745925075582	0.999564777829136	0.998963616886242	0.999945667244020	0.999771817698680	0.999289764133159	0.999731783552309	0.999905022307646	0.999184600351097	0.999641735414211	0.999915340723924	0.999241699462597	0.999899528324758	0.999759490648076	0.999971239632482	0.999701797879845	0.999925431463738	0.999950219797166	0.999877521112975	0.999928311753907	0.999950112039018	0.999903238606264	0.999920031549152	0.999932861723988	0.999356351947533	0.999857221280759	0.999736745596142	0.999946692101293	0.999824367877673	0.999709161664366	0.999838406545682	0.999822896566667	0.999664734709485	0.999991691634657	0.999790890191691	0.999996090142838	0.999930095619956	0.999934513766411	0.999854085972918	0.999820045573930	0.999860818416111	0.999959382995764	0.748869435429486	0.537003910930141	0	0.487883830941817	0.433128718954998	0.388643691729641	0.825631661838706	0.329047656828714	0.397091446415435	0.113178179446186	0.273168799021233	0.909389242839420	0.715179402015379	0.237942587062862	0.732592791297311	0.394582272110021	0.882913077491057	0.584502632059845	0.979817373021780	0.867287761284187	0.923961155748304	0.981236360315359	0.896484970041739	0.915105983881968	0.975913959843625	0.991461357720856	0.984281678571223	0.855629329605350	0.984980533010178	0.961179072024956	0.973676178812799	0.990597618085085	0.982516569388198	0.994937816398435	0.996196670294055	0.995211433784437	0.996060502645071	0.991774330531973	0.992261064539002	0.996515029657241	0.997818636700661	0.992668309266632	0.996812091686098	0.997230144260490	0.997182411855512	0.998221285272266	0.998472588325580	0.999487576725911	0.998672222454188	0.999126236638211	0.998958748731992	0.999396436012157	0.998812924264580	0.999276448344045	0.999588534062019	0.999764449638784	0.999714116488240	0.999591041791388	0.999931563461550	0.999451559676414	0.999779968579435	0.999054490477172	0.999756850229936	0.999681083539263	0.998929760644787	0.999818090046557	0.999323640427317	0.999555506492039	0.999658944270008	0.999948812994917	0.999904255112522	0.999530767395706	0.999854635088103	0.999844490096731	0.999590214589897	0.999937745052480	0.999746434011736	0.999711679155466	0.999605280909221	0.999612153421645	0.999521208682661	0.999927589969801	0.999931540450865	0.999863874143514	0.999915036644396	0.999910448000638	0.999621717623767	0.999939690696209	0.999833773128862	0.999929150955253	0.999759922112159	0.999877024353115	0.999989085645557	0.999712001433680	0.999974601265660	0.999880170847331	0.999866797516873	0.999993106714789	0.999765746755196	0.999828728291861	0.999856869578885	0.362492681307686	0	0.489964076059402	0.442654959124407	0.689626552630853	0.712215336096494	0.357791650556430	0.0109670217796999	0.692437429571240	0.600652259176712	0.575708569439340	0.670732780011965	0.886244840572757	0.228011709325100	0.783166219203953	0.724175461265752	0.837969669877483	0.665026520301370	0.655254618053757	0.862218757845004	0.879247418422142	0.800263120792887	0.958196482670665	0.997135701186860	0.982662961654330	0.953337207021861	0.909510743365728	0.886133266226940	0.973763799206979	0.968456731113346	0.992811355398068	0.991549761437935	0.997367561207487	0.989725429234598	0.981608606493046	0.982128965181339	0.993155515317533	0.995737977937662	0.996067670663204	0.994276039935785	0.992854238526164	0.999112922176259	0.993668357695981	0.999114608613150	0.998897071570144	0.999228734673803	0.998355150533096	0.997992499418127	0.999993421013289	0.999533436013459	0.998957541065316	0.998838080897363	0.999434644192441	0.999871925435650	0.999486471034359	0.999773895684863	0.999401039712196	0.999880089696453	0.999546324493779	0.998865623562096	0.999945228133798	0.999428552855271	0.999909058615964	0.999730000424161	0.999677700957703	0.999749827445652	0.999449115726650	0.999840614444167	0.999789986389742	0.999359618271135	0.999924731495810	0.999732007756362	0.999668210043881	0.999468517785532	0.999935898482731	0.999833478120008	0.999668436382565	0.999992710767159	0.999758913420695	0.999625745296543	0.999904137104664	0.999373371193595	0.999822361885193	0.999984915033933	0.999862572294288	0.999798632364251	0.999967942106140	0.999901986988733	0.999988211784447	0.999847862828367	0.999932844003679	0.999920266858947	0.999947492649113	0.999995217947932	0.999837355378596	0.999944582404573	0.999947830723210	0.999901270408747	0.999885340989311	0.999791098952176	0.999945226457658	0.731834075110938	0.498687975802868	0.574009558644591	0.914806166121024	0.539259287859819	0.695673241476112	0.351284884311444	0.331135893546186	0.667798225924357	0.463961195028908	0.288684350257460	0.262574416001035	0.275338032878205	0.230392220600953	0.996761803874984	0.741377631811139	0.912157985769964	0.753624482869953	0.979061525334458	0.709005864553470	0.949492278117888	0.975622646935757	0.990327768980605	0.982370169231588	0.950001528255703	0.977454790227320	0.982460181578495	0.935078993903257	0.940966111662407	0.993187833980452	0.976202893343513	0.991233893513842	0.991404851104601	0.984081229497060	0.985996700916069	0.991497205024188	0.991072802903314	0.998771756633348	0.997981279499986	0.993639753295495	0.994192603239805	0.993563809830546	0.997963044475754	0.997666470156990	0.996683315791759	0.997645631539569	0.998600588212050	0.999153444938002	0.999829595393139	0.998913153129660	0.998974538562686	0.999097130509076	0.999129781486935	0.999268599876479	0.999465658572127	0.999868707133736	0.999978824037981	0.999859531545621	0.998318341546506	0.999721046966548	0.999801412351867	0.999567816103514	0.999407577661835	0.999334696243320	0.999636199255901	0.999729975647101	0.999943905150260	0.999717735187914	0.999844206674875	0.999971648716230	0.999942113875801	0.999732088636545	0.999591333932224	0.999887417416226	0.999800734960162	0.999872777588696	0.999763942906843	0.999919913751008	0.999631394313753	0.999843147536311	0.999896459642969	0.999904484896821	0.999715664044999	0.999949581331139	0.999957510404071	0.999851361031753	0.999608791534749	0.999807820182626	0.999850697559248	0.999926399147144	0.999814091274046	0.999870436712161	0.999884983671751	0.999770482020500	0.999819691959478	0.999987572465839	0.999958824800364	0.999898533317290	0.999861991063299	0.999880985772679	0.999984387253162	0.568487943025759	0.502436954887660	0.0694288765944507	0.444176280053287	0.445988291157734	0.396122345409254	0.367658443814761	0.353375935341682	0.365399828334664	0.301064820507097	0.485475423518712	0.898842414710502	0.403705880025380	0.996604798687151	0.951070586864642	0.515245304981011	0.664829527868955	0.877612904696691	0.902244918030076	0.852069375136190	0.863233870780532	0.907940136340594	0.909781179914146	0.892163814241874	0.995462606491314	0.909912480249273	0.975772905751829	0.973423084397400	0.987018132665974	0.970875691394396	0.991236669483778	0.967387537610147	0.981064807855608	0.994413299984506	0.968951310066184	0.990384996350969	0.997817600572982	0.993410774516426	0.999882097213410	0.995457605543564	0.999855439868511	0.997513317079176	0.997870863422289	0.998286352689866	0.999340222947948	0.995449134997639	0.998743142808221	0.998429071528170	0.999324686512206	0.998803764724615	0.999304288315778	0.999590051573541	0.999452371541958	0.999726022644427	0.999957658622332	0.998792590755701	0.999261414319180	0.998673745910414	0.999988657674964	0.999597251591265	0.999312580943343	0.999249559613358	0.999665106359372	0.999850622822638	0.999957481832133	0.999856990489890	0.999926634750758	0.999919748675865	0.999463063734202	0.999717374056805	0.999926652534954	0.999674327071771	0.999596682570827	0.999722562135364	0.999896728712780	0.999835451259619	0.999856287706248	0.999898325537983	0.999684039977901	0.999730951518195	0.999696912629877	0.999950628480576	0.999977808723589	0.999870312476965	0.999893123685330	0.999981612739425	0.999680611559687	0.999882057958699	0.999899545603211	0.999891710255660	0.999718599981991	0.999845420506969	0.999903835295957	0.999948621047479	0.999964668805624	0.999856704696192	0.999596925152263	0.999987789821347	0.999772979249149	0.999944465215209	0.999769480401427	0.328633075165607	0.537839769350092	0.477269325558751	0	0.425252521369453	0.188653335985070	0.381483420982821	0.685766608640730	0.496833003922065	0.286988776350252	0.289746351470454	0.271762441592113	0.507059971325041	0.216629272293957	0.700240317085579	0.732634786552011	0.681800934487861	0.924474801664728	0.940142334309990	0.928086602434437	0.911077817432944	0.866128746499415	0.944999072301148	0.904960312859802	0.948391453332880	0.977474553999281	0.997663186572500	0.999479412526012	0.990772670315806	0.934505994652345	0.983839072046652	0.992355765765067	0.998253589131344	0.991723628856488	0.991517559320380	0.995701578231963	0.997540720700216	0.995049201312550	0.995974543089128	0.999604269870584	0.997286888053055	0.998069809164611	0.999857882247195	0.996119463637517	0.997288402962454	0.995895043879635	0.996845021374620	0.996572142754398	0.999969037871365	0.998105407323032	0.999225577838290	0.999310058666589	0.999617946444250	0.997292889146709	0.998718432867042	0.998906841708464	0.999655513159474	0.999481958810069	0.999483552592448	0.999589134130149	0.999584868744243	0.999952923813492	0.999832722912240	0.999574671167252	0.999906066337270	0.999993869974398	0.999037439591727	0.999666706910601	0.999877140086795	0.999837796313274	0.999510447520076	0.999736675510261	0.999850052264951	0.999851995822845	0.999857330025603	0.999272601902180	0.999880135465312	0.999679581525783	0.999941291387007	0.999905748620614	0.999757311089015	0.999853108669936	0.999914731491377	0.999610669396649	0.999945329605278	0.999863045913534	0.999849697249370	0.999878011108804	0.999868837660286	0.999983918127311	0.999927121816443	0.999868083547314	0.999675835924509	0.999998224062668	0.999942352298293	0.999977783218680	0.999856739483297	0.999952663618249	0.999835051627199	0.999975815251633	0.999735498897575	0.569892758034922	0	0.414779169803774	0	0	0	0.358996620640891	0.141547764140680	0	0.550181853874854	0.696531322991728	0.875358714181061	0.888630455091146	0.801827308989813	0.984480245603407	0.919051625227870	0.691207225634136	0.812179496719203	0.901116349216176	0.576325670061810	0.983909385735967	0.950517442565436	0.991515273738797	0.934756480392496	0.993343235292248	0.995785545711675	0.981320436216459	0.957608591119930	0.994823192602711	0.995201490990407	0.974483379816809	0.991476021406461	0.997095069203287	0.988114435125935	0.981038773292222	0.993926064253657	0.999364890350010	0.999067408331316	0.997391866153925	0.994306862054789	0.991363329163338	0.993793563216380	0.995793311927788	0.998991330123102	0.999745949977434	0.999125655588055	0.994002204065119	0.998594532353790	0.999839091982713	0.998198609956521	0.998910652829218	0.998774680433835	0.999194100777733	0.999152874274568	0.998890860664117	0.998046217428567	0.999935662155103	0.999603256861665	0.999547087992911	0.999341099102956	0.999651271259616	0.999122434469592	0.999674629594925	0.999222800750575	0.999518040145095	0.999677206808557	0.999763210919094	0.999861899748996	0.999911325945747	0.999643667282068	0.999973038537749	0.999835143076337	0.999528703618858	0.999379743791601	0.999963411106285	0.999318586578614	0.999669495653877	0.999978454141412	0.999764917962896	0.999934583254640	0.999867297660712	0.999979266984622	0.999938106543579	0.999996036714983	0.999870598550797	0.999919670292334	0.999801620638947	0.999647111423200	0.999967372922231	0.999969009623359	0.999868328091430	0.999913950762864	0.999905922455578	0.999917266197683	0.999909528778861	0.999836425412326	0.999999422087947	0.999800209820131	0.999945941529726	0.999999718146351	0.999961651897909	0.554270642664805	0.524581633845839	0.465609005793742	0	0.585250076490197	0.410043348855816	0.442626005573243	0.328270262895174	0.514451780065694	0.287716921908355	0.943380858125172	0.542567212122142	0.444953467831409	0.235217575963169	0.793621881405710	0.859754451467450	0.922972930346208	0.971929168485609	0.158426309950440	0.959362885726467	0.773215249191288	0.997892621131221	0.966848983012250	0.778004542123701	0.949470089624026	0.965427330024933	0.976788755084094	0.984955765302974	0.979939215674962	0.951129977604111	0.996767221992498	0.992293438142110	0.995240821763817	0.995733577448501	0.995899758241090	0.996619768950837	0.995042460023631	0.995596217363635	0.998660297588899	0.988251891858547	0.997516418528436	0.997821064773673	0.997152503598051	0.995500349550211	0.999336992375127	0.995664437071463	0.999329453362161	0.999814489735244	0.999923707500583	0.998642321698317	0.996030159427094	0.999409488790537	0.999629197286703	0.998084662450907	0.999480760291300	0.998842904186494	0.997694213445310	0.999272133746003	0.999750519864400	0.998852593353103	0.999908747314350	0.999696595902755	0.999507707558958	0.999844537528222	0.999865021651197	0.999769122701137	0.999349426385650	0.999929279424709	0.999545749092361	0.999668283976017	0.999811086443159	0.999936182171000	0.999591181819563	0.999637539172610	0.999733519308540	0.999556461202339	0.999866403546424	0.999917342056359	0.999648350398741	0.999983571661768	0.999628225938690	0.999767138760638	0.999818924153461	0.999915880859926	0.999932628491338	0.999909161612301	0.999918291492471	0.999831831146523	0.999945539935289	0.999947016008306	0.999821137428349	0.999596043346743	0.999920470246995	0.999950600287935	0.999917791466385	0.999953131953355	0.999979808100637	0.999956279098798	0.999898996560864	0.999941318247092	0.999933654007682	0	0	0.0989609427612771	0.0635120616065433	0.432290714947904	0.388741669242420	0.541976296041383	0.337960961041510	0.281778034001388	0.289875711434300	0.771221708345687	0.426488999215898	0.253907223832192	0.231740422710392	0.944758966837478	0.620350531890252	0.721558104675302	0.695950505077865	0.937916691553694	0.711509324851866	0.965212817780946	0.853839538375345	0.673769651756834	0.932587551840185	0.998869424427712	0.943666213117103	0.984999159944069	0.980281875781630	0.972318927821431	0.976706093588491	0.981374153178520	0.987210126530810	0.976553168572423	0.993093799982160	0.984122513588188	0.997556574026419	0.990303715789693	0.999948771130789	0.998472059641937	0.996398735196432	0.999612941673007	0.998164386639641	0.999501621924843	0.998678023908174	0.998876379030620	0.998858417138085	0.997710682273740	0.999159198441200	0.997310713694938	0.999942153213298	0.999640653361960	0.999816816255430	0.999812980326790	0.999204095366382	0.999549822023882	0.999732673201493	0.999726374397852	0.999991350174119	0.999243455847322	0.999196937669709	0.999334977326249	0.999861128621014	0.999727443525565	0.999771380188663	0.999220471130943	0.999909446530948	0.999795004614808	0.999264587790740	0.999900629266919	0.999964897520551	0.999834969540634	0.999599484616729	0.999434347690474	0.999727070654095	0.999965709506092	0.999720239381186	0.999655619052525	0.999957196613986	0.999779982923401	0.999764191350671	0.999597174299793	0.999828369991831	0.999785469431740	0.999863225387363	0.999746433284427	0.999913605623130	0.999831392426193	0.999871332436142	0.999784855077144	0.999842602944600	0.999936780582439	0.999806135259393	0.999857593079744	0.999824541245182	0.999859744830067	0.999946811074304	0.999950004691063	0.999994080673970	0.999946763118080	0.999908854397652	0.999855751711712	0.557090189762671	0.513122232708396	0.622535735698548	0.391755807341625	0.405555686619159	0.466000231836090	0.558319981793347	0.0232368120322775	0.663223114701478	0.299742185721511	0.379573991094103	0.462719912120502	0.788380258823342	0.537567754003326	0.799774478481184	0.999592705767928	0.982578630553100	0.694963069353688	0.954028146181412	0.989629217292582	0.972769094492341	0.878030685706608	0.988445150277969	0.974012684467553	0.953534995155312	0.970443858835992	0.945561900656957	0.893670936439477	0.994152207060380	0.997594901200520	0.993591315195212	0.978373742770135	0.996112523692736	0.999072949488076	0.997536485036440	0.984338484385071	0.994365375220815	0.998092327698627	0.991957639291601	0.996012938720580	0.994120303150207	0.999214934809608	0.999419152214340	0.998691623147345	0.998301484211206	0.999509704942017	0.999024078710476	0.998504422160216	0.998334707735108	0.996948637088052	0.996871613218768	0.999980683325503	0.998643763436099	0.998086166986767	0.997816747593577	0.999808369334548	0.998809978523362	0.999538381191691	0.999717529197612	0.999905062965275	0.999997560600665	0.999388497821317	0.999516166661784	0.999838183408560	0.999723857830061	0.999972018969375	0.999977549901529	0.999963890880705	0.999651052436775	0.999445257894710	0.999273838787384	0.999868492139991	0.999642040345452	0.999516132861509	0.999590504409176	0.999896259507715	0.999957168324079	0.999868930006005	0.999805405089985	0.999799802281912	0.999587739043935	0.999917525513345	0.999788414494853	0.999991596679048	0.999757931131495	0.999836161019036	0.999861772483744	0.999815000433040	0.999884651499933	0.999662347170861	0.999919072810035	0.999981919287629	0.999909625476989	0.999915172783183	0.999915162700253	0.999815737367592	0.999920659998108	0.999828726746331	0.999919261206841	0.999939422977154	0.999965090678776	0.576930980551708	0.534277096877880	0.471001880424422	0.487865884162572	0	0.481223862813518	0.363933182941549	0	0.369649479538141	0.746715842253481	0.281019058266153	0.637736535098878	0.601347717730040	0.227331775575299	0.224135905748761	0.956857672041645	0.914063609527877	0.940138646268520	0.957114025333363	0.814748076866681	0.783659371129393	0.888090667445055	0.952612926717966	0.923908214238732	0.988193950309390	0.936149602708444	0.982396992077113	0.965897860531382	0.975661644849598	0.987170322722926	0.990623433324968	0.991894338986307	0.998264020924357	0.987902782535004	0.997117158105815	0.993284872679315	0.994663088653720	0.984842604951550	0.999552129682794	0.996572914646308	0.995953947572060	0.999323694958353	0.994555074482757	0.995805626500725	0.998599679762830	0.997564600691702	0.992638603260787	0.998971264724795	0.999484858122804	0.999743459732087	0.998695763762606	0.999638310834863	0.999238356887906	0.999740293237450	0.999860775854059	0.999384263202764	0.999326881409784	0.999573164175906	0.999646067245812	0.999706543006316	0.999801590486427	0.999660859886739	0.999570839665518	0.999611446582871	0.999925838927214	0.999493153573730	0.999580649632599	0.999928049991040	0.999687400487998	0.999846908223334	0.999966148243955	0.999853350902012	0.999714239796027	0.999743306700015	0.999921324966110	0.999997887020654	0.999759498960921	0.999964309722503	0.999954586759032	0.999966400119813	0.999932572653595	0.999891575386117	0.999761492375984	0.999908345299769	0.999963114102264	0.999793434585174	0.999996060861465	0.999518424667111	0.999796340020534	0.999976842315323	0.999860792457610	0.999910274343583	0.999853280662110	0.999928005833614	0.999947748707740	0.999912743154936	0.999860872499674	0.999807627200278	0.999816337476195	0.999850956332970	0.999920278190996	0.547280219171682	0.133731657549646	0.497724143277222	0.437109505986903	0.422369456141355	0.440760886406309	0	0.352315056214750	0.910670527440909	0.298910772382493	0.674997545807632	0.990614684695469	0.238908392324782	0.754973555589106	0.221371704445847	0.964655488259302	0.925612172212154	0.976542661666439	0.742894103221200	0.811795834432634	0.909244678444445	0.973177826415332	0.972748134946322	0.994613953424192	0.940771011730182	0.959820146182958	0.978836199259363	0.968236053881718	0.964259220651521	0.985706014114731	0.975714993165977	0.980479044439902	0.987378329102258	0.989917649701333	0.990133935002039	0.994402730052038	0.997720629560442	0.993193612078379	0.995405178537767	0.998997621612825	0.998863181119316	0.996201566136041	0.998742109783709	0.997792813503752	0.997073295031059	0.998861061271425	0.997694695583947	0.996787440975141	0.998131662210219	0.997702338304813	0.998759001250196	0.999473856252236	0.998163337620211	0.998985587563888	0.999537253997409	0.999161240122452	0.999068738285545	0.999429425974522	0.999806114064271	0.999887579316245	0.999627305460856	0.999741726499014	0.999815293751165	0.999640035470334	0.999627197740836	0.999766359122707	0.999779458270346	0.999941356808255	0.999741616777287	0.999899116731200	0.999606120438444	0.999386510711405	0.999847864886326	0.999901492422894	0.999643908570809	0.999794797923066	0.999806091917411	0.999760732348972	0.999839851714104	0.999857969089384	0.999819379176280	0.999933511623551	0.999922862733651	0.999924304336702	0.999823035086736	0.999935455813499	0.999920972471282	0.999928798559541	0.999943091672018	0.999946280161009	0.999902815762971	0.999972805932237	0.999981559699236	0.999864082836018	0.999944725930487	0.999755541378133	0.999823895390236	0.999997721467361	0.999998109286639	0.999914468201994	0.999934024844518	0	0.734739563565449	0.697369392559326	0.0304820047663813	0.421535321554051	0.198246096769889	0.456054122325098	0	0.322995309367680	0.893561661488780	0.269837241986267	0.273138161329056	0.246571159302201	0.888200451731024	0.946749033178475	0.930840110740302	0.185306022895502	0.870396567995819	0.887509890845770	0.976605012318728	0.851662271011890	0.954266669060564	0.983532723239480	0.956902593330152	0.990116057637962	0.994618692400340	0.983327242494497	0.968889225732844	0.970663653045662	0.982481007667974	0.994394120323466	0.977787781773057	0.988974197028304	0.998746121521634	0.990750992378969	0.999825631303537	0.984583703940078	0.999071828489779	0.997149856950210	0.995990948500149	0.998160477522029	0.994842354442532	0.998725426419215	0.998668545394830	0.998692247652265	0.998974361179494	0.999139563252093	0.999384314604449	0.997580289848117	0.998346842193912	0.999998524565708	0.999904335036963	0.999396433556160	0.998985368351082	0.999136179791935	0.999948342776941	0.999987986383043	0.999974926685908	0.999373496553669	0.999931462978725	0.999922989211631	0.999813476908997	0.999766803889049	0.999563620599747	0.999923778823349	0.999716010377611	0.999609404639893	0.999917720020792	0.999989036725867	0.999933097520098	0.999736625373818	0.999948174916760	0.999967813244942	0.999505105565285	0.999869481519613	0.999592875143072	0.999645073929099	0.999646639335543	0.999903431539037	0.999923607123760	0.999881803006030	0.999935259723921	0.999961836643870	0.999735594230310	0.999944073819832	0.999904211299303	0.999769959988778	0.999992834069298	0.999903475917551	0.999901269456402	0.999754532843282	0.999757386608776	0.999805516934721	0.999822382313987	0.999737682904636	0.999872068329960	0.999794520326945	0.999944695947023	0.999919642533382	0.999901673322382	0.999830615452380	0.569368895295834	0.645522129135764	0.155083658402676	0.442072429372924	0.875939801707034	0.388299540641190	0.824443726333794	0	0.517905443239431	0.305904323845094	0.826506586770259	0.260595354989005	0.777241915806166	0.877319868711031	0.225421540286786	0.942571080376110	0.891979521682101	0.648536671403816	0.959580616025623	0.964687741392644	0.890639996590674	0.931263552517027	0.949170777407571	0.847067569599901	0.950738840535287	0.905489702411994	0.991566409946300	0.966248429728168	0.980972517798676	0.966844421072128	0.987907791305751	0.970309000513952	0.986750326363125	0.982031635965991	0.993199252564794	0.994856922032509	0.997907492202477	0.990949230196761	0.998967449029665	0.998832516587194	0.994155674313501	0.996703675019551	0.998028730592143	0.999365885348914	0.999388251753831	0.999928756056427	0.998601159733640	0.999505123140475	0.997783512760892	0.999680483267345	0.997298326363170	0.999715999692447	0.998800002781426	0.999917662828158	0.999559884589513	0.999457712262862	0.999196725578429	0.999338808243560	0.999510420838169	0.998928528568155	0.999273198013041	0.999633446031133	0.999409047719655	0.999871590571871	0.999806257121425	0.999734981421901	0.999460401802371	0.999816793755451	0.999999909423307	0.999958758098203	0.999795019617766	0.999740165500793	0.999542884952279	0.999373900403983	0.999765785030544	0.999900135225331	0.999854351457383	0.999603666187572	0.999980656506506	0.999536144208305	0.999854287723727	0.999624063598577	0.999757316138000	0.999990905966354	0.999959319440283	0.999774226497831	0.999957467774491	0.999982206753105	0.999978287542535	0.999744685634184	0.999769126685514	0.999964064308064	0.999911477838286	0.999968233297177	0.999869673101797	0.999810466067348	0.999999965890720	0.999807103035605	0.999761303630718	0.999944970931717	0.999861980846601	0	0	0.495316813424239	0	0.408204709775394	0.379052249026084	0.407428502296900	0	0.313521059653216	0.407455385223689	0.285266860710159	0.759464039967077	0.835449036306448	0.890083187523020	0.876254607270404	0.800974561279282	0.924205317156062	0.872157224125089	0.899492797890704	0.753090381383392	0.592574429736890	0.941003514174788	0.987477569536651	0.925532966446947	0.910454488131313	0.993967818902182	0.981298628559293	0.996961114639202	0.999974739956982	0.984889695374152	0.980344631314309	0.984097713394425	0.988642265888749	0.998746999594072	0.998099588225630	0.992027008921981	0.999385557698079	0.980377398625056	0.999376624875708	0.994198139782640	0.998285374640783	0.998203978164547	0.998697270091347	0.997339226621126	0.997768005025240	0.999929371398434	0.997430932962838	0.998553338216355	0.998638538918036	0.998424490537671	0.999845007078310	0.998408595457548	0.998410318946279	0.998859474075199	0.998987428773597	0.999835866780893	0.999667196969557	0.999209091808544	0.999159943850198	0.999957640886339	0.999761909146874	0.999096383594283	0.999948078650951	0.999463768745261	0.999598560930394	0.999782391854253	0.999776437802010	0.999522341912953	0.999974157081374	0.999719628957588	0.999910128093542	0.999649900449917	0.999745723311863	0.999699716618462	0.999950561345992	0.999712361369859	0.999654385362193	0.999599557696313	0.999638300314638	0.999765122029053	0.999981437384498	0.999993425508152	0.999905598452339	0.999924586892042	0.999801833164254	0.999817088284212	0.999801475988108	0.999944430080725	0.999957774003964	0.999911920915309	0.999847216196462	0.999924261942814	0.999876041579154	0.999897562839870	0.999941676982815	0.999963853790317	0.999886135674142	0.999903848885806	0.999846354510449	0.999941426865993	0.999921592767163	0.560498631239072	0.453341648259969	0	0.974613414331053	0	0	0.370857921774306	0.531964762291741	0.328270367688528	0.291969047186092	0.799462910710719	0.791744837223438	0.692554630477106	0.987462999551394	0.742288632609597	0.653252129027974	0.949902925517358	0.764283169222829	0.861719301662071	0.764543698506319	0.823407622063143	0.869904517473594	0.943531048383517	0.926007701069358	0.942690971608934	0.959964045669004	0.934814639074885	0.990434423189157	0.980284924738666	0.982506686221853	0.995569283984921	0.998293350038482	0.996446860622393	0.999441303416789	0.990643143373908	0.992119904014112	0.998662984640068	0.994073824209888	0.994159670258975	0.996985438895388	0.998819387291761	0.997047525126161	0.999907766323720	0.998056206500239	0.999640606446812	0.997133973387393	0.999276486357989	0.999801734011404	0.999662144975093	0.999819128772329	0.999191294831961	0.999794976110722	0.997963871679304	0.998916685531832	0.998819653166965	0.999702762193144	0.999328622304611	0.999608048283010	0.998894299607115	0.999519624092713	0.999785179277963	0.999442175554137	0.999027683637951	0.999867566949829	0.999551124718552	0.999825795104608	0.999910637618469	0.999694740749720	0.999502151443798	0.999855475451407	0.999488378876225	0.999772157427375	0.999983040101089	0.999938849688755	0.999923483727333	0.999562897084886	0.999654666984845	0.999687352013195	0.999752746367785	0.999826687076468	0.999974857556389	0.999887555885908	0.999868843290443	0.999841605536908	0.999918769427014	0.999718650444482	0.999881512497667	0.999866663107907	0.999840630445053	0.999812248652953	0.999948341227585	0.999984149687248	0.999941151974044	0.999958672910784	0.999964568127034	0.999814119269131	0.999823805944467	0.999965192954046	0.999884603777652	0.999931540397708	0.999837476297065	0.466565360647205	0	0.978753544739253	0.507604151082644	0.422107820335113	0.410869005038964	0.517049392329949	0.488154746140144	0.947551067914911	0.303923062211973	0.800147881891918	0.512238101051984	0.496564050266554	0.996318326413245	0.867236239043708	0.848426307049542	0.836596820681429	0.576640037730598	0.816778587810389	0.784556166306217	0.922786733548076	0.939096417968425	0.894219681252714	0.999710146258949	0.954638960097308	0.911716598478981	0.939329315973345	0.970539352768567	0.962396814634218	0.971055669512067	0.991219180178687	0.986410815783176	0.999991659014600	0.981165094132101	0.974211820009308	0.984477124398770	0.996644784492119	0.997908323942798	0.995976176908223	0.996653131406262	0.995276601035590	0.996312816481493	0.999595502982550	0.999431657200583	0.999121024196320	0.999406937494449	0.998215973534284	0.999475942515001	0.999529495184434	0.998231270987854	0.997881267357552	0.999310266681725	0.999746216026550	0.999111086632049	0.999928706912585	0.999656399494047	0.999493517264851	0.999795938247961	0.999947812370576	0.999767039305346	0.999658024979268	0.999720897639097	0.999437138447295	0.999913952925454	0.999560456198214	0.999688776495967	0.999818670396563	0.999881542910395	0.999568012649895	0.999862970027878	0.999612066494039	0.999724130044464	0.999908155302721	0.999613094111520	0.999944375829720	0.999838624694022	0.999774318708576	0.999648732998689	0.999875287785598	0.999887608828080	0.999886522910113	0.999931818990635	0.999830116934182	0.999958786441624	0.999855958340718	0.999728134892365	0.999913806725759	0.999993882058502	0.999949568238662	0.999907284722499	0.999973023509773	0.999854520966028	0.999816712327715	0.999821505482207	0.999890157400246	0.999930205743155	0.999980391309479	0.999967711012467	0.999960848042842	0.999987527159979	0.999862707852269	0	0	0	0.462744323761055	0.669595203877463	0	0.734827840503342	0.645921818408135	0.331011225296643	0.677560880275405	0.556220125516810	0.903232502345951	0.555758805428819	0.531283641859886	0.217930991152060	0.793488376139317	0.673666712626389	0.883718389373735	0.881349861567372	0.827440281879485	0.887618389239586	0.974683068348467	0.816526589173687	0.951057947929034	0.977611523638989	0.982755648342646	0.938448832354908	0.983841351148043	0.985327759507479	0.987213094158384	0.971377054775876	0.998864156390682	0.993283432112372	0.994673536192746	0.993564254328801	0.994549541490943	0.986142161255016	0.994078993135401	0.998860498541136	0.999762082327143	0.999543218521316	0.997011506469235	0.993559146929479	0.998963016646508	0.997246150795046	0.999532948872015	0.999291146924284	0.999611000414609	0.999563484512057	0.997578416913903	0.999817011591266	0.998295436400097	0.999426619625979	0.999950345867970	0.999459015819918	0.999253898184653	0.999352326603975	0.999846685951975	0.999922028978648	0.999980732131106	0.999286352443867	0.999721129799931	0.999904230288543	0.999489911100797	0.999508393450857	0.999947641165853	0.999930715031999	0.999274869351086	0.999859857973220	0.999787760968074	0.999891020625005	0.999795154249567	0.999681744026162	0.999703638994661	0.999566058787566	0.999718422178136	0.999988043904566	0.999838457659464	0.999821758063250	0.999790467641535	0.999924693604815	0.999791598247380	0.999676475608391	0.999913757660077	0.999982168185226	0.999984639034839	0.999975926581589	0.999847289146644	0.999978237676010	0.999959026249266	0.999939170434016	0.999975800264762	0.999978986732682	0.999945012837421	0.999770832830706	0.999782125983159	0.999904946538048	0.999990029166214	0.999993555946872	0.999826187736789	0.999991200511125	0.566138683646958	0.0251874985900514	0.470680272444946	0.459319733722969	0.407758623727040	0.656288168556810	0	0.0976395093648775	0.336704265543574	0.637752259456347	0.284620741256154	0.268109658095475	0.395327631534211	0.719393152517663	0.905324295530787	0.201462917160970	0.965959033955775	0.830234907824436	0.974595351969521	0.860957982376384	0.950335952296802	0.970787498754342	0.994366621732702	0.919129439535737	0.964011230024461	0.875527474294656	0.911423287831557	0.970210185533878	0.970345091473637	0.954483998180814	0.995302989073027	0.990376928148845	0.987155300164319	0.984482818111481	0.995290910350282	0.990400797888178	0.999685904845225	0.993393572278630	0.996811716921703	0.999334169204837	0.996841989551978	0.995477530132298	0.999483949501730	0.998697178732829	0.996326335650102	0.997272164063725	0.996468112919694	0.998616604347037	0.998122778361098	0.999364070494674	0.999823796647425	0.998822842488220	0.998717910185787	0.999437849157571	0.999950038867209	0.999632794146060	0.999775309667886	0.999709501911091	0.998417014216126	0.999359776541096	0.999504032740078	0.999257155016130	0.999371436007036	0.999550673130633	0.999877189670188	0.999854133688684	0.999603842342568	0.999646570044619	0.999848756063351	0.999763112774334	0.999297578232787	0.999831614867345	0.999897188462062	0.999873483909720	0.999905930964949	0.999990577682988	0.999636520927760	0.999937619252277	0.999935008911183	0.999888621339762	0.999920637493247	0.999931090003507	0.999878103366802	0.999756202170427	0.999890027744749	0.999960937208447	0.999746820981453	0.999950160612061	0.999969238541688	0.999954728801997	0.999942423507325	0.999961674897946	0.999993188914687	0.999954864680566	0.999730099313451	0.999901118307681	0.999964108989838	0.999984871994618	0.999929216889376	0.999908796410131	0.999908645761992	0	0	0.972637609460587	0.463494449872946	0.413921901833421	0.641446283715570	0.391267473955564	0.956851061855463	0.412825354100974	0.306999606503755	0.565882648014228	0.922819699895756	0.905112292787700	0.986932878800605	0.688452595099745	0.852912578959163	0.746136943404293	0.987591317624276	0.921149078571855	0.871131816299820	0.974295715051282	0.997782109453791	0.871222259069994	0.957610733161282	0.958454681131783	0.964110550882565	0.997724474925628	0.993081397394759	0.970343267423508	0.984707653588158	0.987986505127471	0.982063055084037	0.999438867285842	0.989099847876993	0.991625068186329	0.990248668722187	0.999657088783870	0.993737973480479	0.995911091709663	0.997294580407718	0.991772891922009	0.999323707027808	0.997984676317504	0.996578498911512	0.998510914885980	0.997746388068406	0.998408233566953	0.999334163671005	0.998777562768451	0.999162939553871	0.998153041134907	0.999013264777620	0.999170156876472	0.999567650566975	0.999726714502218	0.999986597154103	0.999524803385702	0.999630870903004	0.999554990646153	0.999795189520436	0.999643176792264	0.999909093724791	0.999622806888172	0.999945067977575	0.999683645640951	0.999907013450624	0.999916728616310	0.999451054321881	0.999820125096350	0.999727392947167	0.999580914799834	0.999948705150381	0.999651760436697	0.999756749776023	0.999938310419042	0.999936779441514	0.999993651487477	0.999635052651750	0.999925675160278	0.999983362744626	0.999887110886726	0.999987530119435	0.999850920367989	0.999941443461911	0.999936593498719	0.999823726360346	0.999974142878394	0.999992998728926	0.999908466783902	0.999989867616934	0.999955289263612	0.999912967693075	0.999904309924980	0.999959705748374	0.999993473009736	0.999892269449381	0.999965483325706	0.999844602959369	0.999894773850915	0.999911357850065	0.999945394308194	0	0.548674697852197	0.215396556703304	0.450884812435516	0.431062447959694	0	0.819357989423153	0.792914733213805	0.571369117755588	0.991128950989462	0.671419312917065	0.683732822648466	0.343299237658654	0.772129172854907	0.823901558324237	0.420444843472135	0.685004157496749	0.926693424987403	0.933389657682679	0.897567753988352	0.791646206973396	0.983134023732765	0.991038767315517	0.974146187205237	0.972523650601186	0.838303673851567	0.967796773497858	0.982600377369186	0.994559377765963	0.963903046531152	0.987541491891138	0.977958301375444	0.988870138201768	0.995021657671332	0.996352599843094	0.987998576109879	0.998821595250172	0.998747078885940	0.998717721644170	0.998607688791974	0.998829095097371	0.999880860001120	0.998053074853784	0.997663923511504	0.995382255416482	0.999437208482933	0.997237619523822	0.999285557264405	0.999316969472272	0.999075016404031	0.999306284907716	0.999983440715537	0.998816416127335	0.999184719764234	0.999878777766571	0.999692643116251	0.999472890303531	0.999138894977992	0.999493813261167	0.999893198083237	0.999224924314763	0.999905813062022	0.999809288018640	0.998668384385390	0.999832755977691	0.999632189870983	0.999959804096922	0.999403091386963	0.999759642325886	0.999828657549031	0.999734155163825	0.999473721171937	0.999950052341278	0.999756385424749	0.999524417242272	0.999610297222195	0.999838806985893	0.999929067639592	0.999704289319342	0.999812868768219	0.999879318480068	0.999592487608998	0.999842594427778	0.999999623756311	0.999740830243114	0.999866382256507	0.999755853600644	0.999999297374245	0.999836152381994	0.999840634070149	0.999875778056856	0.999935043780495	0.999881502685060	0.999976355091896	0.999987080212819	0.999929445375076	0.999864436273942	0.999992169960308	0.999878828349174	0.999936252737791	0.999861271585649	0	0.507876217628282	0.691149701784605	0	0.427760073186945	0.544462214795097	0.278241219333830	0.994769782540884	0.441457192176351	0.746015810320879	0.776555770146919	0.582451365878252	0.873262661554173	0.931541777289188	0.948187180098462	0.933087422553199	0.741460043041297	0.917104500015447	0.925173494188797	0.890742352937433	0.952445063308058	0.989733474855364	0.950733858417376	0.984894280483949	0.933509662780199	0.977364977016421	0.895980876559631	0.993754192742891	0.985783486454878	0.999644780076292	0.992847202405555	0.978654030170153	0.991599153418195	0.997283018522093	0.995464108748383	0.994489883649764	0.997519203629356	0.997450213886623	0.988459988043603	0.994426291799086	0.995776027253687	0.999924630901962	0.998315592864935	0.997450338934419	0.999954964447134	0.998693834214761	0.999516627220095	0.998939341675509	0.999929853318598	0.999317821449595	0.999616010119920	0.999297121465509	0.998801060377694	0.999672330980488	0.999018561044348	0.999936656773560	0.999749229587037	0.999749115643258	0.999940239166975	0.999066351984186	0.999712220218911	0.999226122123281	0.999117468792286	0.999460524856932	0.999707452415958	0.999978948419360	0.999446257799363	0.999909903123743	0.999777899698639	0.999870212418494	0.999490739712320	0.999879177860792	0.999971429878442	0.999990530004976	0.999570145418151	0.999696827568225	0.999800630749939	0.999812183940823	0.999904743327681	0.999959436787259	0.999783277961158	0.999962211394264	0.999966942963526	0.999716794466230	0.999963887792855	0.999947405454213	0.999998701348317	0.999978110566824	0.999909153597558	0.999913205322180	0.999929243600048	0.999954875243119	0.999974248821415	0.999948385869505	0.999784855355271	0.999830927114782	0.999809590331358	0.999976195461284	0.999872318822666	0.999948451806243	0.999910685954982	0.546544180857567	0.436828528534203	0.472774044605388	0.705841901284057	0.0147917629623082	0.413067984857043	0.359949243650197	0.347875810821352	0.588739246476664	0.658152433271031	0.305476139219633	0.395931699856565	0.843205790736557	0.227330397580795	0.769835606958174	0.856524771467413	0.695635276123827	0.785596496587863	0.931922982823963	0.894671141935548	0.934883389442029	0.910816709459661	0.928428647224402	0.961539861432311	0.929946525654698	0.952589712346303	0.957577948677798	0.983894603426473	0.991735998036927	0.972797680567969	0.996051473603292	0.981210351309344	0.987460509872090	0.996643650836952	0.993620364528289	0.991467685893108	0.997179856694623	0.991727869444976	0.997719761659861	0.997587073896091	0.997363990316123	0.997442342102272	0.997895218861113	0.999470173399334	0.995600064359340	0.998946415662493	0.999682824261362	0.998960439349671	0.996998273022831	0.998330020579826	0.999242763450190	0.998887386711603	0.999975271010835	0.999701595959233	0.999907682894754	0.999520406941535	0.999481075202553	0.999263897458721	0.999489216715792	0.999953277365261	0.999946516387531	0.998983269671247	0.999930972778549	0.999948444620060	0.999628811822038	0.999089741324205	0.999780835222118	0.999920416455555	0.999765782398392	0.999570739343787	0.999434986770168	0.999905387126330	0.999833951410956	0.999568819662143	0.999742184757363	0.999927379655458	0.999859486533547	0.999569630197237	0.999805162768301	0.999798561070206	0.999989576235522	0.999807071967279	0.999960275216343	0.999820816056928	0.999834260733278	0.999985304600771	0.999933497247360	0.999907203668407	0.999848472221326	0.999810604971145	0.999758333384944	0.999771072657018	0.999898411275304	0.999886559530241	0.999897126151319	0.999881782961754	0.999903142268687	0.999884632295281	0.999994256168865	0.999958359738140	0.999975844566646	0.559100138115496	0	0.494206781126445	0.453562973807749	0.445285032936343	0.435591527459007	0.997804468537992	0.863497232618314	0.889502661254050	0.855298857648016	0.278899432005729	0.936376772338681	0.510057481244979	0.990379765865609	0.903934982246438	0.939640749397309	0.791680497726040	0.702846781936362	0.667464745901623	0.940554113187538	0.952638988040034	0.891794723208895	0.916977637727129	0.939558725613747	0.954500863028045	0.989553870287378	0.995106762987073	0.985697900327459	0.990280046542946	0.995917812924588	0.983986009724740	0.988630479120137	0.997368413848830	0.998340348572798	0.993152401739793	0.993515509580236	0.993858492840306	0.992077979572253	0.997955325349798	0.999734909089901	0.997043722225184	0.994865510677550	0.996051995635550	0.998612539670561	0.999285112826506	0.999922529373435	0.998168489388828	0.998444325095731	0.999697955557554	0.998347652257815	0.999788815193293	0.999550480007626	0.999871779801195	0.999819099266457	0.999954125122996	0.999825025751635	0.999981911041692	0.999967426994107	0.999346995096760	0.999448420141958	0.999376737536005	0.999684791772632	0.999417991518526	0.999789914107042	0.999780494538374	0.999857628135876	0.999631586557747	0.999711911871294	0.999907426422936	0.999753870033366	0.999750046179554	0.999525987216636	0.999521528203181	0.999771946117203	0.999901919636844	0.999992630372698	0.999973603144624	0.999826539201067	0.999955930655093	0.999950948548826	0.999801349448199	0.999985115967180	0.999784654673598	0.999910425749892	0.999923238204292	0.999896150893678	0.999859703471536	0.999879244819742	0.999697062427263	0.999780979415537	0.999873731367671	0.999847184346269	0.999963976782467	0.999996328849987	0.999858101446215	0.999950032742460	0.999921087945185	0.999928091448678	0.999938883043242	0.999980430353340	0.999883981741246	0	0	0.991958293884338	0	0.212189355236344	0.401792007689609	0.484685635355703	0.838079270105526	0.989407728728390	0.930304200334103	0.534285545423631	0.818971629235657	0.790342701096396	0.797455393980568	0.899898132997836	0.387766848676235	0.957728390951677	0.812723559990969	0.873415315508025	0.886467314521314	0.992851136373661	0.921313431147255	0.963680689112392	0.979379243830376	0.977862905888568	0.941946991145279	0.954267572592126	0.983451836067050	0.985145847652199	0.997554889470242	0.992516855150303	0.998377574384789	0.982795178878093	0.996909905930967	0.990105762981344	0.995534594713515	0.999806037275429	0.994825388132929	0.998487212056375	0.993849171022806	0.996719977163238	0.996730078258769	0.999181261296983	0.999743264428157	0.998602925925201	0.997143684350033	0.999557452642451	0.999725726050157	0.999559277962673	0.999393689838525	0.999610602074864	0.999221472603202	0.999738148438629	0.999517557299948	0.999281303946615	0.999884315609924	0.999890811523136	0.999601007406413	0.999917012996029	0.999808727432732	0.999289111056977	0.999976050045622	0.999421095292947	0.999900318720327	0.999676725840296	0.999711631600939	0.999899336292731	0.999741983482889	0.999921996094409	0.999715845176896	0.999787442484447	0.999903857469796	0.999604249933722	0.999625158255165	0.999904627749385	0.999880032652674	0.999939858769876	0.999713402103636	0.999958194107387	0.999884742399217	0.999837861841975	0.999860924962281	0.999765478027663	0.999936927168682	0.999841994376196	0.999944184750864	0.999935215460788	0.999985025554148	0.999908089413251	0.999950651062224	0.999988800148015	0.999896079820432	0.999945478298828	0.999966924257424	0.999958925474151	0.999903703220693	0.999968448069067	0.999833844729767	0.999942015734201	0.999895286400374	0.999886056093684	0	0.520069347880992	0.499428547453516	0.261402882809749	0.459621834958494	0.839649907947118	0.361750190763396	0.345007669182690	0.397561309984667	0.295207499582269	0.679930481608646	0.266340947861449	0.810190083098864	0.906995655989403	0.692274932214124	0.195901940375134	0.790889778273233	0.957426464045492	0.999237586579201	0.829507226166486	0.997080207429715	0.968058986843854	0.959221185803091	0.890457946135880	0.949550671617677	0.967366818659841	0.977517145624929	0.998077611438110	0.980596101846579	0.976635008752394	0.954297896303818	0.984389429903092	0.985644842178698	0.996763537441267	0.999966856809755	0.984906652079502	0.995396014517669	0.993854604985361	0.997190326480203	0.996652503856313	0.998356567844440	0.999573900864164	0.996974733276249	0.998964119048257	0.999048128963286	0.999935113034693	0.998777782114374	0.999574607841942	0.999589583905724	0.999914246526537	0.999355911900864	0.998967866808534	0.997462893769850	0.999913636358985	0.999623723916044	0.999996652723952	0.999372871923675	0.999757992417830	0.999582765384187	0.999679644370657	0.999781519135838	0.999842301604703	0.999568478808139	0.999916955735133	0.999608101182994	0.999514651327840	0.999976960013945	0.999749988869379	0.999274460625004	0.999734585727253	0.999662789757065	0.999811917006736	0.999865820999028	0.999962700366513	0.999462660352572	0.999921563659684	0.999812790646405	0.999783686249201	0.999953569930785	0.999735891017334	0.999951599230680	0.999983608374024	0.999838095434920	0.999798664886868	0.999579594305036	0.999834063314667	0.999878919595215	0.999796176591174	0.999968400636969	0.999908263804003	0.999935076771951	0.999802187210181	0.999917142696932	0.999989618895351	0.999989619612784	0.999926651508704	0.999983718237068	0.999885656172968	0.999891944297603	0.999918426493041	0.999996280876820	0.569436990543219	0	0.464284588255886	0	0.557354742565282	0.404688544292940	0.420084803872656	0.350324576428529	0.526178626080371	0.309022027018080	0.992107848772825	0.281959712946784	0.243679451910119	0.900792099270917	0.834648328264939	0.826550649345417	0.804994516496046	0.971770033087530	0.853966376349501	0.923254590294265	0.896394296843230	0.859755840095754	0.994272260235134	0.967538927151171	0.910494538832171	0.980381622631246	0.967594783894003	0.987130463305549	0.996087116613515	0.988501486525956	0.979077041543726	0.991912745152237	0.986031804904798	0.995911432451779	0.993827465530431	0.987486166827270	0.996162524796613	0.996748082837062	0.996436767816506	0.995298790331962	0.993149711953464	0.999351741896529	0.999134471244732	0.998621857225114	0.999307961601920	0.998417484806418	0.997045753865988	0.998864294306255	0.998636421902778	0.999759888762480	0.999157229050573	0.999107493541665	0.999797565789916	0.999521084956203	0.999752291829306	0.999766619729318	0.999459483089868	0.999468223583642	0.999684372378311	0.999667749310114	0.999970691517176	0.999892564504791	0.999938951139638	0.999754923915655	0.999810651776045	0.999473621136759	0.999892968881619	0.999976346349997	0.999905856032653	0.999998828043726	0.999840716455451	0.999961739397907	0.999879600006063	0.999746327404752	0.999985339228472	0.999974591192692	0.999868767775633	0.999867537812827	0.999765420331196	0.999887217251482	0.999805668024550	0.999930110916340	0.999926385019006	0.999632615373271	0.999948186254281	0.999820108686803	0.999891570118823	0.999969595136297	0.999964773350749	0.999871768115030	0.999898825164156	0.999912130631731	0.999840259056140	0.999831695248414	0.999886191154070	0.999919225559757	0.999979497808870	0.999931780601419	0.999918865080889	0.999942960824902	0.999825040168366	0.456898575848478	0.520987831741053	0.435757489321775	0.353390996071195	0.416087841143311	0	0.912946147483259	0.352492277407770	0.777353604155916	0.468119382515640	0.945031324699856	0.271856692103840	0.856120572323533	0.946627580239240	0.975083711216171	0.898715855440759	0.766198362436664	0.904037056013402	0.956730203250125	0.990240044387622	0.891958400716853	0.992968508692759	0.944480467520235	0.990088651466837	0.894052218852032	0.976386754243062	0.983048833926351	0.966044554648447	0.984805270242630	0.962606377034094	0.980351607168753	0.983540331040992	0.995747592786860	0.970594914098271	0.997697692908290	0.992672364485138	0.999478010981649	0.995508717763488	0.993676770051954	0.997439975069576	0.995361095925907	0.998898458965833	0.996424425072028	0.998773479082548	0.999477870745146	0.999572928760443	0.999439830719917	0.997649161697484	0.999619396850821	0.999987943032725	0.998666780702448	0.998351738537967	0.999810601107962	0.999491119504737	0.999468693436382	0.999762610910896	0.999513086330356	0.999911820375997	0.999974417163264	0.999951719930919	0.999628188389766	0.999935305708939	0.999458595477751	0.999661533388424	0.999936785998064	0.999812969365277	0.999956802226862	0.999761338308773	0.999753072841538	0.999894582675457	0.999730939676107	0.999985455248405	0.999527118360927	0.999510976948167	0.999630017109117	0.999971541860104	0.999921583995248	0.999938757749537	0.999923064565513	0.999766854205952	0.999681630716161	0.999979804408610	0.999864202745387	0.999924002869063	0.999993337079226	0.999964013548852	0.999793155793099	0.999936412530859	0.999963127707393	0.999966087278777	0.999850152562379	0.999785652123721	0.999965467371555	0.999984879172217	0.999925965080704	0.999869563027001	0.999865617735205	0.999959647177437	0.999983769429090	0.999884459616514	0.999856136127624	0	0	0.329874118091220	0.0890919573049490	0.434975238346806	0.385862453105337	0.916369697404365	0.355176465416476	0.995563074655234	0.300743674507536	0.500532225885067	0.926291730141523	0.981721471075372	0.235854000020603	0.232702999736887	0.869837750193799	0.832321658789387	0.957159616778594	0.861930441364443	0.831746409717783	0.949870577442897	0.977738685757022	0.948841531427355	0.939130312244344	0.987303059604195	0.944822468974262	0.976764806226111	0.996228615352088	0.966321291424651	0.992404626017238	0.999273453302615	0.996128327824930	0.994665711066771	0.987307969246135	0.987825966914311	0.990908320018947	0.991424537096134	0.997087503397010	0.997893472054629	0.998454990421600	0.996202407323694	0.995439569106198	0.997504438680787	0.998473237372809	0.998019218251681	0.999228367090286	0.999182774479031	0.999058270488288	0.999008786801762	0.998682451375365	0.999581844006831	0.998411361403701	0.998280305852504	0.999484275552773	0.999420599158657	0.999500283718358	0.999878128235235	0.998876923567169	0.999911265983796	0.999941103317059	0.999910200914590	0.999177811335567	0.999788271016981	0.999848828848653	0.999389328771790	0.999738981422421	0.999889984394337	0.999925125240357	0.999985203324987	0.999793097045637	0.999723262373948	0.999745725325578	0.999977516622482	0.999706681528292	0.999832849548891	0.999950131609879	0.999666380059571	0.999746521675386	0.999925583966719	0.999933796346466	0.999768193939689	0.999964505932148	0.999939816383077	0.999881537228858	0.999891472528121	0.999871998881145	0.999992485736959	0.999969114638977	0.999863777677454	0.999792430972098	0.999757319288686	0.999970784695631	0.999969981062835	0.999932617224467	0.999983513504758	0.999903951613986	0.999964150150268	0.999967643737969	0.999963016148396	0.999928215174053	0.999944396222244	0	0.498227319703443	0.467639169089977	0.462317894296819	0.424776621663988	0.900948216923771	0.389315698992207	0.928274856975809	0.338653848974380	0.950246500523791	0.971918574179112	0.566574240127595	0.457680166669971	0.538283108520248	0.540031100592112	0.939541366395949	0.939749553696929	0.877308606823208	0.933536574653740	0.778092919072674	0.898732179085036	0.888190607626404	0.979918800219909	0.990053549942720	0.991473423015999	0.945131570097466	0.973590807731569	0.974357954140151	0.990043681398347	0.997646983850391	0.986668550514871	0.990663039890797	0.994866384325756	0.989819892824962	0.990330462503710	0.995375876206466	0.998578950644675	0.993022111715054	0.995558442221393	0.999742453177362	0.996381455098362	0.997692036396744	0.998528204373566	0.998998945536038	0.999905682884788	0.999219229349978	0.999024252012844	0.997487217423613	0.998919085594189	0.999809310657101	0.999732403677708	0.999902336550072	0.998660692378284	0.999470811923546	0.999506991789953	0.998992023807659	0.999719079172132	0.999879865220508	0.999553613518473	0.999735977555194	0.999776288148671	0.999560892633507	0.999774374088541	0.999812552553234	0.999696038064195	0.999960602330713	0.999441846724806	0.999947923868002	0.999815384661079	0.999910432704633	0.999956290913012	0.999882187256214	0.999594965606756	0.999828409191079	0.999760815270317	0.999946995538335	0.999780921397409	0.999876646773895	0.999996722025050	0.999787485759146	0.999754623410886	0.999950087288465	0.999792164988007	0.999951626658591	0.999957388114388	0.999873315767028	0.999973016649940	0.999952383352495	0.999932405983433	0.999978380652512	0.999903510490713	0.999895042976531	0.999754543567824	0.999909154390845	0.999834175466834	0.999964094921084	0.999954361644715	0.999971032688549	0.999931770897333	0.999929948790404	0.999999068920527	0.551794801506410	0.532085775497319	0.0497750057696406	0.469219461173359	0.425131440094765	0.0877846089118225	0.360308156582551	0.832912116069112	0.820693325291603	0.312031514232951	0.272573660620417	0.368409658455060	0.252231323032175	0.899749329259040	0.901117661224088	0.712054091417667	0.809876424087535	0.456851952281865	0.925877750993325	0.917066291353966	0.871736464147194	0.891103147579928	0.965620024439007	0.895283482399671	0.956539128390890	0.986486072927136	0.982103679771109	0.972337895763835	0.972926290717799	0.967606098372266	0.995250823416456	0.980941173385565	0.993512703233379	0.994009724834259	0.990371392937645	0.993182932959091	0.986525436455425	0.995740073740093	0.998819876964662	0.990458404512855	0.997852434595649	0.995352211725397	0.996090080357739	0.999460553621627	0.998405532272728	0.999284933759004	0.998830039675715	0.999219290558480	0.999985024933911	0.998873500586753	0.999064564931896	0.999066772210156	0.998909997537345	0.999932224439583	0.999122004023812	0.999954669439448	0.999680951190309	0.999627612953251	0.999866238812581	0.999783636343076	0.999980936703768	0.999895594490077	0.999758782699052	0.999964512698320	0.999980815540630	0.999783477191037	0.999670497137417	0.999749375463310	0.999649514130917	0.999913625243521	0.999986697833072	0.999857494712400	0.999884849296624	0.999994678456247	0.999969615624605	0.999950144625087	0.999931420508135	0.999998404466298	0.999958812742032	0.999922090346412	0.999682473214546	0.999988814964534	0.999941473908417	0.999928978150881	0.999850527523836	0.999865314840773	0.999926435351293	0.999889866283652	0.999974791136944	0.999877565455119	0.999998727345192	0.999949402710884	0.999845385092119	0.999814372236139	0.999958928808298	0.999983306316711	0.999904178202439	0.999919379519402	0.999926534220615	0.999962704495685	0.999849034479988	0.569488298467048	0.514386275598912	0.472527627107758	0.227354221763296	0.425569227242944	0.993397196155450	0.360469191496017	0.626332272768493	0.609682850637892	0.773348567959799	0.942096099382023	0.268568764812450	0.478207932358518	0.722248355362557	0.872077525974156	0.927122800695662	0.894506974554857	0.673191698886305	0.886496408361112	0.788047492958177	0.955776131372353	0.892599040327574	0.890881144398461	0.941608472996949	0.999473093853955	0.949880692880793	0.971469114456478	0.997181399277274	0.974413108670953	0.981024166947716	0.990325308100586	0.999419982068877	0.979881083928008	0.998305271046168	0.997859103038447	0.998655180464638	0.995497946859369	0.992096834564358	0.998905783965143	0.998887951500897	0.997897943146034	0.995221486018668	0.998545770802139	0.999246418751596	0.998628930841127	0.998997464192198	0.997138691948151	0.997970046512098	0.998712934524283	0.999956583889504	0.999262928868387	0.999925761749706	0.999089186910618	0.999895409407221	0.999415814690347	0.999668209989111	0.999974763274884	0.998901502012981	0.999657187451830	0.999904014902738	0.999553077696174	0.999802331608985	0.999631098401919	0.999999444991699	0.999727407105260	0.999716293896139	0.999674497956150	0.999878924222784	0.999849789994101	0.999969105340167	0.999780587222999	0.999958741339541	0.999918281570578	0.999750349326250	0.999731929325157	0.999869760487394	0.999931403914268	0.999924140726186	0.999949806908053	0.999753516746235	0.999895824994188	0.999994963777281	0.999896178744387	0.999860398487442	0.999739363450629	0.999883058122604	0.999874971805284	0.999773353949034	0.999980203459304	0.999854310234604	0.999998021648477	0.999976498499304	0.999918658797744	0.999924753807083	0.999881793108533	0.999909208327682	0.999953646371311	0.999951869694280	0.999971055349496	0.999975054949319	0.999937925081744	0.906297896678967	0	0.656147589478978	0.186435185151412	0.415004976766125	0.242686855070613	0.381195608702452	0.853929871219318	0.961684305585090	0.670003193407380	0.518168310588118	0.390967346891712	0.858190207120262	0.396821882080931	0.746488979925237	0.948131131978838	0.969567217691047	0.983431144670406	0.989239571273310	0.983621266736553	0.977551322614005	0.967219577978618	0.987825718100210	0.994325478039657	0.986778942356086	0.950952349221569	0.949753590958449	0.986633881743699	0.965118927875150	0.999235098318272	0.967997232334623	0.986916125544675	0.999653328147295	0.997431738238543	0.998168499030759	0.994121958022808	0.998676480922502	0.997643156673488	0.996016916445769	0.999252643249547	0.993033505945930	0.999932257977251	0.998580171297900	0.999159116523223	0.999493603726398	0.998188577469230	0.999391579951821	0.999022986404598	0.998980010591890	0.999252213519874	0.999063655245393	0.998584793976207	0.999533949468715	0.999341390779425	0.999835558015240	0.999769493450577	0.999664573342856	0.999675609236424	0.999898601471577	0.999567982062679	0.999451133554478	0.999684505171155	0.999739827939096	0.999834911410276	0.999400914886860	0.999928441223910	0.999844571770779	0.999908096049425	0.999727569549174	0.999889540412739	0.999868352795376	0.999742701639393	0.999924727690187	0.999745826552700	0.999979604035955	0.999969059475094	0.999844805433655	0.999732692123908	0.999848398659680	0.999936839813969	0.999744505615605	0.999520183679158	0.999960947669852	0.999956872926740	0.999745296399514	0.999909079940879	0.999895114214950	0.999934856329380	0.999887298069693	0.999929750421869	0.999928407017967	0.999949542635215	0.999988080098111	0.999866449130068	0.999960167549831	0.999893691081046	0.999840332536025	0.999969429001118	0.999795127855848	0.999909242742130	0.999992113298140	0.563592777510090	0	0	0.456924972538777	0.577021719907488	0.317601119408422	0.119965014657547	0.338644630987298	0.630545408122175	0.510835445603145	0.631989453647473	0.766849383590663	0.948481610309541	0.963514988356080	0.550691264219324	0.804553881275639	0.717997026943696	0.957524095452449	0.734908041236408	0.860595572489238	0.970974647931131	0.967741623417168	0.922624221036405	0.988282753035268	0.992413303428312	0.924138762334805	0.985990834647063	0.991565339259040	0.997901565384944	0.993686719476562	0.984455360704016	0.983333607955603	0.996397813228954	0.985917400104211	0.996394054715789	0.998160927100192	0.997455634180029	0.988534927545459	0.999608344137134	0.992352762802927	0.995916050370195	0.999445998271061	0.997805590859837	0.999666007905279	0.999926494076645	0.997513116026033	0.998996296013956	0.999801905983008	0.999935920170182	0.999776549908764	0.999676704471110	0.998950612406468	0.998952234796446	0.999481476245207	0.999482459136104	0.998749847893703	0.999936329813318	0.999967982575657	0.999516887131349	0.999196370389822	0.999862146679076	0.999565602600489	0.999732093659266	0.999627700101363	0.999804219450090	0.999954380866836	0.999934051013954	0.999850774335579	0.999892696820646	0.999888764327801	0.999990509201954	0.999871197879198	0.999987669891166	0.999607266163229	0.999970718999322	0.999857754043509	0.999872922869494	0.999964030662049	0.999978929334776	0.999968088461208	0.999920555884966	0.999999572933855	0.999880447045778	0.999563328325028	0.999926598975410	0.999722727996021	0.999770319509593	0.999913517207791	0.999832564567427	0.999967610881194	0.999922491893064	0.999943655927293	0.999989796377544	0.999936301711477	0.999996475030086	0.999911124423919	0.999907440931931	0.999964663220484	0.999917773709889	0.999942493033467	0.999921412862923	0.835839882431067	0.349790548948655	0.819517399562478	0.470962105894371	0.420659467757677	0.465193177754108	0.140984784217619	0.355999619132381	0.310857380237914	0.313640078706481	0.280130340873983	0.627614183244433	0.247684692481599	0.700910368280235	0.924867851190594	0.205658053002497	0.882848292158854	0.644507698607515	0.992093137810892	0.773768275540012	0.988740880583756	0.863711941469045	0.971968899084001	0.970259067965734	0.984044291056299	0.997185975411450	0.994606724729416	0.992368232311404	0.985239918202881	0.976647331045577	0.975609201797870	0.994980512146529	0.992815243707323	0.989562080240747	0.987506794558314	0.996334473458197	0.991876349389173	0.996044795087957	0.994999305278516	0.995377490620575	0.998699723567382	0.998838217549681	0.999625614878741	0.997793568912307	0.999114422131418	0.998139564041528	0.999755728754381	0.999471879655945	0.998646834648662	0.998547432530859	0.999623309169149	0.999409234814080	0.999055972716030	0.999200545865660	0.999878545648935	0.999780804334521	0.999404648394750	0.999745052439410	0.999818491971244	0.999015754473818	0.998998583907745	0.999440213737316	0.999893781224726	0.999824821845966	0.999663016425314	0.999721638655890	0.999740300330651	0.999991252382940	0.999806404452415	0.999918767092116	0.999844514713706	0.999977867542387	0.999980947308866	0.999874861827006	0.999553321951542	0.999960784603351	0.999955367116077	0.999789166750500	0.999788123929978	0.999934315238415	0.999983859403770	0.999907565315084	0.999925578515783	0.999909849452689	0.999758664461375	0.999954715936984	0.999911319871032	0.999798184156276	0.999987922442621	0.999876824348935	0.999839381485465	0.999912477554439	0.999846814910527	0.999974847923012	0.999976049536843	0.999827979777633	0.999910927890090	0.999892350667757	0.999969307976603	0.999959756329630	0.999863470374265	0.703994072605609	0	0.497627288949124	0.461051142190926	0.127891544720986	0.330276885655552	0.378556312727359	0.730815038996364	0.325856758655662	0.999448841132355	0.705938697974270	0.261568797115605	0.952980234930210	0.657571186533285	0.948127339274819	0.202951167138815	0.912822592358468	0.813949109586142	0.996937084863716	0.991285638252429	0.933142812509964	0.988734363884399	0.872699313106916	0.954862025950676	0.982102832271310	0.940571342325010	0.963313500099681	0.999390279171583	0.973781555339403	0.989151191716292	0.995178909857575	0.993390982128946	0.984025738380311	0.997726901029501	0.997132206031957	0.999639647875752	0.992532741151292	0.998633964514220	0.997497797886613	0.996472625666422	0.996803844554126	0.997980399057172	0.997067858633714	0.998113709172944	0.999188366040523	0.998586060130674	0.998589587355891	0.999162050644849	0.999893834509112	0.999140538459944	0.998975717875955	0.999261743641187	0.999715791478939	0.999928619901496	0.999885668060735	0.999649092032666	0.999374006273337	0.999858860622985	0.999757316852684	0.999767405178384	0.999608236142024	0.999994435198197	0.999566975604912	0.998953474614955	0.999900320311311	0.999529442225849	0.999735326807243	0.999830783251995	0.999420849364548	0.999667310758541	0.999728056485690	0.999778386524679	0.999801097678583	0.999803725844311	0.999821656189441	0.999936820618081	0.999900798180250	0.999684401873046	0.999912599604733	0.999944852315664	0.999734161816266	0.999791107809742	0.999956851730916	0.999995174670831	0.999839814271086	0.999945639819132	0.999833362010413	0.999903840865609	0.999873307832228	0.999947694141175	0.999930501342421	0.999967384445400	0.999941336925417	0.999920270912512	0.999891556916731	0.999944509193059	0.999981901055414	0.999860863386984	0.999997223527806	0.999987440029606	0.999919548742272	0.577670685655482	0.667255010735355	0	0.460960779188750	0	0.397197739819385	0.368242048750157	0.354858190699117	0.962678184289567	0.567994757408314	0.914037302151285	0.789211467221016	0.247810217911355	0.228649067612894	0.985923289799691	0.609774337979846	0.908941184128500	0.981970325791747	0.908934423295845	0.846135871316903	0.928301985449298	0.900970087435091	0.992311836861022	0.987842367630447	0.971979090495387	0.995728941610856	0.984770570004221	0.987968963791382	0.978801490097747	0.963266523462891	0.995707636210001	0.990271852907666	0.982660487788670	0.994844563864095	0.989251444814424	0.997032990731690	0.992950020161280	0.996995941029671	0.993881686271366	0.998959296005045	0.999809531198340	0.995787716933275	0.996261851897724	0.998564326832450	0.996593705244477	0.998835738696196	0.999989987344040	0.997256241703625	0.998936774751431	0.998885336984589	0.999258035386104	0.999354453477193	0.999621600947571	0.999912037396802	0.999862702048673	0.999574161578681	0.999815758274668	0.999193418069209	0.999857996073140	0.999834097849402	0.999567151978678	0.999813933218139	0.999756848874726	0.999849150028486	0.999873778579155	0.999844584373848	0.999961163164929	0.999469661986270	0.999856666204576	0.999305590102091	0.999895954501620	0.999874774084149	0.999976770508838	0.999654255809753	0.999581237123527	0.999876581316723	0.999948959488715	0.999948048582099	0.999910739493585	0.999860361088110	0.999700100632945	0.999849698260305	0.999792303182682	0.999896453687190	0.999990907843058	0.999880391076826	0.999912831341525	0.999941343667384	0.999878536545638	0.999941381631021	0.999960679722217	0.999961213876632	0.999922086657153	0.999840745738407	0.999833429553196	0.999892215364976	0.999886679555702	0.999983576018890	0.999916337576683	0.999980708956883	0.999978061332500	0.00294769054347233	0.524329231519417	0.511322143911243	0.179883764994593	0.435426810350205	0.857757686340626	0.336110553913647	0.541175726489277	0.322074804032125	0.699300631747015	0.774885132827685	0.823292353979584	0.993181992207714	0.767718853477912	0.662218644268529	0.988435709916014	0.956677051507403	0.818641563609885	0.848175590965358	0.950383861499932	0.901582612355236	0.976085278904274	0.928745776794082	0.913949015779074	0.965501267091565	0.919824684558261	0.976707135148570	0.953617112717039	0.982655759299837	0.980247306651450	0.998228124045898	0.999085768245954	0.986671782402053	0.997271165525561	0.988676451121098	0.999743724308846	0.999491806876001	0.994386613504622	0.998227001080888	0.996830954591386	0.998057816934631	0.998100839932903	0.997947938600883	0.996822221456497	0.997242144421662	0.998998562626579	0.999096950750405	0.999228023624840	0.998763443423943	0.999410999679440	0.999774067994210	0.999140845391842	0.999878814973080	0.999979916238536	0.999791101443001	0.999783705933235	0.999587324146871	0.999898961202602	0.999428265860500	0.999630176856865	0.999650329886856	0.999791033403892	0.999860439975808	0.999946613698352	0.999783747754644	0.999924982608256	0.999987507383080	0.999827871801873	0.999699986476916	0.999684517597924	0.999879923023841	0.999895627423989	0.999970922233688	0.999881451571310	0.999801733739019	0.999839029994907	0.999988988954804	0.999774759980868	0.999809566921356	0.999969115152885	0.999884064205724	0.999977074827032	0.999830906074206	0.999909878503779	0.999989893737030	0.999975621305788	0.999931828052545	0.999866953892671	0.999958355914788	0.999995684655464	0.999892308535373	0.999886677022243	0.999970996429565	0.999986560650499	0.999999740352805	0.999985736766838	0.999867289605068	0.999873851975216	0.999978800365349	0.999900967162762	0.999954294937423	0.0525956379729740	0	0.474884270175039	0	0.561325609279249	0.379620454265786	0.499058985350188	0.438313151927112	0.293563523675686	0.359591445620876	0.285326575638772	0.711155123002806	0.873830232603723	0.910374300415032	0.779996997749719	0.897144256664676	0.757119029950955	0.897590859514827	0.938717794043334	0.949364958604974	0.883709906730829	0.947927668215317	0.988519571770975	0.949891225681809	0.953032093052555	0.858784533146590	0.999950033727318	0.989589085197711	0.980968023144446	0.992889357917638	0.966182116113627	0.990112819047956	0.999633627703186	0.980225766652473	0.996219289559357	0.995754728060167	0.997223681264244	0.997165869496323	0.993852122712451	0.997293491900021	0.998300668089435	0.997450247618796	0.998350362679352	0.996125724920711	0.998440565328263	0.998846207633941	0.997804435487371	0.999920716022283	0.998401823388643	0.999049364459193	0.998561334517046	0.999609051139286	0.999863305339333	0.999924686074411	0.999588449826789	0.999243383057197	0.999840803585443	0.999696640713219	0.999362726625653	0.999850086682229	0.999717332984911	0.999578181970287	0.999730218936145	0.999522097516468	0.999789247607293	0.999830977262589	0.999592923224448	0.999979622090688	0.999945210263851	0.999819914531862	0.999914066285389	0.999868461815302	0.999867287648768	0.999840346297418	0.999950261094719	0.999963321669315	0.999859826476434	0.999990799352925	0.999746070977070	0.999998533023086	0.999686996475373	0.999862058295129	0.999832303843189	0.999708293005754	0.999809667804239	0.999936948864957	0.999913778712219	0.999957288855573	0.999976114054153	0.999929182708370	0.999850850987941	0.999976626430725	0.999856927822871	0.999977480405124	0.999980614354994	0.999885641560325	0.999989240035378	0.999971433741952	0.999932046556606	0.999971304853996	0.999870870287569	0.571513109602240	0	0.487155816639034	0.453211823648758	0.709377134447393	0.386555539174566	0.270327023659854	0.354142187321563	0.359828294385771	0.775157375214454	0.455020176489314	0.921674798643473	0.848646856701333	0.888830512331398	0.845583706135315	0.641129163587415	0.877323100432538	0.949857685773401	0.742067894249568	0.977420715470005	0.883057288610376	0.956585706077595	0.971932633868712	0.961289626652366	0.982394931932782	0.964128589451272	0.995368293972099	0.951524422128653	0.993372282727304	0.978004255636873	0.993241843642188	0.992995381569111	0.993033488644606	0.992592487465519	0.990881306627053	0.987824304439321	0.999033232942596	0.996671752787852	0.994325398419016	0.996776812725664	0.998492533980255	0.995776833507532	0.998524769604498	0.999021271985387	0.999560224112330	0.998913283344058	0.998943622909916	0.999947427102983	0.998262101684585	0.999214506589977	0.998957012774672	0.998875127116384	0.999882128829384	0.998873109890166	0.999838838800147	0.999298121086990	0.999512306710385	0.999597181336017	0.999840778527753	0.999840675270358	0.999671240743376	0.999812859983183	0.999445340587584	0.999663370173177	0.999399514965564	0.999983462355233	0.999591101003736	0.999765291793778	0.999528531024246	0.999871518888027	0.999899392617364	0.999996538444658	0.999765659988759	0.999974513503779	0.999965577797655	0.999833379027232	0.999991900730409	0.999894158032387	0.999924175372247	0.999846402369889	0.999948387751606	0.999949445635873	0.999894413881950	0.999932545307791	0.999999883965845	0.999977297675633	0.999995650149923	0.999736304137938	0.999985133457777	0.999974818760233	0.999846846144562	0.999939580300242	0.999951992176597	0.999940057232591	0.999897244775785	0.999908309916423	0.999950142054718	0.999923849426275	0.999973911453575	0.999971620831998	0.999982979721400	0.586787996737601	0.522916829318605	0.494005073867201	0.466683542383180	0.407519364134543	0.395718535015951	0.356091088375129	0.353716991361657	0.329354822696480	0.695513684987934	0.290766920211863	0.826238373848945	0.806721227891872	0.710508585965385	0.958515931689409	0.597539734289031	0.753651354853548	0.889160957709451	0.891345818939875	0.975593128079873	0.997926695367598	0.860628196330044	0.952577378947538	0.999297116933127	0.955826554773037	0.960896067739931	0.968872836574168	0.975673508458927	0.986402083317918	0.981998686049432	0.961530429231304	0.999918863791220	0.993195542693380	0.991755916781982	0.998699357425571	0.999724774609959	0.998402324466179	0.998775374340792	0.994518943896655	0.997370811700054	0.999670105137941	0.996273221874427	0.999338046681528	0.999557222313002	0.996594930337071	0.999300161793924	0.999599023515970	0.999876135438561	0.998368226704940	0.999820396910803	0.998387253382863	0.999835394094302	0.999878033274978	0.999047056953913	0.999601667871809	0.999808252232658	0.999558498577534	0.999782236089679	0.999680318390982	0.999647021895494	0.999932059039210	0.999361173520898	0.999966971421216	0.999644320492562	0.999791087355353	0.999909757012958	0.999770932674497	0.999743235041625	0.999942198135395	0.999854628373261	0.999864161304890	0.999958557857738	0.999879618265177	0.999891242898897	0.999851590780219	0.999786447063455	0.999949153556817	0.999960985911223	0.999811588457877	0.999826999815569	0.999731086416271	0.999857847133877	0.999823559431053	0.999812959632416	0.999986696708466	0.999949197498707	0.999960696385125	0.999932199961958	0.999976295623093	0.999926110239547	0.999993391518449	0.999895498871465	0.999846059352684	0.999980577390060	0.999924889414920	0.999988407979968	0.999927574567857	0.999942119178285	0.999912986682071	0.999936607512153	0.999921131426416	0.559104044728656	0.521318495823079	0.176899744080743	0.557071035070616	0.122283988979898	0.411480393496462	0.187421188563801	0.148140389892893	0.929723749996852	0.298467982276157	0.619185439251303	0.875916762003416	0.543495317510082	0.666129378865666	0.845650011575617	0.806196767856705	0.985397603252255	0.834587335981160	0.995603621876087	0.835969118666416	0.898364289551022	0.947298408829532	0.972485156693921	0.990976817979483	0.942817435095183	0.977606915381362	0.989315523741169	0.996592307246579	0.983709242949451	0.993300089994169	0.987544174665617	0.993810192041224	0.991526924890912	0.994515410516972	0.999807888721388	0.998990214370081	0.996520045237663	0.999885217268494	0.993272959820878	0.999549369218236	0.997986178653243	0.996246578369243	0.996469899048592	0.998467224434313	0.998653154110829	0.998315064147068	0.999088998171470	0.998958303073095	0.999616718329071	0.999001370018474	0.999875203135243	0.999848517336821	0.999205596576922	0.999164936373919	0.999640953149111	0.999836258987019	0.999636992311740	0.999713456992195	0.999718552261247	0.999912408309276	0.999253811389905	0.998727773104614	0.999980801498326	0.999757909067161	0.999231198908854	0.999551731872423	0.999773340499276	0.999635684045606	0.999799393685392	0.999889739404063	0.999891495735129	0.999787914288533	0.999994603982422	0.999806143431835	0.999994829185531	0.999854751730372	0.999806752195952	0.999975564151970	0.999884606235820	0.999988883208273	0.999971917925610	0.999894542783368	0.999958394292401	0.999965044723263	0.999901154326361	0.999918290126272	0.999795303907298	0.999939081052533	0.999929050190059	0.999942332576558	0.999978311796853	0.999910927873917	0.999971617754773	0.999901178521047	0.999967340397194	0.999977591391253	0.999978523013615	0.999952081404799	0.999797355021152	0.999947093758545	0.999908551876478	0.577128473332780	0	0.487508744928062	0.537976897599891	0.507465211807589	0.620082603335949	0.794496996647898	0.361490929876504	0.331513672722549	0.312751131933760	0.562343295102995	0.429197676621656	0.687781238727574	0.932271643583515	0.913062317846913	0.892775089809249	0.673984240972883	0.855273805342561	0.988839243891979	0.906697866620421	0.936352024981334	0.892572020383072	0.994493096043567	0.954448067535971	0.980654156351553	0.927372561693173	0.945741984658542	0.998713204946470	0.991792289244424	0.992372803184911	0.982978438376482	0.986837715541928	0.987847809510915	0.995120737278382	0.995630392302033	0.993762737756685	0.986486081015115	0.995351255251808	0.992497182474444	0.999380515236717	0.998438219809846	0.998083149855913	0.999050137113934	0.998207631593557	0.996816341797892	0.997369310928316	0.999933593699758	0.997871691711840	0.998644484116947	0.999891697174613	0.999436941064933	0.999531235683909	0.999717050012620	0.999046205751945	0.999895515101430	0.999700842851836	0.999567757304936	0.999644786511999	0.999634187833751	0.999843176843395	0.999314601103919	0.999885107546117	0.999851767108850	0.999510701601359	0.999701057384917	0.999709050584977	0.999787468085637	0.999701373235429	0.999820453820061	0.999940531644884	0.999753455842655	0.999908999798010	0.999973937018383	0.999972354038397	0.999947877180726	0.999940933346621	0.999961987873795	0.999850187770847	0.999785804005991	0.999952490191122	0.999780516749741	0.999991311942389	0.999952744422349	0.999946718841928	0.999940640742165	0.999877829050666	0.999938098789146	0.999724659576551	0.999984349183146	0.999975357584841	0.999985349155702	0.999939176173010	0.999932704740273	0.999974708443813	0.999976325348777	0.999938132072929	0.999995621867336	0.999960498989085	0.999901495045384	0.999959551471361	0.999934016213153	0.556972451096057	0.526709728385178	0.792582036954925	0	0.681439885021303	0.405461270337273	0.589599727382914	0.354671748691500	0.336475966949961	0.547964964918938	0.273387303261316	0.582381416233904	0.745839527009050	0.773550883997443	0.922993673412006	0.776402495522393	0.952556818582193	0.848954888193272	0.869730732694364	0.916294254014395	0.864072067809536	0.998288004890665	0.858304695010147	0.983590484308170	0.982607663121322	0.924747098866691	0.938281067624299	0.981623784302593	0.996390503244835	0.995970832798767	0.995159895056140	0.983754258298534	0.988480720101022	0.994621899848977	0.995493981741161	0.998804557152397	0.997162313877608	0.994382081455194	0.997672854514450	0.998177538062178	0.996496472199537	0.998225551305613	0.998714926529499	0.998166470339354	0.998584449793324	0.997594554936366	0.998229211612783	0.999078161521126	0.997179522553671	0.999116549830289	0.999565502699530	0.999761208363478	0.999239089529832	0.999489261207137	0.999491834666129	0.999694656858859	0.999811299969632	0.999527617798661	0.999594286734509	0.999739186410391	0.999786810940356	0.999992725329283	0.999886081163435	0.999878087995250	0.999800911024655	0.999940587462091	0.999994736661912	0.999677619949001	0.999768413744497	0.999836953795947	0.999727845211259	0.999751444262638	0.999874869137119	0.999783389110916	0.999852030080940	0.999871869986629	0.999924870838023	0.999896827407966	0.999919971015852	0.999744574267057	0.999849510676032	0.999590999236021	0.999938986381176	0.999905613960435	0.999887310991431	0.999873610900221	0.999847285462659	0.999818687423886	0.999946702633952	0.999859178360198	0.999967430020199	0.999916355335981	0.999951435037708	0.999910711091993	0.999902568049512	0.999914244147351	0.999882189326459	0.999949217099540	0.999934890079835	0.999984623046976	0.999974846799402	0.0966134395561119	0.548484885359255	0.749031631167888	0.633213261935823	0.430237933938562	0.383094620086783	0.357721439387069	0.399875321501831	0.880572229358301	0.551252459824294	0.697660378798434	0.706664049005796	0.689925813296730	0.238273994499215	0.983562388523524	0.796634829622311	0.903766872330988	0.777694558883132	0.928743049996819	0.999210511848655	0.909000138792706	0.937890676591698	0.991884440323356	0.958802453602626	0.990166623896899	0.969271005653984	0.939653753576297	0.972182571448108	0.965874873887038	0.982025806086683	0.998276877091827	0.985478673452133	0.998868827115034	0.996803117811890	0.979010693046038	0.995095817316451	0.998614109170608	0.995290802021213	0.999933419835941	0.998901251752526	0.998631212698458	0.997348612818783	0.995632131804447	0.999045811783139	0.999088853677666	0.999117491260943	0.999235077278023	0.998952347421966	0.999789572650275	0.999518017814747	0.999907480133537	0.999523825514050	0.999623439415321	0.999745875662187	0.999399790084153	0.999532516495024	0.999270657827958	0.999342777400788	0.999801536173099	0.999758418561688	0.999636338188864	0.999712001859456	0.999834815018823	0.999472977800626	0.999508532171698	0.999990218319901	0.999967821514612	0.999738349701449	0.999969944888190	0.999988318896210	0.999849076513874	0.999992507778323	0.999756785946927	0.999535692400617	0.999955358633900	0.999900829830721	0.999724435448397	0.999741835947999	0.999810841823232	0.999984194662439	0.999953282686948	0.999995947104395	0.999938892027820	0.999838864053193	0.999754661624413	0.999947277160296	0.999894618021844	0.999911203474407	0.999961545632504	0.999836091248170	0.999939970473463	0.999819918254378	0.999922759294924	0.999787537411830	0.999886972814447	0.999942052832598	0.999971504890561	0.999850252091184	0.999882586368822	0.999971218860808	0.999943703809539	0.574549708598977	0.713258502531822	0.972555799869854	0.917479313431007	0.296501390229630	0.170163855205443	0.800855911783059	0.529859612284821	0.688700397791062	0.684784752548074	0.286744701769421	0.998276788415879	0.841538690182916	0.695494488142737	0.864483747534160	0.885218720933816	0.664647860074660	0.821796350627389	0.892619400847219	0.975137715249173	0.887520371943529	0.947305112133239	0.905224002751053	0.966651507656990	0.996527021367163	0.985116326680101	0.985047594125751	0.981949615881805	0.977443703264111	0.984046741341163	0.999938130334983	0.995646823928523	0.993655780084590	0.995475896160832	0.999240658631931	0.995030950768378	0.997057400079668	0.999231235676129	0.999185365770034	0.992752664174932	0.996306732396145	0.998394463067881	0.997954489575189	0.997487321380168	0.999079724951913	0.998172153680793	0.997927247745751	0.998581661415694	0.998874491495453	0.999546556998076	0.999539622149504	0.999943764680876	0.999536043724026	0.999717372307815	0.999466384992919	0.999731953643245	0.999735435056320	0.999861161607891	0.999545194334380	0.999893210055831	0.999689404170353	0.999921035337430	0.999923839537974	0.999815442008886	0.999466939543085	0.999764670427947	0.999691561405889	0.999840069075908	0.999766819075391	0.999987927764864	0.999883075959007	0.999990463359144	0.999882612878498	0.999774191149945	0.999689561517991	0.999914719422177	0.999904078761935	0.999950664488654	0.999889235296059	0.999899172622904	0.999787988817495	0.999875105892946	0.999724041547466	0.999925311521308	0.999852107314730	0.999920673378750	0.999975518386581	0.999950785420661	0.999990680208826	0.999928217311630	0.999883201559958	0.999972830962463	0.999901925919997	0.999910362465058	0.999996854570025	0.999966268381940	0.999926527953194	0.999903569276268	0.999943683197694	0.999897774088784	0.999987313530750	0.571523092554643	0	0.305573766378031	0	0.438059647462076	0.582010294705415	0.805041543007053	0.473588296364629	0.468957096911239	0.704374573017589	0.905563020544050	0.900726100087744	0.942558813354952	0.735457157314193	0.860791471456774	0.696132815403077	0.834958985502722	0.936041869077028	0.979765822874238	0.954733194099049	0.949315878566395	0.903544915142461	0.917164521447104	0.953528008211101	0.981791334406870	0.984858157356044	0.990657203766778	0.948441787761600	0.989779581295247	0.984100058714828	0.997150524194240	0.994706185021274	0.980578053402112	0.999206050552688	0.994160063662289	0.997840516900434	0.996706225778516	0.992347619470997	0.992711631274724	0.998877879002588	0.998106421795075	0.999381794924423	0.999152989657486	0.998759557885235	0.998880718785948	0.999176300948519	0.999415699063775	0.999685298620802	0.998813541429243	0.999773563446877	0.999569732856448	0.999199764629044	0.999081949473447	0.999347089015844	0.999725327467036	0.998599788930521	0.999387078868189	0.998923290142977	0.999921548976876	0.999612561583360	0.999908164066989	0.999635986856543	0.999885801418021	0.998793869386421	0.999808829223500	0.999903341210702	0.999890175313357	0.999934095116543	0.999966855981596	0.999984262832822	0.999981780150225	0.999914133956860	0.999840173393631	0.999922262693939	0.999819894722430	0.999878819722738	0.999941988760751	0.999996654687005	0.999983846493815	0.999867579868283	0.999978142466631	0.999949360699659	0.999968196235421	0.999855385090939	0.999967057425769	0.999861806003587	0.999923785800809	0.999852998784095	0.999861193464693	0.999948080857478	0.999885917527140	0.999947089058100	0.999972277221152	0.999944014227977	0.999953133248792	0.999953508868829	0.999929777196764	0.999946648481717	0.999929684490896	0.999799228929405	0.999944450798663	0.677139033750207	0	0.730925290611509	0.518920547697537	0.866496269308743	0.671215524933870	0.778240045353694	0.426478725701342	0.915754455011400	0.595705432745725	0.971974651807748	0.923370862758431	0.570133075130824	0.231678400116472	0.895909991056321	0.970158831462485	0.877504871907414	0.842178458480657	0.943248784187446	0.990353649074910	0.912175858845631	0.906391219081400	0.969808056519328	0.939310519418299	0.938912176471149	0.981457893937238	0.942860751830879	0.998183829453062	0.970529431364996	0.994703820172200	0.977168406050896	0.999022142167146	0.999974772561692	0.981728559919516	0.987848618738592	0.998336239555425	0.999088717274274	0.991434564191793	0.999831493422877	0.998481432795569	0.999968849237769	0.995611987654767	0.999003033466089	0.997376105667873	0.997721041080576	0.999524680086693	0.999938686469855	0.999647343750713	0.999394750439125	0.999357880259802	0.997703076150647	0.999383291285605	0.999780484475143	0.999493919195990	0.999312002413483	0.999385160790621	0.999383375708982	0.999488213000960	0.999844534463525	0.999974095005462	0.999984361811413	0.999809198690553	0.999806443673793	0.999771922334784	0.999871524423006	0.999717860832122	0.999820794331294	0.999948082937649	0.999985648124084	0.999821229562082	0.999737073536925	0.999957181395114	0.999978288286973	0.999776688333524	0.999985114613628	0.999907351696773	0.999939331635770	0.999854124418646	0.999969862506009	0.999924898683695	0.999892169633564	0.999999609496778	0.999940061349361	0.999894566872181	0.999877709033808	0.999914208898864	0.999943085464256	0.999961126149107	0.999950546238531	0.999922667339073	0.999987516354913	0.999981932387179	0.999912903745452	0.999970267469797	0.999987651151661	0.999997774277552	0.999843424676869	0.999960920334313	0.999983208638381	0.999909704792159	0.999996192464765	0.611881756493355	0.533213658493255	0.564320651166892	0.475398562138118	0.735106085414443	0.807182691421747	0.872893565925325	0.0435246820618300	0.599586694990891	0.660814702146968	0.578752381680566	0.809389115960818	0.928440149127608	0.779357047705076	0.774472761370445	0.802992856605046	0.877248352594155	0.775755914258436	0.975183637834959	0.902778421821629	0.941123807344212	0.978062951115239	0.972648986080859	0.921650215393715	0.986138445111297	0.937396061608421	0.980888755242755	0.983118344840879	0.984959127730619	0.990140245880858	0.994885663854986	0.993252961234984	0.982184032664948	0.995759552013774	0.987154100110446	0.996638243611427	0.995975404736913	0.993300278574387	0.999316389872439	0.999324239793119	0.999965246489907	0.998435578117336	0.999317758186037	0.996964502217397	0.999669511400539	0.999754006468923	0.999734222924470	0.999972314185242	0.999024660212357	0.998948189902951	0.999345905788458	0.999364668485746	0.999408498791404	0.999944783243427	0.999648915800037	0.999390742013055	0.999377801217823	0.999676046544405	0.999710356762029	0.999490739283229	0.999936827395752	0.999580172775843	0.999750284305633	0.999735171342502	0.999563773577898	0.999694955330914	0.999822702783089	0.999757976114645	0.999827778022298	0.999949708298896	0.999938591720880	0.999804229677912	0.999800807951333	0.999903113663132	0.999770029162542	0.999865294103491	0.999933532243085	0.999866737292193	0.999908365650565	0.999897290632188	0.999934316621359	0.999878923028122	0.999933032150205	0.999975295184948	0.999986286660532	0.999944585455651	0.999899335130476	0.999920486570777	0.999873943438812	0.999937875921261	0.999999097565890	0.999813779122736	0.999929453385538	0.999942233218175	0.999994204986440	0.999990361960835	0.999920404074978	0.999997763087095	0.999969745181764	0.999963021733481	0.999890518427255	0.426394940940817	0.655084225341702	0.0584089686237710	0.179397711592571	0.0113071568785867	0.405941535316193	0.943553946202803	0.484836433330854	0.833727232166055	0.863954722781615	0.285980639913905	0.679863773515038	0.941068538243550	0.742221151171380	0.639178610951232	0.922819909591178	0.994529222649001	0.967798898629562	0.912792082689564	0.971531476548486	0.857150153107382	0.940922693965484	0.958577382192546	0.979909981513535	0.972314991329898	0.993727496458647	0.976906571425928	0.973507094740259	0.982316112739249	0.995822875464956	0.983366629749697	0.991088117948893	0.998377971551384	0.998513411856436	0.999434405053913	0.989666828520319	0.994650575769383	0.997680048059039	0.998122918396987	0.999555115657464	0.999344695120236	0.997351847687887	0.998416103203012	0.997808284774113	0.999135032795669	0.998789022211212	0.997585788739559	0.999543349859196	0.999991172854031	0.999948318010881	0.999294771778394	0.999281992128727	0.999648972123615	0.999830535553711	0.999878374906562	0.999018394363832	0.999538170656178	0.999275324831939	0.999687090291690	0.999849996675814	0.999612170351687	0.999802613053861	0.999966051978048	0.999898425840718	0.999638547909385	0.999527080055406	0.999488703846325	0.999981159264836	0.999948502969719	0.999909681791673	0.999941316690480	0.999653798950690	0.999852212251292	0.999811252295677	0.999704132767466	0.999906824641165	0.999997840124585	0.999963616703160	0.999986336886099	0.999977830034922	0.999907829565272	0.999812681422242	0.999758497629026	0.999892985406720	0.999881564103445	0.999804603275232	0.999852676020266	0.999980893541681	0.999988431878705	0.999927454841520	0.999923469736904	0.999911924006698	0.999936192898009	0.999986553188117	0.999927521858735	0.999914443388919	0.999993671847077	0.999912617568527	0.999989573684635	0.999945530557073	0.999922818361644];
    fitacc = reshape(fitacc,101,107); % fitting accuracy as a function of lifetime and SNR
    % figure; contourf(snrsim,ltsim,fitacc); set(gca,'YScale','log');

    % Match SNR and lifetime
    [~,snrind] = min(abs(SNR-snrsim));
    [~,ltind] = min(abs(LIFETIME-ltsim));

    % Get accuracy and convert to error
    currentaccuracy = fitacc(ltind,snrind);
    currenterror = 1-currentaccuracy; % relative error
    ERROR = currenterror*LIFETIME; % absolute error
    return
end

% Normalize a signal
function NORMSIG = normsweep(SIGNAL)
%%% This function background-subtracts and normalizes a vector.
    
    NORMSIG = (SIGNAL-min(SIGNAL))./(max(SIGNAL)-min(SIGNAL));
    return
end

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

% Fitting cell-silent BonFIRE spectra with a Fano fit for NDR-TPA
function [NX,NY,SOLN,FWHM] = pkfitnt(XVAL,YVAL)
%%% This function fits y(x) to a Gaussian with a linear background and a
% Fano resonance.
    
    % Ensure column vectors
    if length(XVAL) == width(XVAL)
        XVAL = XVAL.';
    end
    if length(YVAL) == width(YVAL)
        YVAL = YVAL.';
    end
    % Make sure lengths match
    if length(XVAL) ~= length(YVAL)
        if length(XVAL) < length(YVAL) % YVAL is longer
            YVAL = YVAL(1:length(XVAL));
        else % XVAL is longer
            XVAL = XVAL(1:length(YVAL));
        end
    end

    % Set up fitting function
    fn = @(x) x(1)+x(2)*XVAL+... % linear baseline
        x(3).*exp(-x(4).*(XVAL-x(5)).^2)+... % Gaussian (BonFIRE)
        x(6).*((x(7).*x(8)+XVAL-x(5)).^2)./(x(8).^2+(XVAL-x(5)).^2)+... % Fano (NDR-TPA)
        -1.*YVAL;
    
    % Set options
    opt=optimoptions(@lsqnonlin);
    opt.MaxFunctionEvaluations = 1e6; % let the fit run longer
    opt.MaxIterations = 1e6; % let the fit run longer
    opt.FunctionTolerance = abs(max(YVAL)-min(YVAL)).*1e-12; % make the fit more accurate
    opt.OptimalityTolerance = abs(max(YVAL)-min(YVAL)).*1e-12; % make the fit more accurate
    opt.StepTolerance = abs(max(YVAL)-min(YVAL)).*1e-12; % make the fit more precise
    opt.Display = 'off'; % silence console output
    
    [~,maxind] = max(YVAL);

    % Set initial guesses
    x0 = [-1 (YVAL(end)-YVAL(1))./(XVAL(end)-XVAL(1))  ...
        0.6148 0.0096 XVAL(maxind) ...
        -1 .44 10];
    lb = [-Inf 0 ...
        0  0.0025 -Inf ... % 0.0025 --> max FWHM = 33.3 cm-1
        -1 0 0];
    ub = [Inf Inf ...
        Inf 0.04 Inf ... % 0.04 --> min FWHM = 8.3 cm-1
        0 3 10];
    
    % Calculate fit
    SOLN = lsqnonlin(fn,x0,lb,ub,opt);
    FWHM = 2*sqrt(2*log(2))*(sqrt(1./(2*SOLN(4))));
    
    % Calculate function for plotting
    NX = min(XVAL):0.1:max(XVAL);
    NY = SOLN(1)+SOLN(2)*NX+SOLN(3).*exp(-SOLN(4).*(NX-SOLN(5)).^2)+...
        SOLN(6).*((SOLN(7).*SOLN(8)+NX-SOLN(5)).^2)./(SOLN(8).^2+(NX-SOLN(5)).^2);
    return
end