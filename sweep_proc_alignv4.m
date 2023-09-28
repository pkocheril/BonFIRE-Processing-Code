%%% Batch processing and alignment of 2D temporal sweep data
%%% initially written from sweep_processv10
%%% v2 - added IR sweep processing, heaviside function fitting, multi-peak
% fitting, configurable number of beating frequencies, writing out powers,
% cleaned up lifetime logic, improved .txt handling, automated DFG/idler
% assignment
%%% v3 - padding signal for fitting, fitting negative peaks, SPCM + PMT
%%%%%%%%%%%%
%%% v4 - updated fitting based on fit_temp_sweepv1 (big overhaul!)

% Initialize
cd '/Users/pkocheril/Documents/Caltech/Wei Lab/Data/2023_09_28'
clear; clc; close all;

% Configuration options - check before running!
loadprevious = 0; % 0 = new analysis, 1 = load previous, [] = auto-detect
testrun = 1; % 0 = process all, 1 = test run with a few files,
% 2 = examine a single file
testfilesperfolder = 1; % number of files to test per folder (default 1)
targetfolders = 1:13; % indices of folders to process, [] = dialog
filetype = 2; % 1 = .txt, 2 = .tif (solution), 3 = .tif (image)
filenameconv = 2; % 0 = other,
% 1 = [pumpWL]-[pumppower]-[signalWL]-[IRpower]-[ND]-[PMTgain]-[modfreq]-[PMTBW]-[etc],
% 2 = [etc]_[size]_[idler]_[DFG]_[power]_[channel],
% 3 = [IRWN]_[probeWL]_[etc]_[size]_[idler]_[DFG]_[power]_[channel],
% 4 = [probeWL]_[IRWN]_[etc]_[size]_[idler]_[DFG]_[power]_[channel]
t0pos = []; % specify t0 position (mm), [] = autofind
pairedDCAC = 0; % 2 = SPCM + PMT CH2, 1 = PMT CH1 + CH2, 0 = only PMT CH2
writeprocyn = 1; % 1 = write batch processed files, 0 = not
writefigsyn = 1; % 1 = write figure files, 0 = not
powernormyn = 0; % 0 = no normalization, 1 = normalize by IR power,
% 2 = normalize by probe and IR powers, 3 = 2 + PMT gain correction
tempmod = 0; % 1 = temporal modulation in place, 0 = no beamsplitters
ND = 2; % ND filter strength in probe path (script will update if possible)
PMT = 1; % PMT gain (norm to 1, script will update if possible)
prbpowerset = 250; % probe power in mW (script will update if possible)
IRpowerset = 70; % IR power in mW (script will update if possible)
normIRpower = 50; % IR power on-sample to normalize to (default is 50)
normprobepower = 1.5; % probe power on-sample to normalize to (default 1.5)
trimlastT = 0; % how many points to remove from end of Tlist (default 0)
trimfirstT = 0; % points to remove from start of Tlist (default 0)
basefittype = 3; % 0 = no baseline fit, 1 = linear, 2 = exp, 3 = exp+line
floatbase = 0.1; % fraction to float baseline coeffs (0.1 -> +/- 10%)
cutlow = 0.05; % lower baseline fit cutoff, (default 5th %ile)
cuthigh = 0.95; % upper baseline fit cutoff, (default 95th %ile)
ltfittype = []; % [] = auto-choose, 0 = no fitting, 1 = Gaussian*monoexp,
% 2 = Gaussian*biexp, 3 = Gauss*stretchexp
setpulsewidth = []; % define pulse width (ps) in fit, [] = float
ltmin = [];  % define minimum lifetime (ps) in fit, [] = default (0.1)
ltmax = []; % define maximum lifetime (ps) in fit, [] = default (100)
tuningrate = 0.033; % PicoEmerald signal t0 tuning rate (ps/nm)
minprobe = []; maxprobe = []; % probe range for aligning tD, [] = auto-find
annotatefigure = 1; % 1 = annotate figure (unfinished, leave as 1)
troubleshoot = 0; % 0 = default, 1 = show everything along the way

D = pwd; % get current directory
S = dir(fullfile(D,'*')); % search current directory
N = setdiff({S([S.isdir]).name},{'.','..'}); % subfolders of D

if isempty(loadprevious) == 1 % auto-detect previous analysis
    ymllist = dir(fullfile(D,'*.yml')); % look for .yml files
    if isempty(ymllist) == 1 % no .yml files
        loadprevious = 0;
    else % found a .yml file
        if length(ymllist) == 1 % found 1 .yml file
            loadprevious = 1;
            testrun = 0; % not a test run if loading previous analysis
        else
            error("Multiple YMLs found - can't load previous analysis")
        end
    end
end

if troubleshoot == 1
    testrun = 2;
end

if testrun == 2 % to examine a single file - choose folder
    [indx,~] = listdlg('PromptString',{'Select one folder.'},...
    'SelectionMode','single','ListString',N);
    targetfolders = indx;
end

if isempty(targetfolders) == 1 % dialog to select folders
    [indx,~] = listdlg('PromptString',{'Select folders to process.'},...
    'SelectionMode','multiple','ListString',N);
    targetfolders = indx;
end

if testrun > 0  % setup test run conditions
    visibility = 'on'; % show figures
    writeprocyn = 0; % don't write batch processed files
    writefigsyn = 0; % don't write figure files
    if length(targetfolders) >= 3 % folder logic for test run
        folders = targetfolders(1:3);
    else
        folders = targetfolders;
    end
    if powernormyn > 2
        powernormyn = 2;
    end
else % not a test run
    visibility = 'off'; % hide figures
    if isempty(targetfolders) == 0 % if target folders are specified
        folders = targetfolders;
    else
        folders = 1:numel(N); % all folders
    end
end

% Column indices for master array
indct = 1; % "indc" = column index
indccorrsig = indct+1; % 2
indcvalue = indccorrsig+1; % 3 - individual values
indcrawsig = indcvalue+1; % 4
indcsigbase = indcrawsig+1; % 5
indcsnrsweep = indcsigbase+1; % 6
indctint = indcsnrsweep+1; % 7
indcsigint = indctint+1; % 8
indcfitcurve = indcsigint+1; % 9
indcfitval = indcfitcurve+1; % 10
indcncsig = indcfitval+1; % 11
indctalign = indcncsig+1; % 12
indcsigalign = indctalign+1; % 13
indcrawDC = indcsigalign+1; % 14
indccorrDC = indcrawDC+1; % 15
indcDCbase = indccorrDC+1; % 16
% Row indices for individual values in column 3
indrIRWN = 1; % "indr" = row index
indrprobe = indrIRWN+1; % 2
indrsigpeak = indrprobe+1; % 3
indrpeaksd = indrsigpeak+1; % 4
indrintegsig = indrpeaksd+1; % 5
indrintegsd = indrintegsig+1; % 6
indrsnr = indrintegsd+1; % 7
indrsbr = indrsnr+1; % 8
indrnoise = indrsbr+1; % 9
indrbasesd = indrnoise+1; % 10
indrintsnr = indrbasesd+1; % 11
indrbleach = indrintsnr+1; % 12
indrtlist = indrbleach+1; % 13
indrlifetime1 = indrtlist+1; % 14
indrlifetime2 = indrlifetime1+1; % 15
indrlt1lt2 = indrlifetime2+1; % 16
indrnfiles = indrlt1lt2+1; % 17
indralignlength = indrnfiles+1; % 18
indrIRpower = indralignlength+1; % 19
indrprbpower = indrIRpower+1; % 20
indrDCpeak = indrprbpower+1; % 21
indrDCpeaksd = indrDCpeak+1; % 22
indrintegDC = indrDCpeaksd+1; % 23
indrDCintegsd = indrintegDC+1; % 24
indrDCsnr = indrDCintegsd+1; % 25
indrDCsbr = indrDCsnr+1; % 26
indrDCnoise = indrDCsbr+1; % 27
indrDCbleach = indrDCnoise+1; % 28 - minimum value of maxTlength

if loadprevious == 0 % run new analysis
    if filetype == 3 % for image data
        writefigsyn = 0; % don't write individual figures
        if powernormyn > 2 % don't do photobleach correction
            powernormyn = 2;
        end
    end
    
    % Pre-loop to figure out dimensions of master array
    maxTlength = 0; maxfilesinfolder = 0; % counter variables
    for ii = folders % pre-loop to figure out dimensions of master array
        if filetype == 1  % look for .txt files in subfolders
            T = dir(fullfile(D,N{ii},'*.txt'));
        else % look for CH2 tifs
            T = dir(fullfile(D,N{ii},'*CH2*.tif'));
        end
        C = {T(~[T.isdir]).name}; % all data files
        if numel(C) > maxfilesinfolder
            maxfilesinfolder = numel(C);
        end
        if isempty(T) == 0 % skips subfolders without data files
            for jj = 1:numel(C) % jj = file number in subfolder ii
                F = fullfile(D,N{ii},C{jj}); % current file
                if filetype == 1 % load data from .txt
                    data = importdata(F); % load data
                    x = data(:,1); % save delay position as x
                else % load data from Tlist.txt and .tif
                    Tlistname = fullfile(D,N{ii},'Tlist.txt'); % current Tlist
                    x = importdata(Tlistname); % load Tlist as x
                end
                if x(1) == x(end) && trimlastT == 0 % if x isn't single-valued
                    trimlastT = 1; % trimming is needed
                end
                x = x(trimfirstT+1:end-trimlastT); % cleaned up point trimming
                currentTlength = length(x);
                % Check if raw data work for convolution fitting
                deltax = zeros(length(x)-1,1);
                for i=1:length(x)-1
                    deltax(i) = round(1e2*(x(i+1)-x(i)))/1e2;
                end
                xspacing = unique(deltax);
                if length(xspacing) == 1 && mod(length(x),2) == 1
                    % data ok - do nothing
                else % add extra length for interpolation
                    currentTlength = 2*currentTlength-1;
                end
                if currentTlength > maxTlength
                    maxTlength = currentTlength;
                end
            end
        end
    end
    
    if 10*maxTlength-1 < indrDCbleach % make sure master array is more than large enough
        maxTlength = indrDCbleach;
    else
        maxTlength = 10*maxTlength-1; % interpolate for tD compensation
    end
    
    if testrun == 2 % to examine a single file - choose file
        [indx2,~] = listdlg('PromptString',{'Select a file to process.'},...
        'SelectionMode','single','ListString',C);
        targetfile = indx2;
    end
    
    % Make master array to store all data (row,column,folder,file)
    if pairedDCAC > 0
        master = NaN(maxTlength,indcDCbase,length(folders),maxfilesinfolder);
    else  % remove 3 columns for DC raw, corrected, & baseline fit
        master = NaN(maxTlength,indcDCbase-3,length(folders),maxfilesinfolder);
    end
    master_size = size(master);

    if powernormyn > 0 % prep for power normalization
        IRpower = []; prbpower = [];
        if tempmod == 1 % pre-load beamsplitter response curve if needed
            BSdata = importdata('/Users/pkocheril/Documents/Caltech/Wei Lab/Spreadsheets/BP145B2/BP145B2.csv');
        end
    end
    
    % Loop through all files and analyze
    for ii = folders % ii = subfolder number
        if filetype == 1  % look for .txt files in subfolders
            T = dir(fullfile(D,N{ii},'*.txt'));
            infofilesraw = [];
            infofiles = [];
        else % look for CH2 tifs
            T = dir(fullfile(D,N{ii},'*CH2*.tif'));
            infofilesraw = dir(fullfile(D,N{ii},'*IRsweep*.txt'));
            infofiles = {infofilesraw(~[infofilesraw.isdir]).name}; % info files
        end
        C = {T(~[T.isdir]).name}; % all data files
        if numel(C) == numel(infofiles)
            infofound = 1;
        else
            infofound = 0;
        end
        if testrun > 0 % limit test run files per folder
            totalfilesC = testfilesperfolder;
        else
            totalfilesC = numel(C);
        end
        if testrun == 2 % load only the chosen file
            totalfilesC = targetfile;
            startjj = targetfile;
        else
            startjj = 1;
        end
        if isempty(T) == 0 % skip subfolders without data files
            for jj = startjj:totalfilesC % jj = file number in subfolder ii
                F = fullfile(D,N{ii},C{jj}); % current file
                IRpower = []; prbpower = [];
                
                % Filename parsing
                folderinfo = split(F,"/"); % split path into folders
                if filenameconv == 1 % parse .txts
                    fileinfo = split(folderinfo(end),"-"); % split filename at "-"s
                    subfolderinfo = split(N{ii}," "); % split subfolder by spaces
                    prbWLstrfull = char(fileinfo(1)); % extract pump WL - "760m"
                    prbWLstr = prbWLstrfull(1:end-1); % remove "m" suffix - "760"
                    prbWL = sscanf(prbWLstr,'%f'); % convert to double
                    prbpowerstrfull = char(fileinfo(2)); % extract pump power
                    prbpowerstr = prbpowerstrfull(1:end-2); % remove "mW"
                    prbpower = sscanf(prbpowerstr,'%f');
                    signalWLstrfull = char(fileinfo(3)); % extract signal WL
                    signalWLstr = signalWLstrfull(1:end-2); % remove "nm"
                    signalWL = sscanf(signalWLstr,'%f');
                    IRpowerstrfull = char(fileinfo(4)); % extract IR power
                    IRpowerstr = IRpowerstrfull(1:end-2); % remove "mW"
                    IRpower = sscanf(IRpowerstr,'%f');
                    DFGWL = 1/((2/signalWL) - (1/1031.2)); % get DFG WL (nm)
                    idlerWL = 1/((1/1031.2)-(1/signalWL));
                    DFGWN = 1e7/DFGWL;
                    idlerWN = 1e7/idlerWL;
                    if idlerWL < 2600 % nm
                        IRWN = DFGWN;
                    else
                        IRWN = idlerWN;
                    end
                    NDstrfull = char(fileinfo(5)); % extract ND
                    NDstr = NDstrfull(3:end); % remove "ND"
                    ND = sscanf(NDstr,'%f');
                    PMTstrfull = char(fileinfo(6)); % extract PMT gain
                    PMTstr = PMTstrfull(4:end); % remove "PMT"
                    PMT = sscanf(PMTstr,'%f');
                    modfreqstrfull = char(fileinfo(7)); % extract modulation freq
                    modfreqstr = modfreqstrfull(1:end-3); % remove "kHz" or "MHz"
                    if modfreqstrfull(end-2) == "k" % kHz
                        modfreq = 1e3*sscanf(modfreqstr,'%f');
                    else % MHz
                        modfreq = 1e6*sscanf(modfreqstr,'%f');
                    end
                    PMTBWstrfull = char(fileinfo(8)); % extract PMT bandwidth
                    PMTBWstr = PMTBWstrfull(1:end-3); % remove "kHz"
                    if PMTBWstrfull(end-2) == "k" % kHz
                        PMTBW = 1e3*sscanf(PMTBWstr,'%f');
                    else % MHz
                        PMTBW = 1e6*sscanf(PMTBWstr,'%f');
                    end
                end
                if filenameconv == 2 % parse default .tifs
                    fileinfo = split(folderinfo(end),"_"); % split at "_"s
                    channel = char(fileinfo(end)); % 'CH2'
                    IRpowermeterstrfull = char(fileinfo(end-1)); % 'P0.4409'
                    IRpowermeterstr = IRpowermeterstrfull(2:end); % '0.4409'
                    IRpowermeter = sscanf(IRpowermeterstr,'%f'); % 0.4409
                    DFGWNstrfull = char(fileinfo(end-2)); % 'DFG3797.9'
                    DFGWNstr = DFGWNstrfull(4:end); % '3797.9'
                    DFGWN = sscanf(DFGWNstr,'%f'); % 3797.9
                    idlerWNstrfull = char(fileinfo(end-3)); % 'Idler2949.8'
                    idlerWNstr = idlerWNstrfull(6:end); % '2949.8'
                    idlerWN = sscanf(idlerWNstr,'%f'); % 2949.8
                    %imgdim = char(fileinfo(end-4)); % '5X5'
                    prbWL = 1; % <-- unknown from filename
                    idlerWL = 1E7/idlerWN;
                    signalWL = 1/((1/1031.2)-(1/idlerWL));
                    if IRpowermeter ~= 0
                        IRpower = IRpowermeter*300; % assuming 300 mW scale
                    end
                    if idlerWL < 2600 % nm
                        IRWN = DFGWN;
                    else
                        IRWN = idlerWN;
                    end
                end
                if infofound == 1 % parse info from IRsweep.txt files
                    infofilename = fullfile(D,N{ii},infofiles{jj});
                    infooptions = detectImportOptions(infofilename);
                    infooptions = setvaropts(infooptions,'Var1','InputFormat','MM/dd/uuuu');
                    infofiletable = readtable(infofilename,infooptions);
                    infofile = table2array(infofiletable(:,3:end)); % cut date/time
                    xinitial = infofile(1);
                    xfinal = infofile(2);
                    xsteps = infofile(3);
                    yinitial = infofile(4);
                    yfinal = infofile(5);
                    ysteps = infofile(6);
                    zinitial = infofile(7);
                    zfinal = infofile(8);
                    zsteps = infofile(9);
                    tinitial = infofile(10);
                    tfinal = infofile(11);
                    tsteps = infofile(12);
                    if round(infofile(13)) ~= 0 % for frequency doubling
                        prbWL = round(infofile(13));
                    else
                        prbWL = 515.6;
                    end
                    if round(infofile(14)) ~= 0
                        prbpower = infofile(14); % set to prbpowerset later if needed
                    end
                    ND = infofile(15);
                    idlerWN = round(infofile(16));
                    DFGWN = round(infofile(17));
                    if infofile(18) > 0
                        IRpower = infofile(18);
                    end
                    idlerWL = 1e7/idlerWN;
                    if idlerWL < 2600 % nm
                        IRWN = DFGWN;
                    else
                        IRWN = idlerWN;
                    end
                    dwelltime = infofile(19);
                end
                if filenameconv == 3 % parse probe sweep .tifs
                    subfolderinfo = split(N{ii}," "); % split subfolder at " "s
                    fileinfo = split(folderinfo(end),"_"); % split filename at "_"s
                    IRWNstrfull = char(fileinfo(1)); % extract IRWN - "1598cm-1"
                    IRWNstr = IRWNstrfull(1:end-4); % remove "cm-1" suffix - "1598"
                    IRWN = sscanf(IRWNstr,'%f'); % convert to double - 1598
                    prbstrfull = char(fileinfo(2)); % probe WL - "780 (10)"
                    if length(prbstrfull) >= 6 % "780 (10)" --> 790 nm
                        prbsplit = split(prbstrfull," "); % "780" "(10)"
                        if length(prbsplit) > 1
                            prbWL1char = char(prbsplit(1)); % '780'
                            prbWL1 = sscanf(prbWL1char,'%f'); % 780
                            prbWL2char = char(prbsplit(2)); % '(10)' or 'finalcheck'
                            if length(prbWL2char) >= 7 % 'finalcheck'
                                prbWL2 = 0;
                            else % '(10)'
                                prbWL2str = prbWL2char(2:end-1); % '10'
                                prbWL2 = sscanf(prbWL2str,'%f'); % 10
                            end
                            if totalfilesC <= 130
                                prbWL2 = 2*prbWL2; % step size 2 nm
                            end
                            prbWL = prbWL1 + prbWL2; % 790
                        else % "780finalcheck"
                            prbstr = prbstrfull(1:3); % '780'
                            prbWL = sscanf(prbstr,'%f'); % 780
                        end
                    else % "780" or "780nm"
                        prbstr = prbstrfull(1:3); % '780'
                        prbWL = sscanf(prbstr,'%f'); % 780
                    end
                    IRpowermeterstrfull = char(fileinfo(end-1)); % 'P0.4409'
                    IRpowermeterstr = IRpowermeterstrfull(2:end); % '0.4409'
                    IRpowermeter = sscanf(IRpowermeterstr,'%f'); % 0.4409
                    DFGWNstrfull = char(fileinfo(end-2)); % 'DFG3797.9'
                    DFGWNstr = DFGWNstrfull(4:end); % '3797.9'
                    DFGWN = sscanf(DFGWNstr,'%f'); % 3797.9
                    idlerWNstrfull = char(fileinfo(end-3)); % 'Idler2949.8'
                    idlerWNstr = idlerWNstrfull(6:end); % '2949.8'
                    idlerWN = sscanf(idlerWNstr,'%f'); % 2949.8
                    %imgdim = char(fileinfo(end-4)); % '5X5'
                    idlerWL = 1E7/idlerWN;
                    signalWL = 1/((1/1031.2)-(1/idlerWL));
                    if IRpowermeter ~= 0
                        IRpower = IRpowermeter*300; % assuming 300 mW scale
                    end
                end
                if filenameconv == 4 % parse IR sweep .tifs (unnecessary if IR checked)
                    subfolderinfo = split(N{ii}," "); % split subfolder at " "s
                    fileinfo = split(folderinfo(end),"_"); % split filename at "_"s
                    prbWLstrfull = char(fileinfo(1)); % extract probe WL - "780nm"
                    prbWLstr = prbWLstrfull(1:end-2); % remove "nm" suffix - "780"
                    prbWL = sscanf(prbWLstr,'%f'); % convert to double - 780
                    IRWNstrfull = char(fileinfo(2)); % IRWN - "1000cm-1 (10)"
                    IRWNsplit = split(IRWNstrfull," "); % "1000cm-1" "(10)"
                    IRWN1char = char(IRWNsplit(1)); % '1000cm-1'
                    IRWN1chartrim = IRWN1char(1:end-4); % "1000"
                    IRWN1 = sscanf(IRWN1chartrim,'%f'); % 1000
                    IRWN2char = char(IRWNsplit(2)); % "(10)"
                    IRWN2chartrim = IRWN2char(2:end-1); % "10"
                    IRWN2 = sscanf(IRWN2chartrim,'%f'); % 10
                    IRWN = IRWN1+10*IRWN2; % 1100
                    IRpowermeterstrfull = char(fileinfo(end-1)); % 'P0.4409'
                    IRpowermeterstr = IRpowermeterstrfull(2:end); % '0.4409'
                    IRpowermeter = sscanf(IRpowermeterstr,'%f'); % 0.4409
                    DFGWNstrfull = char(fileinfo(end-2)); % 'DFG3797.9'
                    DFGWNstr = DFGWNstrfull(4:end); % '3797.9'
                    DFGWN = sscanf(DFGWNstr,'%f'); % 3797.9
                    idlerWNstrfull = char(fileinfo(end-3)); % 'Idler2949.8'
                    idlerWNstr = idlerWNstrfull(6:end); % '2949.8'
                    idlerWN = sscanf(idlerWNstr,'%f'); % 2949.8
                    %imgdim = char(fileinfo(end-4)); % '5X5'
                    idlerWL = 1E7/idlerWN;
                    signalWL = 1/((1/1031.2)-(1/idlerWL));
                    if IRpowermeter ~= 0
                        IRpower = IRpowermeter*300; % assuming 300 mW scale
                    end
                end

                if isempty(IRpower) == 1 % if no IR power found
                    IRpower = IRpowerset; % set to power from config
                end
                if isempty(prbpower) == 1 % if no probe power found
                    prbpower = prbpowerset; % set to power from config
                end
    
                % Load data
                if filetype == 1 % load data from .txt
                    data = importdata(F); % load data
                    x = data(:,1); % save delay position as x
                    if pairedDCAC > 0
                        DC = data(:,2); % save channel 1 as DC
                        signal = data(:,3); % save channel 2 as signal
                    else
                        if width(data) >= 3
                            signal = data(:,3); % save channel 2 as signal
                        else
                            signal = data(:,2);
                        end
                        DC = zeros(height(signal),width(signal));
                    end
                    if x(1) == x(end) && trimlastT == 0 % if x isn't single-valued
                        trimlastT = 1; % trimming is needed
                    end
                    x = x(trimfirstT+1:end-trimlastT);
                    DC = DC(trimfirstT+1:end-trimlastT);
                    signal = signal(trimfirstT+1:end-trimlastT);
                    % No standard deviation data present - write as zeros
                    sigsds = zeros(height(signal),width(signal));
                    if pairedDCAC > 0
                        DCsds = zeros(height(signal),width(signal));
                    end
                    imagesize = [1 1];
                else % load Tlist and .tif
                    Tlistname = fullfile(D,N{ii},'Tlist.txt'); % current Tlist
                    x = importdata(Tlistname); % load Tlist
                    data = double(tiffreadVolume(F));
                    imagesize = size(data);
                    if pairedDCAC > 0
                        if pairedDCAC == 2 % SPCM + CH2
                            F_DC = F;
                            F_DC(end-6:end+1) = "SPCM.tif";
                            DCdata = double(tiffreadVolume(F_DC));
                        else % pairedDCAC = 1, PMT CH1 + CH2
                            F_DC = F;
                            F_DC(end-4) = "1";
                            DCdata = double(tiffreadVolume(F_DC));
                        end
                    else % no paired DCAC
                        DCdata = zeros(height(data),width(data),length(x));
                    end
                    if x(1) == x(end) && trimlastT == 0 % if x isn't single-valued
                        trimlastT = 1; % trimming is needed
                    end
                    x = x(trimfirstT+1:end-trimlastT);
                    data = data(:,:,trimfirstT+1:end-trimlastT);
                    DCdata = DCdata(:,:,trimfirstT+1:end-trimlastT);
                    if length(x) ~= length(data)
                        warning('Tlist mismatch - trimming to correct')
                        if length(x) > length(data)
                            x = x(1:length(data));
                        else
                            data = data(:,:,1:length(x));
                        end             
                    end
                    if filetype == 2 % average & stdev of XY data in .tif
                        signal = squeeze(mean(data, [1 2]));
                        sigsds = squeeze(std(data, 0, [1 2]));
                        DC = squeeze(mean(DCdata, [1 2]));
                        DCsds = squeeze(std(DCdata, 0, [1 2]));
                        imagesize = size(signal);
                    else % image data - make arrays to write outputs into
                        r2array = zeros([imagesize(1) imagesize(2)]);
                        lt1array = r2array; lt2array = r2array;
                        ssresidarray = r2array; FWHMarray = r2array;
                    end
                end
                
                if testrun > 0 % do test run in the middle of an image
                    startiii = round(imagesize(1)/2);
                    startjjj = round(imagesize(2)/2);
                    imagesize(1) = (2-1)+startiii; 
                    imagesize(2) = (1-1)+startjjj;
                else
                    startiii = 1; startjjj = 1;
                end

                if filetype ~= 3 % solution data
                    imagesize(1) = 1; imagesize(2) = 1;
                    startiii = 1; startjjj = 1;
                end
    
                for iii=startiii:imagesize(1)
                    for jjj=startjjj:imagesize(2)
                        if filetype == 3 % load pixel data
                            signal = squeeze(data(iii,jjj,:));
                            sigsds = zeros(height(signal),width(signal));
                            if pairedDCAC > 0
                                DC = squeeze(DCdata(iii,jjj,:));
                                DCsds = zeros(sigsds);
                            else
                                DC = sigsds; DCsds = sigsds;
                            end
                        end
                        if troubleshoot == 1
                            figure; errorbar(x,signal,sigsds,'o'); title('Signal vs position');
                        end
                
                        % Converting to time
                        cmmps = 299792458*1E3*1E-12; % c in mm/ps
                        t = zeros([length(x),1]);
                        if isempty(t0pos) == 1 % guess t0 by peak
                            indexoffset = floor(length(x)/6); % using an offset to ignore sharp decay at start
                            [sigmax, t0index] = max(signal(indexoffset:end)); % estimating t0 by peak
                            t(:) = (x(:)-x(t0index+indexoffset-1))*4/cmmps; % time vector in ps
                        else % defined t0 position
                            t(:) = (x(:)-t0pos)*4/cmmps; % time vector in ps
                        end

                        if powernormyn > 0 % power normalization
                            if isempty(IRpower) == 1 % if no IR power found
                                IRpower = IRpowerset; % set to power from config
                            end
                            if isempty(prbpower) == 1 % if no probe power found
                                prbpower = prbpowerset; % set to power from config
                            end
                            signal = signal*normIRpower/IRpower;
                            sigsds = sigsds*normIRpower/IRpower;
                            DC = DC*normIRpower/IRpower;
                            DCsds = DCsds*normIRpower/IRpower;

                            if powernormyn > 1 % 2 = IR and probe
                                prbpower = 0.25*prbpower; % loss from pinhole, 250 mW -> 62.5 mW
                                if tempmod == 1 % correct for temporal modulation
                                    transmittance = BSdata.data(prbWL,1)/100; % picoEmerald is p-polarized (horizontal)
                                    reflectance = BSdata.data(prbWL,4)/100; % CSV in %
                                    prbpower = transmittance*reflectance*prbpower; % 62.5 mW -> 15.6 mW
                                end
                                prbpower = prbpower/(10^(ND)); % ND1: 15.6 mW -> 1.56 mW
                                signal = signal*(normprobepower/prbpower);
                                sigsds = sigsds*(normprobepower/prbpower);
                                if pairedDCAC > 0
                                    DC = DC*(normprobepower/prbpower);
                                    DCsds = DCsds*(normprobepower/prbpower);
                                end
                                if powernormyn > 2 % 3 = PMT gain correction
                                    gaincurve = importdata('/Users/pkocheril/Documents/Caltech/Wei Lab/Data/2023_04_28/PMT Gain Data/gain_calibration_curve.txt');
                                    gainfactor = gaincurve(PMT,2);
                                    signal = signal/gainfactor;
                                    sigsds = sigsds/gainfactor;
                                    if pairedDCAC > 0
                                        DC = DC/gainfactor;
                                        DCsds = DCsds/gainfactor;
                                    end
                                end
                            end
                        end
                        
                        % Setup baseline fitting
                        ilow = ceil(length(t)*cutlow); % using ceiling to round up
                        ihigh = ceil(length(t)*cuthigh);
                        tbase = [t(2:ilow); t(ihigh:end)]; % trimmed t
                        sigbase = [signal(2:ilow); signal(ihigh:end)]; % trimmed signal
                        sbr = max(signal(ilow:ihigh))/mean(signal(ihigh:end));
                        if pairedDCAC > 0
                            DCbase = [DC(2:ilow); DC(ihigh:end)]; % trimmed DC
                            DCsbr = max(DC(ilow:ihigh))/mean(DC(ihigh:end));
                        end
                        
                        % Baseline fitting
                        if basefittype == 0 % no baseline fitting, set basecurve to 0
                            basecurve = zeros(height(t),width(t));
                            bleachrate = 0;
                            if pairedDCAC > 0
                                DCbasecurve = zeros(height(t),width(t));
                                DCbleachrate = 0;
                                DCsbr = 0;
                            end
                            basefit = zeros(1,5); lbbase = basefit; ubbase = basefit;
                        else
                            basefn = @(b) b(1)+b(2)*exp(-b(3)*(tbase-b(4)))+b(5)*tbase-sigbase; % error function
                            lbbase = [-Inf -Inf -Inf -Inf -Inf];
                            ubbase = [Inf Inf Inf Inf Inf];
                            if basefittype == 1
                                lbbase(2:4) = 0; ubbase(2:4) = 0;
                            end
                            if basefittype == 2
                                lbbase(5) = 0; ubbase(5) = 0;
                            end
                            baseopt=optimoptions(@lsqnonlin);
                            baseopt.MaxFunctionEvaluations = 1e6; % let the fit run longer
                            baseopt.MaxIterations = 1e6; % let the fit run longer
                            baseopt.FunctionTolerance = 1e-12; % make the fit more accurate
                            baseopt.OptimalityTolerance = 1e-12; % make the fit more accurate
                            baseopt.StepTolerance = 1e-12; % make the fit more precise
                            baseopt.Display = 'off'; % silence console output
                            basegs = [min(sigbase) 0.1*(max(sigbase)-min(sigbase)) 0 0 0];
                            basefit = lsqnonlin(basefn,basegs,lbbase,ubbase,baseopt);
                            basecurve = basefit(1)+basefit(2)*exp(-basefit(3)*(t-basefit(4)))+...
                                basefit(5)*t;
                            if basefit(2) > 0 && basefit(3) > 0 
                                bleachrate = basefit(3); % ps-1
                            else % bleach rate doesn't have physical meaning
                                bleachrate = 0;
                            end
                            if pairedDCAC > 0 % fit DC baseline
                                DCbasefn = @(b) b(1)+b(2)*exp(-b(3)*(tbase-b(4)))+b(5)*tbase-DCbase;
                                DCbgs = [min(DCbase) 0.1*(max(DCbase)-min(DCbase)) 0 0 0];
                                DCbasefit = lsqnonlin(DCbasefn,DCbgs,[],[],baseopt);
                                DCbasecurve = DCbasefit(1)+DCbasefit(2)*exp(-DCbasefit(3)*(t-DCbasefit(4)))+...
                                    DCbasefit(5)*t;
                                if DCbasefit(2) > 0 && DCbasefit(3) > 0
                                    DCbleachrate = DCbasefit(3); % ps-1
                                else
                                    DCbleachrate = 0;
                                end
                            end
                            if troubleshoot == 1
                                % Check baseline fitting
                                figure; plot(t,signal,'o',t,basecurve,'-'); title('Baseline fitting');
                            end
                        end

                        % Calculate peak height and integrated signal
                        corrsig = signal - basecurve;
                        % Determine if signal peak is negative or positive
                        abscorrsig = abs(corrsig);
                        flipcorrsig = -1*corrsig;
                        corrabsdist = sum((corrsig-abscorrsig).^2);
                        flipabsdist = sum((flipcorrsig-abscorrsig).^2);
                        if corrabsdist < flipabsdist
                            [sigpeak,maxindex] = max(corrsig); % corrsig is closer to abscorrsig
                            peaksign = 1; % peak is positive
                        else
                            [sigpeak,maxindex] = min(corrsig); % flipcorrsig is closer to abscorrsig
                            peaksign = -1; % peak is negative
                        end
                        peaksd = sigsds(maxindex);
                        integsig = sum(corrsig);
                        integsd = mean(sigsds,'all');
                        corrsigbase = [corrsig(2:ilow); corrsig(ihigh:end)];
                        noise = range(corrsigbase);
                        basesd = std(corrsigbase);
                        snr = sigpeak/noise;
                        snrsweep = corrsig/noise;
                        intsnr = sum(snrsweep);
                
                        % Calculate t spacing
                        deltat = zeros(length(t)-1,1);
                        for i=1:length(t)-1
                            deltat(i) = round(1e2*(t(i+1)-t(i)))/1e2; % round to nearest 0.01 ps (had issues without rounding)
                        end
                        tspacing = unique(deltat); % length 1 <-> evenly spaced t
                        
                        if isempty(ltfittype) == 1 % auto-choose lifetime fitting
                            emptytype = 1;
                            if IRWN > 2100 && IRWN < 2400
                                ltfittype = 1;
                            else
                                ltfittype = 2;
                            end
                        else
                            emptytype = 0;
                        end
                
                        if ltfittype > 0 % lifetime fitting
                            % Check t spacing (even) and length (odd)
                            if length(tspacing) == 1  && mod(length(t),2) == 1
                                tint = t; sigint = signal; % data are good for fitting
                            else % interpolate for convolution
                                tint = linspace(t(1),t(end),length(t)*2-1).';
                                sigint = spline(t, signal, tint);
                                if troubleshoot == 1
                                    figure; plot(tint,sigint,'o',t,signal,'o'); title('Interpolation')
                                end
                            end
                
                            if max(tint) < 30 % auto-pad signal if tint doesn't go past 30 ps
                                padsignal = 1;
                            else
                                padsignal = 0;
                            end
                
                            % Prepare for signal padding
                            tintspace = tint(2)-tint(1);
                            originallength = length(tint);
                            tlong = linspace(tint(1),tint(end)+tintspace*100,length(tint)+100).';
                            baselong = basefit(1)+basefit(2)*exp(-basefit(3)*(tlong-basefit(4)))+...
                                basefit(5)*tlong;
                            if basefittype > 0
                                siglong = baselong;
                                siglong(1:length(sigint)) = sigint;
                            else
                                if peaksign == 1
                                    siglong(length(sigint)+1:length(tlong)) = min(sigint(end-5:end));
                                else
                                    siglong(length(sigint)+1:length(tlong)) = max(sigint(end-5:end));
                                end
                            end

                            % Signal padding for short-tailed data
                            if padsignal == 1
                                tint = tlong;
                                sigint = siglong;
                                if troubleshoot == 1
                                    figure; plot(tint,sigint,'-o',t,signal,'o'); title('Padding');
                                end
                            end
                
                            tc = tint(1:2:end); % odd values of tint --> convolution input
                            if powernormyn > 0
                                IRFamp = peaksign*IRpower/1e3; % IRF amplitude (~ IR power in W)
                            else
                                IRFamp = peaksign*50/1e3; % default to 50 mW
                            end
                
                            fun = @(r) r(1)+r(2)*exp(-r(3)*(tint-r(4)))+r(5)*tint+ ... % baseline
                                conv(exp(-((tc-r(6))./(r(7)/(2*sqrt(log(2))))).^2), ...
                                IRFamp*heaviside(tc).*(r(8)*exp(-(tc)./r(9))+r(10)*exp(-((tc)./(r(11))).^r(12)))) - sigint;
                            
                            % Initial guesses
                            r0 = [basefit,0,4,... % basecoefs, IRF center (ps), IRF width (ps)
                                0.02,2,0.02,10,1]; % amp 1, τ1 (ps), amp 2, τ2 (ps), β (stretch)
                
                            % Float baseline coeffs
                            newlbbase = basefit; newubbase = basefit;
                            for i=1:length(lbbase)
                                if abs(basefit(i)) < 1e-6
                                    newlbbase(i) = 0;
                                    newubbase(i) = 0;
                                else
                                    newlbbase(i) = min((1-floatbase)*basefit(i),(1+floatbase)*basefit(i));
                                    newubbase(i) = max((1-floatbase)*basefit(i),(1+floatbase)*basefit(i));
                                end
                            end
                
                            % Lower and upper bounds
                            lb = [newlbbase,-10,2,... % basecoefs, IRF center (ps), IRF min (ps)
                                0,0.1,... % amp 1, τ1min (ps)
                                0,0.1,... % amp 2, τ2min (ps)
                                1e-2]; % β min
                            ub = [newubbase,20,9,... % basecoefs, IRF center (ps), IRF max (ps)
                                900,100,... % amp 1, τ1max (ps)
                                900,100,... % amp 2, τ2max (ps)
                                1e2]; % β max

                            % Implement fit type options
                            if ltfittype == 1 % Gauss*monoexp
                                lb(10) = 0; ub(10) = 0; % first exp only
                            end
                            if ltfittype == 2 % Gauss*biexp
                                lb(12) = 1; ub(12) = 1; % no stretching
                            end
                            if ltfittype == 3 % Gauss*stretchexp
                                lb(8) = 0; ub(8) = 0; % second exp only
                            end
                
                            % Parameter constraints from config
                            if isempty(setpulsewidth) == 0 % pulse width was specified
                                lb(7) = setpulsewidth*sqrt(2); ub(7) = setpulsewidth*sqrt(2);
                            end
                            if isempty(ltmin) == 0 % minimum lifetime was specified
                                lb(9) = ltmin; lb(11) = ltmin;
                            end
                            if isempty(ltmax) == 0 % maximum lifetime was specified
                                ub(9) = ltmax; ub(11) = ltmax;
                            end
                
                            % Run lsqnonlin - minimizes error function by adjusting r
                            options=optimoptions(@lsqnonlin);
                            options.MaxFunctionEvaluations = 1e6; % let the fit run longer
                            options.MaxIterations = 1e6; % let the fit run longer
                            options.FunctionTolerance = 1e-12; % make the fit more accurate
                            options.OptimalityTolerance = 1e-12; % make the fit more accurate
                            options.StepTolerance = 1e-12; % make the fit more precise
                            options.Display = 'off'; % silence console output
                            fitval = lsqnonlin(fun,r0,lb,ub,options); % fit coeffs
                
                            fitcurve = fitval(1)+fitval(2)*exp(-fitval(3)*(tint-fitval(4)))+fitval(5)*tint+ ... % baseline
                                conv(exp(-((tc-fitval(6))./(fitval(7)/(2*sqrt(log(2))))).^2), ...
                                IRFamp*heaviside(tc).*(fitval(8)*exp(-(tc)./fitval(9))+fitval(10)*exp(-((tc)./(fitval(11))).^fitval(12))));
                            
                            % Trim padded signal to original length
                            tint = tint(1:originallength);
                            sigint = sigint(1:originallength);
                            fitcurve = fitcurve(1:originallength);
                            
                            % Calculate pulse width from fit
                            IRFwidth = fitval(7); % ps
                            FWHM = IRFwidth/sqrt(2); % ps

                            % Biexponential lifetime assignments
                            lifetime1 = min(fitval(9),fitval(11)); % ps
                            lifetime2 = max(fitval(9),fitval(11)); % ps
                            if abs(fitval(9)) < abs(fitval(11)) % lifetime amplitude ratio
                                lt1lt2 = fitval(8)/fitval(10); % A1/A2
                            else % means fitval(9) is longer (τ2)
                                lt1lt2 = fitval(10)/fitval(8); % still A1/A2
                            end

                            if ltfittype == 1 % Gauss*monoexp
                                lifetime1 = fitval(9); lifetime2 = 0;
                            end

                            if ltfittype == 3 % Gauss*stretchexp
                                lifetime1 = fitval(11); lifetime2 = 0;
                                lt1lt2 = fitval(12); % stretching factor
                            end
                
                            % Calculate residuals and ssresid
                            resid = sigint-fitcurve;
                            ssresid = sum(resid.^2); % sum of squares of residuals
                            tss = sum((sigint-mean(sigint)).^2); % total sum of squares
                            r2 = 1-(ssresid/tss); % coefficient of determination (R^2)
                            if r2 < 0 % remove nonsense R^2 values
                                r2 = 0;
                            end
                
                            % FFT of residuals
                            ts = tint*1e-12; % time, s
                            dt = abs(ts(1)-ts(2)); % time step, s
                            fs = 1/dt; % freq BW, Hz
                
                            resfft = fft(resid); % raw FFT
                            resfftshift = fftshift(resfft); % centered FFT
                            powerfft = abs(resfftshift).^2/length(resid); % power spectrum
                            freqHz = (-length(resid)/2+1/2:length(resid)/2-1/2)*(fs/length(resid)); % FFT frequency vector in Hz
                            wn = freqHz/29979245800; % FFT freq in cm-1
                            
                            if troubleshoot == 1
                                % Check fit
                                figure; tiledlayout(3,1); nexttile([1 1]); plot(tint,resid); nexttile([2 1]); plot(t,signal,'o',tint,fitcurve); title('Fitting')
                                % Check FFT
                                figure; plot(wn,powerfft); title('FFT');
                            end
                        else % no lifetime fitting
                            tint = t; sigint = signal; fitcurve = signal;
                            lifetime1 = 0; lifetime2 = 0; lt1lt2 = 0;
                            fitval = zeros(12,1); r2 = 0; FWHM = 0;
                        end

                        if prbWL < 961 && prbWL > 699 % interpolate for tD compensation
                            alignlength = ceil(((t(end)-t(1))/(tuningrate/2))+1); % want talignspacing >= tuningrate/2
                            talign = linspace(t(1),t(end),alignlength).';
                            sigalign = spline(t, corrsig, talign);
                            talignspacing = talign(2)-talign(1); % evenly spaced

                            % tD compensation
                            if isempty(minprobe) == 1
                                minprobe = 750;
                            end
                            if isempty(maxprobe) == 1
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
                            talign = t; sigalign = signal;
                        end
                        
                        if filetype == 3
                            ssresidarray(iii,jjj) = ssresid;
                            r2array(iii,jjj) = r2;
                            FWHMarray(iii,jjj) = FWHM;
                            lt1array(iii,jjj) = lifetime1;
                            lt2array(iii,jjj) = lifetime2;
                        end
                        
                        % Make rounded strings for figure annotation
                        lt1 = string(round(10*lifetime1)/10); % round to nearest 0.1 ps
                        lt2 = string(round(10*lifetime2)/10);
                        sr2 = string(sprintf('%0.3g',r2)); % round to 3 sig figs
                        slt1lt2 = string(sprintf('%0.3g',lt1lt2));
                        
                        if annotatefigure > 0 % make figure annotation
                            annot = {'IR = '+string(IRWN)+' cm-1','probe = '+...
                                string(prbWL)+' nm'};
                            if ltfittype > 0
                                annotdim = [0.47 0.25 0.2 0.27]; % [posx posy sizex sizey]
                            else
                                annotdim = [0.35 0.6 0.2 0.2];
                            end
                            if pairedDCAC > 0
                                sDCbleachrate = string(sprintf('%0.2g',DCbleachrate)); % round to 2 sig figs
                                annot(length(annot)+1) = {'bleach = '+...
                                    sDCbleachrate+' ps^{-1}'};
                            end
                            if ltfittype == 1
                                annot(length(annot)+1) = {'τ_{mono} = '+lt1+...
                                    ' ps, r^2 = '+sr2};
                            end
                            if ltfittype == 2
                                annot(length(annot)+1) = {'τ_1 = '+lt1+' ps, τ_2 = '+...
                                    lt2+' ps, A_{1}/A_{2} = '+slt1lt2+', r^2 = '+sr2};
                            end
                            if ltfittype == 3
                                annot(length(annot)+1) = {'τ_{mono} = '+lt1+...
                                    ' ps, β = '+string(fitval(12))+', r^2 = '+sr2};
                            end
                            if ltfittype > 0 && r2 < 0.7
                                annot(end) = {'Poor fit, r^2 = ' + sr2};
                            end
                        end
                        
                        % Prepare output filename(s)
                        if filetype > 1 && pairedDCAC > 0
                            outname = string(F(1:end-8)) + "_DCAC";
                        else
                            outname = string(F(1:end-4));
                        end
                
                        if ltfittype > 0 % prepare figure legend
                            legend1 = {'Data','Baseline data','Baseline','Fit'};
                        else
                            legend1 = {'Data','Baseline data','Baseline'};
                        end
                        if pairedDCAC > 0 % make figure
                            corrDC = DC - DCbasecurve;
                            [DCpeak,DCmaxindex] = max(corrDC);
                            DCpeaksd = DCsds(DCmaxindex);
                            integDC = sum(corrDC);
                            DCintegsd = mean(DCsds,'all');
                            corrDCbase = [corrDC(2:ilow); corrDC(ihigh:end)];
                            DCnoise = range(corrDCbase);
                            DCsnr = DCpeak/DCnoise;
                            if pairedDCAC == 2
                                ylabelDC = 'SPCM signal (cpms)';
                            else
                                ylabelDC = 'DC signal (AU)';
                            end
                
                            % Plot DC and AC side-by-side
                            fig = figure('visible',visibility); 
                            if ltfittype > 0
                                tiledlayout(2,3);
                            else
                                tiledlayout(1,2);
                            end
                            nexttile([1 1]); 
                            hold on; errorbar(t,DC,DCsds,'o');
                            plot(tbase,DCbase,'o',t,DCbasecurve,'-'); hold off;
                            xlabel('Time (ps)'); ylabel(ylabelDC);
                            legend('Data','Baseline data','Baseline');
                            nexttile([1 1]); 
                            if ltfittype > 0
                                plot(tint,resid); xlabel('Time (ps)'); 
                                ylabel('Residuals (AU)'); legend('Residuals');
                                nexttile([1 1]);
                                plot(wn,powerfft); xlabel('Frequency (cm^{-1})'); xlim([0 Inf]);
                                ylabel('|FFT|^{2}'); legend('FFT');
                                nexttile([1 3]);
                            end
                            hold on; errorbar(t,signal,sigsds,'o');
                            plot(tbase,sigbase,'go',t,basecurve,'-'); 
                            if ltfittype > 0
                                plot(tint,fitcurve,'-')
                            end
                            hold off;
                            xlabel('Time (ps)'); ylabel('AC signal (AU)');
                            legend(legend1);
                            annotation('textbox',annotdim,'String',annot);
                            if writefigsyn == 1
                                saveas(fig,outname+'_proc.png')
                            end
                        else % plot only signal (AC)
                            fig = figure('visible',visibility);
                            if ltfittype > 0
                                tiledlayout(3,2); nexttile([1 1]);
                                plot(tint,resid); xlabel('Time (ps)'); 
                                ylabel('Residuals (AU)'); legend('Residuals');
                                nexttile([1 1]); plot(wn,powerfft); 
                                xlabel('Frequency (cm^{-1})'); xlim([0 Inf]);
                                ylabel('|FFT|^{2}'); legend('FFT');
                                nexttile([2 2]);
                            end
                            hold on; errorbar(t,signal,sigsds,'o')
                            plot(tbase,sigbase,'go',t,basecurve,'-'); 
                            if ltfittype > 0
                                plot(tint,fitcurve,'-')
                            end
                            hold off;
                            xlabel('Time (ps)'); ylabel('AC signal (AU)');
                            legend(legend1);
                            annotation('textbox',annotdim,'String',annot);
                            if writefigsyn == 1
                                saveas(fig,outname+'_proc.png')
                            end
                        end
    
                        if filetype ~= 3 % solution data -> write master array
                            master(1:length(t),indct,ii,jj) = t(:);
                            master(1:length(corrsig),indccorrsig,ii,jj) = corrsig(:);
                            master(indrIRWN,indcvalue,ii,jj) = IRWN;
                            master(indrprobe,indcvalue,ii,jj) = prbWL;
                            master(indrsigpeak,indcvalue,ii,jj) = sigpeak;
                            master(indrpeaksd,indcvalue,ii,jj) = peaksd;
                            master(indrintegsig,indcvalue,ii,jj) = integsig;
                            master(indrintegsd,indcvalue,ii,jj) = integsd;
                            master(indrsnr,indcvalue,ii,jj) = snr;
                            master(indrsbr,indcvalue,ii,jj) = sbr;
                            master(indrnoise,indcvalue,ii,jj) = noise;
                            master(indrbasesd,indcvalue,ii,jj) = basesd;
                            master(indrintsnr,indcvalue,ii,jj) = intsnr;
                            master(indrbleach,indcvalue,ii,jj) = bleachrate;
                            master(indrtlist,indcvalue,ii,jj) = length(t);
                            master(indrlifetime1,indcvalue,ii,jj) = lifetime1;
                            master(indrlifetime2,indcvalue,ii,jj) = lifetime2;
                            master(indrlt1lt2,indcvalue,ii,jj) = lt1lt2;
                            master(indrnfiles,indcvalue,ii,jj) = totalfilesC;
                            master(indralignlength,indcvalue,ii,jj) = length(talign);
                            master(indrIRpower,indcvalue,ii,jj) = IRpower;
                            master(indrprbpower,indcvalue,ii,jj) = prbpower;
                            master(1:length(signal),indcrawsig,ii,jj) = signal(:);
                            master(1:length(basecurve),indcsigbase,ii,jj) = basecurve(:);
                            master(1:length(snrsweep),indcsnrsweep,ii,jj) = snrsweep(:);
                            master(1:length(tint),indctint,ii,jj) = tint(:);
                            master(1:length(sigint),indcsigint,ii,jj) = sigint(:);
                            master(1:length(fitcurve),indcfitcurve,ii,jj) = fitcurve(:);
                            master(1:length(fitval),indcfitval,ii,jj) = fitval(:);
                            master(1:length(corrsig),indcncsig,ii,jj) = corrsig(:)/sigpeak;
                            master(1:length(talign),indctalign,ii,jj) = talign;
                            master(1:length(sigalign),indcsigalign,ii,jj) = sigalign;
                
                            if pairedDCAC > 0
                                master(1:length(DC),indcrawDC,ii,jj) = DC(:);
                                master(1:length(corrDC),indccorrDC,ii,jj) = corrDC(:);
                                master(1:length(DCbasecurve),indcDCbase,ii,jj) = DCbasecurve(:);
                                master(indrDCpeak,indcvalue,ii,jj) = DCpeak;
                                master(indrDCpeaksd,indcvalue,ii,jj) = DCpeaksd;
                                master(indrintegDC,indcvalue,ii,jj) = integDC;
                                master(indrDCintegsd,indcvalue,ii,jj) = DCintegsd;
                                master(indrDCsnr,indcvalue,ii,jj) = DCsnr;
                                master(indrDCsbr,indcvalue,ii,jj) = DCsbr;
                                master(indrDCnoise,indcvalue,ii,jj) = DCnoise;
                                master(indrDCbleach,indcvalue,ii,jj) = DCbleachrate;
                            end

                            % To keep importdata from trimming NaN rows
                            if isnan(master(end,end,ii,jj))
                                master(end,end,ii,jj) = 69420;
                            end
                            
                            if writeprocyn == 1 % write files
                                temp = master(:,:,ii,jj);
                                writematrix(temp,outname+'_proc.dat',...
                                    'FileType','text')
                            end
                        end
                    end
                end
                if writeprocyn == 1 & filetype == 3 % write processed image tif
                    Save_tiff(outname+'_proc.tif',lt1array,...
                        lt2array,ssresidarray,FWHMarray,r2array)
                end
                prbpower = []; % clear powers for next file
                IRpower = []; % (otherwise, normalization can compound)
                if emptytype == 1
                    ltfittype = []; % empty for next file
                end
            end
        end
    end
    
    if writeprocyn == 1 && filetype ~= 3 % write batch file
        writematrix(master,'batch_align_flattenedv4.md','FileType','text') % backup, not very useful
        writematrix(size(master),'master_size.yml','FileType','text')
    end
    if testrun == 0 % close background figures if not a test run
        close all;
    end
else % reload previously processed data
    master_size = importdata('master_size.yml');
    master = NaN(master_size(1),master_size(2),master_size(3),master_size(4));
    D = pwd; % get current directory
    S = dir(fullfile(D,'*')); % search current directory
    N = setdiff({S([S.isdir]).name},{'.','..'}); % subfolders of D
    for ii = folders
        T = dir(fullfile(D,N{ii},'*_proc.dat'));
        C = {T(~[T.isdir]).name}; % all proc files
        if isempty(T) == 0
            for jj = 1:numel(C)
                F = fullfile(D,N{ii},C{jj});
                master(:,:,ii,jj) = importdata(F);
            end
        end
    end
end


%% Post-batch analysis
clc; close all;
clearvars -except maste* ind*
figvis = 'off';

% -Master array is 4D (x,y,ii,jj); ii,jj are folder,file
% x = rows (individual t points); y = columns
% -To pull columns: squeeze(master(1:TL),[indexcolumn],[folder],1:NF))
% TL = master(indrtlist,indcvalue,[folder],1)
% NF = master(indrnfiles,indcvalue,[folder],1)
% (1:NF optional - can also just use :)
% -To pull specific values: squeeze(master(indr...,indcvalue,[folder],1:NF))
% -Use 'rmmissing' to remove NaN values
% -Colorbar default for contour maps:
%cb=colorbar;
%cb.Label.String='Corrected signal (AU)';
%cb.Label.Rotation=270; cb.Label.VerticalAlignment = "bottom";

% % Red-white-blue gradient for contour map
% fineness = 100;
% map1 = [linspace(0,1,fineness) ones(1,fineness-1)];
% map3 = flip(map1);
% map2 = min(map1,map3);
% map = [map1.' map2.' map3.'];
% figure; % contour plot for w
% contourf(conc,rad,w)
% colormap(map);
% cb = colorbar;

% Default contour code
% analyzed = 3:7;
% 
% time = NaN(length(analyzed),max(master(indrtlist,indcvalue,:,:),[],"all"));
% wIR = NaN(max(master(indrnfiles,indcvalue,:,:),[],"all"),length(analyzed));
% csig = NaN(height(wIR),width(time),length(analyzed));
% 
% for i=1:length(analyzed)
%     time(i,1:master(indrtlist,indcvalue,analyzed(i),1)) = squeeze(master(1:master(indrtlist,indcvalue,analyzed(i),1),indct,analyzed(i),1)).';
%     wIR(1:master(indrnfiles,indcvalue,analyzed(i),1),i) = squeeze(master(indrIRWN,indcvalue,analyzed(i),1:master(indrnfiles,indcvalue,analyzed(i),1)));
%     csig(1:master(indrnfiles,indcvalue,analyzed(i),1),1:master(indrtlist,indcvalue,analyzed(i),1),i) = squeeze(master(1:master(indrtlist,indcvalue,analyzed(i),1),indccorrsig,analyzed(i),1:master(indrnfiles,indcvalue,analyzed(i),1))).';
% end
% 
% blankrectangle = [0.78 0.78 0.4 0.4];
% 
% for i=1:length(analyzed)
%     t1 = time(i,:); w1 = wIR(:,i);
%     [x1,x2] = meshgrid(t1,w1);
%     figure;%('visible',figvis);
%     tiledlayout(4,4,'TileSpacing','compact','Padding','compact');
%     nexttile([1 3]); % sum vs time
%     plot(t1(1:master(indrtlist,indcvalue,analyzed(i),1)),sum(csig(1:master(indrnfiles,indcvalue,analyzed(i),1),1:master(indrtlist,indcvalue,analyzed(i),1),i),1)); 
%     xlim([min(t1) max(t1)]); 
%     xlabel('Time delay (ps)'); ylabel('Sum (AU)');
%     nexttile([3 3]);
%     contourf(x1,x2,csig(:,:,i)); cb=colorbar;
%     cb.Label.String='Corrected signal (AU)';
%     cb.Label.Rotation=270; cb.Label.VerticalAlignment = "bottom";
%     xlabel('Time delay (ps)'); ylabel('ω_{IR} (cm^{-1})');
%     xlim([min(t1) max(t1)]); 
%     ylim([min(w1) max(w1)]);
%     nexttile([1 1]); xticks([]); yticks([]); % blank square
%     annotation('rectangle',blankrectangle,'Color',[1 1 1],'FaceColor',[1 1 1]); % box to cover
%     nexttile([3 1]); % sum vs freq
%     plot(sum(csig(1:master(indrnfiles,indcvalue,analyzed(i),1),1:master(indrtlist,indcvalue,analyzed(i),1),i),2),w1(1:master(indrnfiles,indcvalue,analyzed(i),1))); 
%     ylim([min(w1) max(w1)]);
%     xlabel('Sum (AU)'); ylabel('ω_{IR} (cm^{-1})');
% end




%% Functions

function Save_tiff(name, A, B, C, D, E)
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

