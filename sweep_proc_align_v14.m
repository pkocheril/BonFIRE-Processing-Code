%% Batch processing and alignment of 2D temporal sweep data
%%% initially written from sweep_processv10
%%% v2 - added IR sweep processing, heaviside function fitting, multi-peak
% fitting, configurable number of beating frequencies, writing out powers,
% cleaned up lifetime logic, improved .txt handling, automated DFG/idler
% assignment
%%% v3 - padding signal for fitting, fitting negative peaks, SPCM + PMT
%%%%%%%%%%%%
%%% v4 - updated fitting based on fit_temp_sweepv1 (big overhaul!)
%%% v5 - added pulse width to summary on PNGs, added import from .raw
%%% v6 - improved default contour plotting, added targetfolders writeout
%%% v7 - bug fixes, combined additional info writeout into one file,
% automated file type recognition
%%% v8 - auto-sorted contours, auto-Tlist generation
%%% v9 - bug fixes
%%% v10 - updated fitting constraints, added lifetimes to default analysis,
% updated picoEmerald tuning compensation, added DFG tuning compensation,
% added per-folder pulse width specification, added loading from batch.dat,
% added support for Windows file system, made smoother fit curves, 
% improved info file support, added auto-baseline fit type selection
%%% v11 - housekeeping, added chopping to IR power normalization
%%% v12 - updated contours
%%%%%%%%%%%%
%%% v13 - fixed time delay spacing (big change!), updated default plots,
% added FFT to master array
%%% v14 - cleaned up power corrections, bug fixes, added SRR and manual
% fitting, added more default vectors to post-batch analysis

% Initialize
%cd '/Users/pkocheril/Documents/Caltech/Wei Lab/Data/2024_02_14/'
clear; clc; close all;

% Main configuration options
loadprevious = 0; % 0 = new analysis, 1 = load previous, [] = auto-detect
runtype = 0; % 0 = process all, 1 = test run with a few files,
% 2 = examine a single file
targetfolders = []; % indices of folders to process, [] = dialog
t0pos = []; % specify t0 position (mm), [] = autofind

% Additional configuration options
ltfittype = []; % [] = auto-choose, 0 = no fitting, 1 = Gaussian*monoexp,
% 2 = Gaussian*biexp, 3 = Gauss*stretchexp, 4 = Gauss*strbiexp
basefittype = []; % [] = auto-choose, 0 = no baseline fit, 1 = linear, 
% 2 = exponential, 3 = exponential+linear
writeprocyn = 1; % 1 = write batch processed files, 0 = not
writefigsyn = 1; % 1 = write figure files, 0 = not
powernormyn = 0; % 0 = no normalization, 1 = normalize by IR power,
% 2 = normalize by probe and IR powers, 3 = 2 + PMT gain correction
setpulsewidth = []; % define pulse width (ps) in fit, [] = float
floatbase = 0.1; % fraction to float baseline coeffs (0.1 -> +/- 10%)
cutlow = 0.05; % lower baseline fit cutoff, (default 5th %ile)
cuthigh = 0.8; % upper baseline fit cutoff, (default 95th %ile)
verbose = 2; % 2 = all info on figure, 1 = regular figure annotation, 
% 0 = no figure annotation
troubleshoot = 0; % 0 = default, 1 = show everything along the way

% Even more configuration options (largely obsolete)
tempmod = 0; % 1 = temporal modulation in place, 0 = no beamsplitters
ND = 0; % ND filter strength in probe path (auto-update if possible)
PMT = 1; % PMT gain (norm to 1, auto-update if possible)
prbpowerset = 300; % probe power in mW (auto-update if possible)
IRpowerset = 70; % IR power in mW (auto-update if possible)
normIRpower = 50; % IR power on-sample to normalize to (default is 50)
normprobepower = 1; % probe power on-sample to normalize to (default 1)
pinholeyn = 0; % 0 = bypassed, 1 = in place
trimlastT = 0; % how many points to remove from end of Tlist (default 0)
trimfirstT = 0; % points to remove from start of Tlist (default 0)
ltmin = [];  % define minimum lifetime (ps) in fit, [] = default (0.1)
ltmax = []; % define maximum lifetime (ps) in fit, [] = default (100)
tuningrate = []; % PicoEmerald signal t0 tuning rate (ps/nm), [] = auto
minprobe = []; maxprobe = []; % probe range for aligning tD, [] = auto-find
pairedDCAC = []; % [] = auto-assign, 2 = SPCM + PMT CH2, 
% 1 = PMT CH1 + CH2, 0 = only PMT CH2
filetype = []; % [] = auto-assign, 1 = .txt, 2 = .tif (solution), 
% 3 = .tif (image), 4 = .raw (solution), 5 = .raw (image)
datatype = []; % [] = auto-assign, 0 = solution, 1 = image
filenameconv = []; % [] = auto-assign, 0 = other,
% 1 = [pumpWL]-[pumppower]-[signalWL]-[IRpower]-[ND]-[PMTgain]-[modfreq]-[PMTBW]-[etc],
% 2 = [etc]_[size]_[idler]_[DFG]_[power]_[channel],
% 3 = [IRWN]_[probeWL]_[etc]_[size]_[idler]_[DFG]_[power]_[channel],
% 4 = [probeWL]_[IRWN]_[etc]_[size]_[idler]_[DFG]_[power]_[channel]
testfilesperfolder = 1; % number of files to test per folder (default 1)

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
            runtype = 0; % not a test run if loading previous analysis
        else
            error("Multiple YMLs found - can't load previous analysis")
        end
    end
end

if troubleshoot == 1
    runtype = 2;
end

if runtype == 2 % to examine a single file - choose folder
    [indx,~] = listdlg('PromptString',{'Select one folder.'},...
    'SelectionMode','single','ListString',N);
    targetfolders = indx;
end

if isempty(targetfolders) == 1 % dialog to select folders
    [indx,~] = listdlg('PromptString',{'Select folders to process.'},...
    'SelectionMode','multiple','ListString',N);
    targetfolders = indx;
end

if runtype > 0  % setup test run conditions
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
indcpowernormsig = indcrawsig+1; % 5
indcsigbase = indcpowernormsig+1; % 6
indcsnrsweep = indcsigbase+1; % 7
indctint = indcsnrsweep+1; % 8
indcsigint = indctint+1; % 9
indcfitcurve = indcsigint+1; % 10
indcfitval = indcfitcurve+1; % 11
indcncsig = indcfitval+1; % 12
indctalign = indcncsig+1; % 13
indcsigalign = indctalign+1; % 14
indcfftx = indcsigalign+1; % 15
indcffty = indcfftx+1; % 16
indcrawDC = indcffty+1; % 17
indccorrDC = indcrawDC+1; % 18
indcDCbase = indccorrDC+1; % 19
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
indrpulsewidth = indrprbpower+1; % 21
indrDCpeak = indrpulsewidth+1; % 22
indrDCpeaksd = indrDCpeak+1; % 23
indrintegDC = indrDCpeaksd+1; % 24
indrDCintegsd = indrintegDC+1; % 25
indrDCsnr = indrDCintegsd+1; % 26
indrDCsbr = indrDCsnr+1; % 27
indrDCnoise = indrDCsbr+1; % 28
indrDCbleach = indrDCnoise+1; % 29 - minimum value of maxTlength

if loadprevious == 0 % run new analysis
    for passes = 1:2
        for ii = folders % loop through folders
            % On both passes - examine files in folders
            if isempty(filetype) == 1 % figure out data type on a per-folder basis
                autotype = 1;
                pairedDCAC = [];
            else
                autotype = 0;
            end
    
            if autotype == 1 % figure out file type
                checkraw = dir(fullfile(D,N{ii},'*CH2.raw')); % check for raw files
                checkraw = {checkraw(~[checkraw.isdir]).name};
                if isempty(datatype) == 1
                    if isempty(checkraw)
                        warning('No raw files found - data assumed to be a solution sample. Specify data type and re-run if needed.');
                        datatype = 0;
                    else
                        lookforFOV = strfind(checkraw(1),'FOV'); % check first file in folder
                        if isempty(lookforFOV) == 0 % file name doesn't have 'FOV'
                            datatype = 0; % solution data
                        else % file name has 'FOV'
                            datatype = 1; % image data
                        end
                    end
                end
                if numel(checkraw) > 0 % raw files present
                    if datatype == 1 % image data
                        filetype = 5;
                    else
                        filetype = 4; % solution data
                    end
                    if isempty(pairedDCAC) == 1
                        CH1raw = dir(fullfile(D,N{ii},'*CH1.raw')); % check for CH1 files
                        CH1raw = {CH1raw(~[CH1raw.isdir]).name};
                        SPCMraw = dir(fullfile(D,N{ii},'*SPCM.raw')); % check for SPCM files
                        SPCMraw = {SPCMraw(~[SPCMraw.isdir]).name};
                        if numel(CH1raw) == numel(checkraw)
                            pairedDCAC = 1; % paired CH1 and CH2 data
                        else
                            if numel(SPCMraw) == numel(checkraw)
                                pairedDCAC = 2;
                            else
                                pairedDCAC = 0;
                            end
                        end
                    end
                else % no raw files
                    checktif = dir(fullfile(D,N{ii},'*CH2.tif')); % check for tif files
                    checktif = {checktif(~[checktif.isdir]).name};
                    if numel(checktif) > 0 % tif files present
                        if datatype == 1 % image data
                            filetype = 3;
                        else
                            filetype = 2; % solution data
                        end
                        if isempty(pairedDCAC) == 1
                            CH1tif = dir(fullfile(D,N{ii},'*CH1.tif')); % check for CH1 files
                            CH1tif = {CH1tif(~[CH1tif.isdir]).name};
                            SPCMtif = dir(fullfile(D,N{ii},'*SPCM.tif')); % check for SPCM files
                            SPCMtif = {SPCMtif(~[SPCMtif.isdir]).name};
                            if numel(CH1tif) == numel(checktif)
                                pairedDCAC = 1; % paired CH1 and CH2 data
                            else
                                if numel(SPCMraw) == numel(checkraw)
                                    pairedDCAC = 2;
                                else
                                    pairedDCAC = 0;
                                end
                            end
                        end
                    else % no tif or raw files
                        checktxt = dir(fullfile(D,N{ii},'*.txt'));
                        checktxt = {checktxt(~[checktxt.isdir]).name};
                        if numel(checktxt) == 0
                            error('No data files found');
                        else % txt files present
                            filetype = 1;
                        end
                    end
                end
            end
            if datatype == 1 % for image data
                writefigsyn = 0; % don't write individual figures
            end
            if passes == 1 % 1st pass: pre-loop to figure out dimensions of master array
                maxTlength = 0; maxfilesinfolder = 0; % counter variables
                if filetype == 1  % look for .txt files in subfolders
                    T = dir(fullfile(D,N{ii},'*.txt'));
                else
                    if filetype < 4 % look for CH2 tifs
                        T = dir(fullfile(D,N{ii},'*CH2.tif'));
                    else % look for CH2 raws
                        T = dir(fullfile(D,N{ii},'*CH2.raw'));
                    end
                end
                C = {T(~[T.isdir]).name}; % all data files
                if numel(C) > maxfilesinfolder
                    maxfilesinfolder = numel(C);
                end
                if isempty(T) == 0 % skips subfolders without data files
                    for jj = 1:numel(C) % jj = file number in subfolder ii
                        F = fullfile(D,N{ii},C{jj}); % current file
                        if filetype == 1 % load Tlist from .txt
                            data = importdata(F); % load data
                            x = data(:,1); % save delay position as x
                        else % load Tlist from Tlist.txt
                            Tlistname = fullfile(D,N{ii},'Tlist.txt'); % current Tlist
                            if isfile(Tlistname) % check if Tlist exists
                                x = importdata(Tlistname); % load Tlist
                                trueTlist = 1;
                            else % guess a Tlist
                                trueTlist = 0;
                                x = linspace(176,180.05,82);
                            end
                        end
                        if x(1) == x(end) && trimlastT == 0 % if x isn't single-valued
                            trimlastT = 1; % trimming is needed
                        end
                        x = x(trimfirstT+1:end-trimlastT); % point trimming
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
                if autotype == 1
                    filetype = [];
                end

                if 10*maxTlength-1 < indrDCbleach % make sure master array is large enough
                    maxTlength = indrDCbleach;
                else
                    maxTlength = 10*maxTlength-1; % interpolate for tD compensation
                end
                
                if runtype == 2 % to examine a single file - choose file
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
                    if tempmod == 1 % pre-load beamsplitter response curve if needed
                        BSdata = importdata('/Users/pkocheril/Documents/Caltech/Wei Lab/Spreadsheets/BP145B2/BP145B2.csv');
                    end
                end
            end

            if passes == 2 % 2nd pass: full analysis
                if filetype == 1  % look for .txt files in subfolders
                    T = dir(fullfile(D,N{ii},'*.txt'));
                    infofilesraw = [];
                    infofiles = [];
                else % look for CH2 images
                    infofilesraw = dir(fullfile(D,N{ii},'*IRsweep*.txt'));
                    if isempty(infofilesraw) == 1 % if no IRsweeps
                        infofilesraw = dir(fullfile(D,N{ii},'*XYZT*.txt')); % look for XYZT info
                        if isempty(infofilesraw) == 1
                            infofilesraw = dir(fullfile(D,N{ii},'*Z.txt')); % look for Z stack info
                            if isempty(infofilesraw) == 1
                                infofilesraw = dir(fullfile(D,N{ii},'*T.txt')); % look for T stack info
                            end
                        end
                    end
                    infofiles = {infofilesraw(~[infofilesraw.isdir]).name}; % info files
                    if filetype < 4 % tifs
                        T = dir(fullfile(D,N{ii},'*CH2.tif')); 
                    else % raws
                        T = dir(fullfile(D,N{ii},'*CH2.raw'));
                    end
                end
                C = {T(~[T.isdir]).name}; % all data files
                if numel(infofiles) >= numel(C) % can have extra infofiles during scanning
                    infofound = 1;
                else
                    infofound = 0;
                end
                if runtype > 0 % limit test run files per folder
                    totalfilesC = testfilesperfolder;
                else
                    totalfilesC = numel(C);
                end
                if runtype == 2 % load only the chosen file
                    totalfilesC = targetfile;
                    startjj = targetfile;
                else
                    startjj = 1;
                end
                if isempty(T) == 0 % skip subfolders without data files
                    for jj = startjj:totalfilesC % jj = file number in subfolder ii
                        F = fullfile(D,N{ii},C{jj}); % current file

                        % Empty parameters for new file
                        IRpower = []; prbpower = [];

                        if isempty(filenameconv) == 1
                            autoconv = 1;
                        else
                            autoconv = 0;
                        end
        
                        % Filename parsing
                        macseparator = count(F,'/');
                        winseparator = count(F,'\');
                        if macseparator > winseparator % MacOS/UNIX
                            delimiter = "/";
                        else % Windows
                            delimiter = "\";
                        end
                        folderinfo = split(F,delimiter); % split path into folders
                        
                        if autoconv == 1 % guess file naming convention
                            if filetype == 1 % .txt
                                lookforhyphens = count(folderinfo(end),'-');
                                if numel(lookforhyphens) > 5
                                    filenameconv = 1;
                                else
                                    filenameconv = 0;
                                end
                            else % .tif or .raw
                                filename = char(folderinfo(end));
                                if strcmp(filename(5:9),'cm-1_') % filename ####cm-1_...
                                    filenameconv = 3; % probe sweep
                                else
                                    filenameconv = 2; % not a probe sweep
                                end
                            end
                        end
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
                            DFGWN = 1e7/DFGWL;
                            idlerWN = 1e7*((1/1031.2)-(1/signalWL));
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
                        else % image file (.raw or .tif)
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
                            imgdim = char(fileinfo(end-4)); % '5X5'
                            imgdimsplit = split(string(imgdim),"X");
                            xsteps = imgdimsplit(1); % "5"
                            ysteps = imgdimsplit(end); % "5"
                            xsteps = sscanf(xsteps,'%f');
                            ysteps = sscanf(ysteps,'%f');
                            prbWL = 1; % <-- unknown from filename
                            if IRpowermeter ~= 0
                                IRpower = IRpowermeter*300; % assuming 300 mW scale
                            end
                        end
                        if infofound == 1 % parse info from IRsweep.txt files
                            infofilename = fullfile(D,N{ii},infofiles{jj});
                            infoarray = importdata(infofilename);
                            infofile = infoarray.data;
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
                                    PMT = infofile(22); % PMT gain
                                    PMTBW = 1e3*infofile(23); % PMT bandwidth (kHz in file)
                                    modfreq = 1e3*infofile(24); % IR modulation frequency (kHz in file)
                                    pinholeyn = infofile(25); % 0 = bypassed, 1 = in place
                                    prbdichroic = infofile(26);
                                    pmtcmosmirror = infofile(27);
                                    pmtfilter = infofile(28);
                                    spcmfilter = infofile(29);
                                    cmosfilter = infofile(30);
                                end
                            else % if not extended info - guess additional parameters
                                if verbose == 2
                                    verbose = 1; % can't print out experimental parameters
                                end
                            end
                        end
                        % Assign IRWN by idler wavelength
                        idlerWL = 1e7/idlerWN;
                        signalWL = 1/((1/1031.2)-(1/idlerWL));
                        if idlerWL < 2600 % nm
                            IRWN = DFGWN;
                        else
                            IRWN = idlerWN;
                        end
                        if filenameconv == 3 % probe sweep
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
                        end
                        if filenameconv == 4 % IR sweep (deprecated)
                            IRWNstrfull = char(fileinfo(2)); % IRWN - "1000cm-1 (10)"
                            IRWNsplit = split(IRWNstrfull," "); % "1000cm-1" "(10)"
                            IRWN1char = char(IRWNsplit(1)); % '1000cm-1'
                            IRWN1chartrim = IRWN1char(1:end-4); % "1000"
                            IRWN1 = sscanf(IRWN1chartrim,'%f'); % 1000
                            IRWN2char = char(IRWNsplit(2)); % "(10)"
                            IRWN2chartrim = IRWN2char(2:end-1); % "10"
                            IRWN2 = sscanf(IRWN2chartrim,'%f'); % 10
                            IRWN = IRWN1+10*IRWN2; % 1100
                        end
                        if isempty(IRpower) == 1 % if no IR power found
                            IRpower = IRpowerset; % set to power from config
                        end
                        if isempty(prbpower) == 1 % if no probe power found
                            prbpower = prbpowerset; % set to power from config
                        end

                        % Power corrections
                        if DFGWN == IRWN
                            IRpower = IRpower/2; % lose 50% power to chopper (DFG only)
                        end
                        if pinholeyn == 1
                            prbpower = 0.27*prbpower;
                        end
                        prbpower = prbpower/(10^(ND));
            
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
                        else % load Tlist and image
                            Tlistname = fullfile(D,N{ii},'Tlist.txt'); % current Tlist
                            if isfile(Tlistname) % check if Tlist exists
                                x = importdata(Tlistname); % load Tlist
                                trueTlist = 1;
                            else % dialog to generate a Tlist
                                trueTlist = 0;
                                prompt = {'Specify start position (mm).','Specify end position (mm).','Specify step spacing (mm).'};
                                dlgtitle = 'Tlist Not Found';
                                dims = [1 45; 1 45; 1 45];
                                definput = {'176','180.05','0.05'};
                                Tlistanswer = inputdlg(prompt,dlgtitle,dims,definput);
                                Tlistanswer = str2double(Tlistanswer);
                                numberofpoints = 1+(Tlistanswer(2)-Tlistanswer(1))/Tlistanswer(3);
                                x = linspace(Tlistanswer(1),Tlistanswer(2),numberofpoints);
                            end
                            if filetype < 4 % .tif
                                data = double(tiffreadVolume(F));
                            else % .raw
                                imgid = fopen(F);
                                imgvector = fread(imgid,'real*8'); % load raw file as a vector
                                matrixarea = length(imgvector)/length(x);
                                if matrixarea == xsteps*ysteps
                                    %'yes';
                                else
                                    if mod(sqrt(matrixarea),1) == 0
                                        xsteps = sqrt(matrixarea); ysteps = sqrt(matrixarea);
                                    else
                                        if filetype == 3 || filetype == 5
                                            if trueTlist == 1
                                                warning('Mismatch between dimensions of .raw file and filename.')
                                            else
                                                warning('Mismatch in dimensions of file - check Tlist specification and file dimensions.')
                                            end
                                        end
                                    end
                                end
                                data = NaN(xsteps,ysteps,length(x));
                                for i=1:length(x) % convert raw vector to array
                                    startpoint = xsteps*ysteps*i-(xsteps*ysteps-1); 
                                    endpoint = xsteps*ysteps*i;
                                    imgvectpiece = imgvector(startpoint:endpoint);
                                    for j=1:ysteps
                                        tempstart = xsteps*j-(xsteps-1);
                                        tempend = xsteps*j;
                                        data(j,:,i) = imgvectpiece(tempstart:tempend);
                                    end
                                end
                            end
                            imagesize = size(data);
                            if pairedDCAC > 0
                                if pairedDCAC == 2 % SPCM + CH2
                                    F_DC = F;
                                    if filetype < 4 % .tif
                                        F_DC(end-6:end+1) = "SPCM.tif";
                                        DCdata = double(tiffreadVolume(F_DC));
                                    else % .raw
                                        F_DC(end-6:end+1) = "SPCM.raw";
                                        DCid = fopen(F_DC);
                                        DCimgvector = fread(DCid,'real*8'); % load raw file as a vector
                                        DCdata = NaN(xsteps,ysteps,length(x));
                                        for i=1:length(x) % convert raw vector to array
                                            DCstartpoint = xsteps*ysteps*i-(xsteps*ysteps-1); 
                                            DCendpoint = xsteps*ysteps*i;
                                            DCimgvectpiece = imgvector(DCstartpoint:DCendpoint);
                                            for j=1:ysteps
                                                DCtempstart = xsteps*j-(xsteps-1);
                                                DCtempend = xsteps*j;
                                                DCdata(j,:,i) = DCimgvectpiece(DCtempstart:DCtempend);
                                            end
                                        end
                                    end
                                else % pairedDCAC = 1, PMT CH1 + CH2
                                    F_DC = F;
                                    F_DC(end-4) = "1";
                                    if filetype < 4 % .tif
                                        DCdata = double(tiffreadVolume(F_DC));
                                    else % .raw
                                        DCid = fopen(F_DC);
                                        DCimgvector = fread(DCid,'real*8'); % load raw file as a vector
                                        DCdata = NaN(xsteps,ysteps,length(x));
                                        for i=1:length(x) % convert raw vector to array
                                            DCstartpoint = xsteps*ysteps*i-(xsteps*ysteps-1); 
                                            DCendpoint = xsteps*ysteps*i;
                                            DCimgvectpiece = DCimgvector(DCstartpoint:DCendpoint);
                                            for j=1:ysteps
                                                DCtempstart = xsteps*j-(xsteps-1);
                                                DCtempend = xsteps*j;
                                                DCdata(j,:,i) = DCimgvectpiece(DCtempstart:DCtempend);
                                            end
                                        end
                                    end
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
                            if filetype == 2 || filetype == 4 % avg & stdev of XY data (|| = OR)
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
                        
                        if runtype > 0 % do test run in the middle of an image
                            startiii = round(imagesize(1)/2);
                            startjjj = round(imagesize(2)/2);
                            imagesize(1) = (2-1)+startiii; 
                            imagesize(2) = (1-1)+startjjj;
                        else
                            startiii = 1; startjjj = 1;
                        end
        
                        if filetype == 1 || filetype == 2 || filetype == 4 % solution data
                            imagesize(1) = 1; imagesize(2) = 1;
                            startiii = 1; startjjj = 1;
                        end
            
                        for iii=startiii:imagesize(1)
                            for jjj=startjjj:imagesize(2)
                                if filetype == 3 || filetype == 5 % load pixel data
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
                                    indexoffset = floor(length(x)/10); % using an offset to ignore sharp decay at start
                                    [sigmax, t0index] = max(signal(indexoffset:end)); % estimating t0 by peak
                                    t(:) = (x(:)-x(t0index+indexoffset-1))*4/cmmps; % time vector in ps
                                else % defined t0 position
                                    t(:) = (x(:)-t0pos)*4/cmmps; % time vector in ps
                                end

                                nonnormsig = signal;
                                if powernormyn > 0 % power normalization
                                    signal = signal./IRpower;
                                    sigsds = sigsds./IRpower;
                                    DC = DC./IRpower;
                                    DCsds = DCsds./IRpower;
        
                                    if powernormyn > 1 % 2 = IR and probe
                                        if tempmod == 1 % correct for temporal modulation
                                            transmittance = BSdata.data(prbWL,1)/100; % picoEmerald is p-polarized (horizontal)
                                            reflectance = BSdata.data(prbWL,4)/100; % CSV in %
                                            prbpower = transmittance*reflectance*prbpower; % 62.5 mW -> 15.6 mW
                                        end
                                        signal = signal./prbpower;
                                        sigsds = sigsds./prbpower;
                                        DC = DC./prbpower;
                                        DCsds = DCsds./prbpower;
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
                                if isempty(basefittype) == 1
                                    emptybase = 1;
                                    frontmean = mean(signal(1:ilow));
                                    tailmean = mean(signal(ihigh:end));
                                    fronttaildiff = 2*(frontmean-tailmean)/(frontmean+tailmean);
                                    if fronttaildiff > 0.03
                                        basefittype = 3;
                                    else
                                        basefittype = 1;
                                    end
                                else
                                    emptybase = 0;
                                end
                                if basefittype >= 2 % if exponential baseline
                                    startbase = 1; % start at first point
                                else % if linear baseline
                                    startbase = 2; % start at second point
                                end
                                tbase = [t(startbase:ilow); t(ihigh:end)]; % trimmed t
                                sigbase = [signal(startbase:ilow); signal(ihigh:end)]; % trimmed signal
                                sbr = max(signal(ilow:ihigh))/mean(signal(ihigh:end));
                                if pairedDCAC > 0
                                    DCbase = [DC(startbase:ilow); DC(ihigh:end)]; % trimmed DC
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
                                    lbbase = [-Inf -Inf 0 -Inf -Inf];
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
                                    if troubleshoot == 0
                                        baseopt.Display = 'off'; % silence console output
                                    end
                                    basegs = [min(sigbase) 0.1*(max(sigbase)-min(sigbase)) 0.1 0 0];
                                    basefit = lsqnonlin(basefn,basegs,lbbase,ubbase,baseopt);
                                    basecurve = basefit(1)+basefit(2)*exp(-basefit(3)*(t-basefit(4)))+...
                                        basefit(5)*t;
                                    if basefit(2) > 0 && basefit(3) > 0 
                                        bleachrate = basefit(3); % ps-1
                                    else % bleach rate doesn't have physical meaning
                                        bleachrate = 0;
                                    end
                                    if pairedDCAC > 0 % fit DC baseline - always exp+linear
                                        DCbasefn = @(b) b(1)+b(2)*exp(-b(3)*(tbase-b(4)))+b(5)*tbase-DCbase;
                                        DCbgs = [min(DCbase) 0.1*(max(DCbase)-min(DCbase)) 0.1 0 0];
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
                                if pairedDCAC > 0
                                    corrDC = DC - DCbasecurve;
                                    [DCpeak,DCmaxindex] = max(corrDC);
                                    DCpeaksd = DCsds(DCmaxindex);
                                    integDC = sum(corrDC);
                                    DCintegsd = mean(DCsds,'all');
                                    corrDCbase = [corrDC(2:ilow); corrDC(ihigh:end)];
                                    DCnoise = range(corrDCbase);
                                    DCsnr = DCpeak/DCnoise;
                                end

                                % Calculate t spacing
                                deltat = zeros(length(t)-1,1);
                                for i=1:length(t)-1
                                    deltat(i) = round(1e2*(t(i+1)-t(i)))/1e2; % round to nearest 0.01 ps (had issues without rounding)
                                end
                                tspacing = unique(deltat); % length 1 <-> evenly spaced t
                                
                                if isempty(ltfittype) == 1 % auto-choose lifetime fitting
                                    emptytype = 1;
                                    if IRWN > 1680 && IRWN < 2400 % triple bond or carbonyl
                                        ltfittype = 1; % single exponential
                                    else
                                        % if IRWN > 1293 && IRWN < 1315 && prbWL > 699 && prbWL < 760
                                        %     ltfittype = 4; % biexp (exp + strexp)
                                        % else
                                        %     ltfittype = 2; % biexponential
                                        % end
                                        ltfittype = 2; % biexponential
                                    end
                                else
                                    emptytype = 0;
                                end
                        
                                % Lifetime fitting
                                if ltfittype > 0
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
                        
                                    if max(tint) < 100 % auto-pad signal if tint doesn't go past 100 ps
                                        padsignal = 1;
                                    else
                                        padsignal = 0;
                                    end

                                    % Prepare for signal padding
                                    tintspace = tint(2)-tint(1);
                                    originallength = length(tint); % should be an odd number
                                    if mod(originallength-1,4) == 0
                                        padlength = round(originallength/2)+1;
                                    else
                                        padlength = round(originallength/2);
                                    end
                                    tlong = linspace(tint(1),tint(end)+tintspace*padlength,length(tint)+padlength).';
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

                                    % Prepare convolution input (same spacing, half length, centered around zero)
                                    tc = linspace(0.5*min(tint),0.5*max(tint),0.5*(length(tint)+1)).';
                                    if powernormyn > 0
                                        IRFamp = peaksign*IRpower/1e3; % IRF amplitude (~ IR power in W)
                                    else
                                        IRFamp = peaksign*50/1e3; % default to 50 mW
                                    end
                        
                                    fun = @(r) r(1)+r(2)*exp(-r(3)*(tint-r(4)))+r(5)*tint+ ... % baseline
                                        conv(exp(-((tc-r(6))./(r(7)/(2*sqrt(log(2))))).^2), ...
                                        IRFamp*heaviside(tc).*(r(8)*exp(-(tc)./r(9))+r(10)*exp(-((tc)./(r(11))).^r(12)))) - sigint;
                                    
                                    % Initial guesses
                                    r0 = [basefit,0,2.5,... % basecoefs, IRF center (ps), IRF width (ps)
                                        0.02,1,0.02,10,1]; % amp 1, τ1 (ps), amp 2, τ2 (ps), β (stretch)
                        
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
                                    lb = [newlbbase,-10,1*sqrt(2),... % basecoefs, IRF center (ps), IRF min (ps)
                                        0,0.1,... % amp 1, τ1min (ps)
                                        0,0.1,... % amp 2, τ2min (ps)
                                        1e-2]; % β min
                                    ub = [newubbase,20,3*sqrt(2),... % basecoefs, IRF center (ps), IRF max (ps)
                                        Inf,100,... % amp 1, τ1max (ps)
                                        Inf,100,... % amp 2, τ2max (ps)
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
                                        if length(setpulsewidth) == 1 % single pulse width
                                            lb(7) = setpulsewidth*sqrt(2); ub(7) = setpulsewidth*sqrt(2);
                                        else % multiple pulse widths
                                            if length(setpulsewidth) == numel(N) % if one per folder
                                                lb(7) = setpulsewidth(ii)*sqrt(2); ub(7) = setpulsewidth(ii)*sqrt(2);
                                            else
                                                lb(7) = setpulsewidth(1)*sqrt(2); ub(7) = setpulsewidth(1)*sqrt(2);
                                                warning('Number of specified pulse widths does not match total number of folders - using first specified only.')
                                            end
                                        end
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
                                    if troubleshoot == 0
                                        options.Display = 'off'; % silence console output
                                    end
                                    fitval = lsqnonlin(fun,r0,lb,ub,options); % fit coeffs

                                    fitvector = fitval(1)+fitval(2)*exp(-fitval(3)*(tint-fitval(4)))+fitval(5)*tint+ ... % baseline
                                        conv(exp(-((tc-fitval(6))./(fitval(7)/(2*sqrt(log(2))))).^2), ...
                                        IRFamp*heaviside(tc).*(fitval(8)*exp(-(tc)./fitval(9))+fitval(10)*exp(-((tc)./(fitval(11))).^fitval(12))));

                                    % Brute-forced smoother fit curve
                                    tfit = linspace(tint(1),tint(originallength),originallength*10-9).'; % smoother t vector for 
                                    fitcurve = spline(tint,fitvector,tfit);

                                    if troubleshoot == 1
                                        figure; tiledlayout(2,1); nexttile([1 1]); 
                                        plot(tint,sigint,'o',tint,fitvector,'-'); 
                                        nexttile([1 1]); plot(tint,sigint-fitvector); title('Fitting without crop')
                                    end

                                    % Trim padded signal to original length
                                    tint = tint(1:originallength);
                                    sigint = sigint(1:originallength);
                                    fitvector = fitvector(1:originallength);
                                    
                                    % Calculate pulse width from fit
                                    IRFwidth = fitval(7); % ps
                                    FWHM = IRFwidth/sqrt(2); % ps
        
                                    % Biexponential lifetime assignments
                                    lifetime1 = min(fitval(9),fitval(11)); % ps
                                    lifetime2 = max(fitval(9),fitval(11)); % ps
                                    % Lifetime amplitude ratio
                                    if abs(fitval(9)) < abs(fitval(11)) % if fitval(9) is shorter (τ1)
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
                                    resid = sigint-fitvector;
                                    srr = sigpeak./std(resid(2:end));
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
                                        figure; tiledlayout(3,1); nexttile([1 1]); plot(tint,resid); nexttile([2 1]); plot(t,signal,'o',tint,fitvector); title('Fitting')
                                        % Check FFT
                                        figure; plot(wn,powerfft); title('FFT of residuals');
                                    end
                                else % no lifetime fitting
                                    tint = t; sigint = signal; fitvector = signal;
                                    resid = signal; wn = t; powerfft = signal;
                                    lifetime1 = 0; lifetime2 = 0; lt1lt2 = 0;
                                    fitval = zeros(12,1); r2 = 0; FWHM = 0;
                                end
        
                                if prbWL < 961 && prbWL > 699 % interpolate for tD compensation
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
                                lt1 = string(sprintf('%0.3g',lifetime1)); % round to 3 sig figs
                                lt2 = string(sprintf('%0.3g',lifetime2));
                                sr2 = string(sprintf('%0.5g',r2));
                                slt1lt2 = string(sprintf('%0.3g',lt1lt2));
                                pulsewidth = string(sprintf('%0.3g',FWHM));
                                
                                if verbose > 0 % make figure annotation
                                    annot = {'ω_{IR}/2πc = '+string(IRWN)+' cm^{-1}, λ_{probe} = '+...
                                        string(prbWL)+' nm, SNR = '+string(snr)};
                                    if ltfittype == 1
                                        annot(length(annot)+1) = {'SRR = '+string(srr)+', τ_{mono} = '+lt1+...
                                            ' ps, τ_{p} = '+pulsewidth+' ps, r^2 = '+sr2};
                                    end
                                    if ltfittype == 2
                                        annot(length(annot)+1) = {'SRR = '+string(srr)+', τ_1 = '+lt1+' ps, τ_2 = '+...
                                            lt2+' ps, τ_{p} = '+pulsewidth+...
                                            ' ps, A_{1}/A_{2} = '+slt1lt2+', r^2 = '+sr2};
                                    end
                                    if ltfittype == 3
                                        annot(length(annot)+1) = {'SRR = '+string(srr)+', τ_{mono} = '+lt1+...
                                            ' ps, τ_{p} = '+pulsewidth+' ps, β = '+...
                                            string(fitval(12))+', r^2 = '+sr2};
                                    end
                                    if ltfittype == 4
                                        annot(length(annot)+1) = {'SRR = '+string(srr)+', τ_1 = '+lt1+' ps, τ_2 = '+...
                                            lt2+' ps, τ_{p} = '+pulsewidth+' ps, β = '+string(fitval(12))+...
                                            ', A_{1}/A_{2} = '+slt1lt2+', r^2 = '+sr2};
                                    end
                                    if ltfittype > 0 && r2 < 0.7
                                        annot(end) = {'Poor fit, r^2 = ' + sr2};
                                    end
                                    if verbose == 2 % print out all experimental parameters available
                                        if pinholeyn == 1
                                            phlabel = 'pinhole in place';
                                        else
                                            phlabel = 'bypassed pinhole';
                                        end
                                        annot(length(annot)+1) = {'PMT gain '+string(PMT)+...
                                            ', PMT BW '+string(PMTBW/1e3)+' kHz, '+...
                                            phlabel+', DM'+string(prbdichroic)+...
                                            ', PMT BP'+string(pmtfilter)};
                                    end
                                else
                                    annot = {};
                                end
                                
                                % Prepare output filename(s)
                                if filetype > 1 && pairedDCAC > 0
                                    outname = string(F(1:end-8)) + "_DCAC";
                                else
                                    outname = string(F(1:end-4));
                                end

                                % Make labels for figure
                                basestring = strings(4,1);
                                basestring(1) = 'Baseline'; basestring(2) = 'Linear baseline';
                                basestring(3) = 'Exponential baseline'; basestring(4) = 'Exp + linear baseline';
                                dclegend = {'Data','Baseline data','Exp + linear baseline'};
                                aclegend = {'Data','Baseline data',basestring(basefittype+1)};
                                if ltfittype > 0
                                    aclegend(4) = {'Fit'};
                                end
                                if powernormyn > 0 % prepare y-axis normalization label
                                    normlabel = 'Norm. ';
                                    if powernormyn == 1
                                        pmtunitlabel = '(V/mW)';
                                        spcmunitlabel = '(cpms/mW)';
                                    else % powernormyn == 2
                                        pmtunitlabel = '(V/mW^{2})';
                                        spcmunitlabel = '(cpms/mW^{2})';
                                    end
                                else
                                    normlabel = 'Raw ';
                                    pmtunitlabel = '(V)';
                                end
                                if pairedDCAC > 0
                                    if pairedDCAC == 2
                                        ylabelDC = string(normlabel)+'SPCM signal '+string(spcmunitlabel);
                                    else
                                        ylabelDC = string(normlabel)+'DC signal '+string(pmtunitlabel);
                                    end
                                end
                                ylabelAC = string(normlabel)+'AC signal '+string(pmtunitlabel);

                                % Make colors for figure
                                dccolor = [0.5 0.5 0.5];
                                dcbasecolor = [0.9 0.1 0.1];
                                dcfitcolor1 = [0.7 0.4 0.1];
                                rescolor = [0.5 0.5 0.5];
                                accolor = [0.5 0.5 0.5];
                                acbasecolor = [0.2 0.8 0.8];
                                acfitcolor1 = [0.2 0.8 0.5];
                                acfitcolor2 = [0.2 0.2 0.8];

                                % Make figure
                                fig = figure('visible',visibility);
                                tiledlayout(3,2+(pairedDCAC > 0)); nexttile([1 2]); % add +2*(pairedDCAC > 2) for 4-channel plotting
                                if pairedDCAC > 0 % plot DC
                                    hold on; errorbar(t,DC,DCsds,'o','LineWidth',2,'Color',dccolor);
                                    plot(tbase,DCbase,'o','Color',dcbasecolor,'LineWidth',2);
                                    plot(t,DCbasecurve,'-','Color',dcfitcolor1,'LineWidth',2);
                                    xlabel('Time delay (ps)'); ylabel(ylabelDC); xlim([min(t) max(t)]);
                                    legend(dclegend); lgdc = legend; lgdc.EdgeColor = [1 1 1];
                                    ax = gca; ax.FontSize = 12; hold off;
                                    nexttile([1 1]); xticks([]); yticks([]);
                                    printoutdim = [0.67 0.7 0.3 0.25];
                                else
                                    xticks([]); yticks([]); printoutdim = [0.05 0.7 0.9 0.25];
                                end
                                % Plot AC
                                nexttile([2 2+(pairedDCAC > 0)]); hold on; xlim([min(t) max(t)]);
                                xlabel('Time delay (ps)'); ylabel(ylabelAC);
                                errorbar(t,signal,sigsds,'o','Color',accolor,'LineWidth',2);
                                plot(tbase,sigbase,'o','Color',acbasecolor,'LineWidth',2); 
                                plot(t,basecurve,'-','Color',acfitcolor1,'LineWidth',2); 
                                if ltfittype > 0 % plot fit
                                    plot(tfit,fitcurve,'-','Color',acfitcolor2,'LineWidth',2);
                                    yyaxis right; % plot residuals on secondary axis
                                    plot(tint,resid,'LineWidth',1); aclegend(5) = {'Residuals'};
                                    ylabel('Residuals (AU)'); legend('Residuals');
                                end
                                legend(aclegend); lgac = legend; lgac.EdgeColor = [1 1 1];
                                ax = gca; ax.FontSize = 12;
                                annotation('rectangle',printoutdim,'Color',[1 1 1],'FaceColor',[1 1 1]); % box to cover
                                annotation('textbox',printoutdim,'String',annot);
                                if writefigsyn == 1
                                    saveas(fig,outname+'_proc.png')
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
                                    master(indrpulsewidth,indcvalue,ii,jj) = FWHM;
                                    master(1:length(signal),indcpowernormsig,ii,jj) = signal(:);
                                    master(1:length(signal),indcrawsig,ii,jj) = nonnormsig(:);
                                    master(1:length(basecurve),indcsigbase,ii,jj) = basecurve(:);
                                    master(1:length(snrsweep),indcsnrsweep,ii,jj) = snrsweep(:);
                                    master(1:length(tint),indctint,ii,jj) = tint(:);
                                    master(1:length(sigint),indcsigint,ii,jj) = sigint(:);
                                    master(1:length(fitvector),indcfitcurve,ii,jj) = fitvector(:);
                                    master(1:length(fitval),indcfitval,ii,jj) = fitval(:);
                                    master(1:length(corrsig),indcncsig,ii,jj) = corrsig(:)/sigpeak;
                                    master(1:length(talign),indctalign,ii,jj) = talign;
                                    master(1:length(sigalign),indcsigalign,ii,jj) = sigalign;
                                    master(1:length(wn),indcfftx,ii,jj) = wn;
                                    master(1:length(powerfft),indcffty,ii,jj) = powerfft;
                        
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
                                        tempmaster = master(:,:,ii,jj);
                                        writematrix(tempmaster,outname+'_proc.dat',...
                                            'FileType','text')
                                    end
                                end
                            end
                        end
                        if writeprocyn == 1 && filetype == 3 % write processed image tif
                            Save_tiff(outname+'_proc.tif',lt1array,...
                                lt2array,ssresidarray,FWHMarray,r2array)
                        end
                        if emptytype == 1
                            ltfittype = []; % empty for next file
                        end
                        if emptybase == 1
                            basefittype = []; % empty for next file
                        end
                        if autoconv == 1
                            filenameconv = []; % empty for next file
                        end
                    end
                end
                if autotype == 1 && ii ~= folders(end)
                    filetype = []; % empty for next folder
                end
            end
        end
    end
    if writeprocyn == 1 && filetype ~= 3 % write batch file
        writematrix(master,'batch.dat','FileType','text') % backup, not very useful
        writematrix([size(master) targetfolders],'add_info.yml','FileType','text'); % additional info
    end
    if runtype == 0 % close background figures if not a test run
        close all;
    end
else % reload previously processed data
    addinfo = importdata('add_info.yml'); % load info from previous analysis
    master_size = addinfo(1:4);
    master = NaN(master_size(1),master_size(2),master_size(3),master_size(4));
    targetfolders = addinfo(5:end);
    D = pwd; % get current directory
    S = dir(fullfile(D,'*')); % search current directory
    N = setdiff({S([S.isdir]).name},{'.','..'}); % subfolders of D
    emptydatcount = 0;
    for ii = folders
        T = dir(fullfile(D,N{ii},'*_proc.dat'));
        C = {T(~[T.isdir]).name}; % all proc files
        if isempty(T) == 0
            for jj = 1:numel(C)
                F = fullfile(D,N{ii},C{jj});
                tempimport = importdata(F);
                master(1:height(tempimport),1:width(tempimport),ii,jj) = tempimport;
            end
        else
            emptydatcount = emptydatcount+1;
        end
    end
    if emptydatcount == length(folders) % no proc.dat files found
        checkdat = dir('batch.dat'); % try to load from batch.dat
        warning('Loading from batch.dat - master array may have errors.');
        if isempty(checkdat) == 0 % found batch.dat
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

% %% Manual fitting
% % Do a test run, then come here
% 
% % Prepare for signal padding
% tintspace = tint(2)-tint(1);
% originallength = length(tint); % should be an odd number
% if mod(originallength-1,4) == 0
%     padlength = round(originallength/2)+1;
% else
%     padlength = round(originallength/2);
% end
% tlong = linspace(tint(1),tint(end)+tintspace*padlength,length(tint)+padlength).';
% baselong = basefit(1)+basefit(2)*exp(-basefit(3)*(tlong-basefit(4)))+...
%     basefit(5)*tlong;
% if basefittype > 0
%     siglong = baselong;
%     siglong(1:length(sigint)) = sigint;
% else
%     if peaksign == 1
%         siglong(length(sigint)+1:length(tlong)) = min(sigint(end-5:end));
%     else
%         siglong(length(sigint)+1:length(tlong)) = max(sigint(end-5:end));
%     end
% end
% 
% % Signal padding for short-tailed data
% if padsignal == 1
%     tint = tlong;
%     sigint = siglong;
%     if troubleshoot == 1
%         figure; plot(tint,sigint,'-o',t,signal,'o'); title('Padding');
%     end
% end
% 
% tc = linspace(0.5*min(tint),0.5*max(tint),0.5*(length(tint)+1)).';
% fun = @(r) r(1)+r(2)*exp(-r(3)*(tint-r(4)))+r(5)*tint+ ... % baseline
%     conv(exp(-((tc-r(6))./(r(7)/(2*sqrt(log(2))))).^2), ...
%     IRFamp*heaviside(tc).*(r(8)*exp(-(tc)./r(9))+...
%     r(10)*exp(-((tc)./(r(11))).^r(12))...
%     )) - sigint;
% 
% % Initial guesses
% r0 = [basefit,0,2.5,... % basecoefs, IRF center (ps), IRF width (ps)
%     0.02,1,0.02,10,1,... % amp 1, τ1 (ps), amp 2, τ2 (ps), β (stretch)
%     ]; 
% 
% % Lower and upper bounds
% lb = [newlbbase,-10,1,... % basecoefs, IRF center (ps), IRF min (ps)
%     0,0.01,... % amp 1, τ1min (ps)
%     0,0.01,... % amp 2, τ2min (ps)
%     1e-2...  % β min
%     ];
% ub = [newubbase,20,3,... % basecoefs, IRF center (ps), IRF max (ps)
%     Inf,100,... % amp 1, τ1max (ps)
%     Inf,100,... % amp 2, τ2max (ps)
%     1e2...  % β max
%     ];
% 
% % Run lsqnonlin - minimizes error function by adjusting r
% options=optimoptions(@lsqnonlin);
% options.MaxFunctionEvaluations = 1e6; % let the fit run longer
% options.MaxIterations = 1e6; % let the fit run longer
% options.FunctionTolerance = 1e-12; % make the fit more accurate
% options.OptimalityTolerance = 1e-12; % make the fit more accurate
% options.StepTolerance = 1e-12; % make the fit more precise
% if troubleshoot == 0
%     options.Display = 'off'; % silence console output
% end
% fitval = lsqnonlin(fun,r0,lb,ub,options); % fit coeffs
% fitvector = fitval(1)+fitval(2)*exp(-fitval(3)*(tint-fitval(4)))+fitval(5)*tint+ ... % baseline
%     conv(exp(-((tc-fitval(6))./(fitval(7)/(2*sqrt(log(2))))).^2), ...
%     IRFamp*heaviside(tc).*(fitval(8)*exp(-(tc)./fitval(9))+...
%     fitval(10)*exp(-((tc)./(fitval(11))).^fitval(12))...
%     ));
% 
% fitval
% % Trim padded signal to original length
% tint = tint(1:originallength);
% sigint = sigint(1:originallength);
% fitvector = fitvector(1:originallength);
% 
% % Calculate residuals and ssresid
% resid = sigint-fitvector;
% srr = sigpeak./std(resid(2:end));
% ssresid = sum(resid.^2); % sum of squares of residuals
% tss = sum((sigint-mean(sigint)).^2); % total sum of squares
% r2 = 1-(ssresid/tss); % coefficient of determination (R^2)
% if r2 < 0 % remove nonsense R^2 values
%     r2 = 0;
% end
% 
% figure; tiledlayout(3,1); nexttile([2 1]); hold on;
% plot(tint,sigint,'o','LineWidth',2);
% plot(tint,fitvector,'-','LineWidth',2);
% ylabel('Signal (AU)');
% hold off; nexttile([1 1]); hold on;
% plot(tint,resid);
% xlabel('Time (ps)'); ylabel('Residuals (AU)')

%% Post-batch analysis
% -Master array is 4D (x,y,ii,jj); ii,jj are folder,file
% x = rows (individual t points); y = columns
% -To pull columns: squeeze(master(1:TL),[indexcolumn],[folder],1:NF))
% TL = master(indrtlist,indcvalue,[folder],1)
% NF = master(indrnfiles,indcvalue,[folder],1)
% (1:NF optional - can also just use :)
% -To pull specific values: squeeze(master(indr...,indcvalue,[folder],1:NF))
% -Use 'rmmissing' to remove NaN values

clc; close all;
clearvars -except targetfolders maste* ind*
figvis = 'on';
savecontours = 0; % 1 = save contours, 0 = not

% Default contour code
subset = targetfolders;

% Set up arrays
time = NaN(length(subset),max(master(indrtlist,indcvalue,:,:),[],"all"));
timealign = NaN(length(subset),max(master(indralignlength,indcvalue,:,:),[],"all"));
prb = NaN(max(master(indrnfiles,indcvalue,:,:),[],"all"),length(subset));
wIR = prb; lt1 = prb; lt2 = prb; ltr = prb; pkht = prb; pulsewidth = prb;
csig = NaN(height(wIR),width(time),length(subset));
rawsig = csig;
csiga = NaN(height(prb),width(timealign),length(subset));

% Pulling data
for i=1:length(subset)
    tlength = master(indrtlist,indcvalue,subset(i),1);
    alength = master(indralignlength,indcvalue,subset(i),1);
    nfiles = master(indrnfiles,indcvalue,subset(i),1);
    time(i,1:tlength) = squeeze(master(1:tlength,indct,subset(i),1)).';
    timealign(i,1:alength) = squeeze(master(1:alength,indctalign,subset(i),1)).';
    prb(1:nfiles,i) = squeeze(master(indrprobe,indcvalue,subset(i),1:nfiles));
    wIR(1:nfiles,i) = squeeze(master(indrIRWN,indcvalue,subset(i),1:nfiles));
    lt1(1:nfiles,i) = squeeze(master(indrlifetime1,indcvalue,subset(i),1:nfiles));
    lt2(1:nfiles,i) = squeeze(master(indrlifetime2,indcvalue,subset(i),1:nfiles));
    ltr(1:nfiles,i) = squeeze(master(indrlt1lt2,indcvalue,subset(i),1:nfiles));
    csig(1:nfiles,1:tlength,i) = squeeze(master(1:tlength,indccorrsig,subset(i),1:nfiles)).';
    rawsig(1:nfiles,1:tlength,i) = squeeze(master(1:tlength,indcrawsig,subset(i),1:nfiles)).';
    csiga(1:nfiles,1:alength,i) = squeeze(master(1:alength,indcsigalign,subset(i),1:nfiles)).';
    pkht(1:nfiles,i) = squeeze(master(indrsigpeak,indcvalue,subset(i),1:nfiles));
    pulsewidth(1:nfiles,i) = squeeze(master(indrpulsewidth,indcvalue,subset(i),1:nfiles));
end

ltavg = (ltr.*lt1+lt2)./(ltr+1);

% Plotting
for i=1:length(subset)
    if length(unique(rmmissing(prb(:,i)))) == 1 % no probe tuning --> IR sweep
        xstring = 'Time delay (ps)';
        ystring = 'ω_{IR}/2πc (cm^{-1})';
        w1 = wIR(:,i);
        w2 = prb(:,i);
        nonswept = string(round(w2(1)))+' nm';
        sumfreq = 1e7./w2+w1;
        t1 = time(i,:); 
        sigmatrix = csig(:,:,i);
        % Purple-white-green gradient for contour map
        startcolor = [0.5 0 0.5]; % purple
        endcolor = [0.2 0.8 0.2]; % green
        %contourlines = [min(sigmatrix,[],"all") 0 logspace(log10(min(abs(sigmatrix),[],"all")),log10(0.3*max(abs(sigmatrix),[],"all")),10) linspace(0.35*max(abs(sigmatrix),[],"all"),max(abs(sigmatrix),[],"all"),5)];
    else % probe sweep
        xstring = 'Compensated time delay (ps)';
        ystring = 'λ_{probe} (nm)';
        w1 = prb(:,i);
        w2 = wIR(:,i);
        nonswept = string(round(w2(1)))+' cm^{-1}';
        sumfreq = 1e7./w1+w2;
        t1 = timealign(i,:);
        sigmatrix = csiga(:,:,i);
        % Blue-white-orange gradient for contour map
        startcolor = [0 0 1]; % blue
        endcolor = [1 0.5 0]; % orange
        %contourlines = linspace(min(sigmatrix,[],"all"),max(sigmatrix,[],"all"),20);
    end

    if length(rmmissing(w1)) > 1
        % Auto-sorting
        [t1,tind] = sort(t1);
        [w1,wind] = sort(w1);
        [x1,x2] = meshgrid(t1,w1);
        sigmatrix = sigmatrix(wind,tind);
        logsigmatrix = log(abs(sigmatrix));

        % Color mapping
        map1 = [linspace(startcolor(1),1,100) linspace(1,endcolor(1),100)];
        map2 = [linspace(startcolor(2),1,100) linspace(1,endcolor(2),100)];
        map3 = [linspace(startcolor(3),1,100) linspace(1,endcolor(3),100)];
        map = [map1.' map2.' map3.'];

        % Plotting
        contour = figure('visible',figvis);
        tiledlayout(4,4,'TileSpacing','compact','Padding','compact');
        %tiledlayout(3,6,'TileSpacing','compact','Padding','compact');

        % Time 1D cross-section
        %nexttile([1 2]);
        nexttile([1 3]);
        %plot(t1,max(sigmatrix,[],1),'Color',startcolor,'LineWidth',2); ylabel('Peak (AU)');
        plot(t1,sum(sigmatrix(1:length(rmmissing(w1)),:)/length(rmmissing(w1)),1),'Color',startcolor,'LineWidth',2); ylabel('Mean (AU)');
        xlim([min(t1) max(t1)]); xlabel(xstring);

        % Normal contour
        %nexttile([2 2]);
        nexttile([3 3]);
        contourf(x1,x2,sigmatrix); colormap(map); cb=colorbar;
        cb.Label.String='Corrected signal at '+nonswept+' (AU)';
        cb.Label.Rotation=270; cb.Label.VerticalAlignment = "bottom";
        xlabel(xstring); ylabel(ystring);
        xlim([min(t1) max(t1)]); ylim([min(w1) max(w1)]);

        % Blank square (top-right)
        nexttile([1 1]);
        xticks([]); yticks([]);
        annotation('rectangle',[0.78 0.78 0.4 0.4],'Color',[1 1 1],'FaceColor',[1 1 1]); % box to cover

        % Frequency 1D cross-section
        %nexttile([2 2]); 
        nexttile([3 1]);
        plot(max(sigmatrix,[],2),w1,'Color',endcolor,'LineWidth',2); xlabel('Peak (AU)');
        %plot(sum(sigmatrix(:,1:length(rmmissing(t1))),2)/length(rmmissing(t1)),w1,'Color',endcolor,'LineWidth',2); xlabel('Mean (AU)');
        ylim([min(w1) max(w1)]); ylabel(ystring);
        if savecontours == 1
            saveas(contour,'Probe contours/'+nonswept+'_sweep.png')
        end
    end
end

% Default lifetime comparison
ltlegend = strings(length(subset),1);
figure('visible',figvis);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
% τ1
nexttile([1 2]); hold on;
for i=1:length(subset)
    redamt = exp((i-length(subset))*2/length(subset));
    greenamt = ((-4/length(subset)^2)*(i-0.5*length(subset))^2+1)^2;
    blueamt = exp(-2*i/length(subset));
    linecolor = [redamt greenamt blueamt];
    plot(sumfreq,lt1(:,i),'Color',linecolor,'LineWidth',2);
    ltlegend(i) = string(wIR(1,i))+' cm{-1}';
end
xticks([]); ylabel('τ_{1} (ps)');
hold off; xlim([min(sumfreq) max(sumfreq)]); ylim([0 5]); 
%xlim([13695 14747]); 
legend(ltlegend);
% τ2
nexttile([1 2]); hold on;
for i=1:length(subset)
    redamt = exp((i-length(subset))*2/length(subset));
    greenamt = ((-4/length(subset)^2)*(i-0.5*length(subset))^2+1)^2;
    blueamt = exp(-2*i/length(subset));
    linecolor = [redamt greenamt blueamt];
    plot(sumfreq,lt2(:,i),'Color',linecolor,'LineWidth',2);
end
xlabel('ω_{IR}+ω_{probe} (cm^{-1})'); ylabel('τ_{2} (ps)');
hold off; xlim([min(sumfreq) max(sumfreq)]); ylim([0 10]); 
%xlim([13695 14747]); ylim([5 8]);
legend(ltlegend);
% A1/A2
nexttile([2 1]); hold on;
for i=1:length(subset)
    redamt = exp((i-length(subset))*2/length(subset));
    greenamt = ((-4/length(subset)^2)*(i-0.5*length(subset))^2+1)^2;
    blueamt = exp(-2*i/length(subset));
    linecolor = [redamt greenamt blueamt];
    plot(sumfreq,ltr(:,i),'Color',linecolor,'LineWidth',2);
end
set(gca,'Yscale','log');
xlabel('ω_{IR}+ω_{probe} (cm^{-1})'); ylabel('A_{1}/A_{2}');
xlim([min(sumfreq) max(sumfreq)]); ylim([1e-2 1e2]); 
%xlim([13695 14747]); ylim([1e-2 1e2]); 
legend(ltlegend);


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

