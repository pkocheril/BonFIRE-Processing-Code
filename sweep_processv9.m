%%% Batch processing of temporal sweep data (written for 2D sweeps)
%%% initially written from batchfit_temporal_v4
%%% v2 - added writing to master array (4D data hypercube)
%%% v3 - added full filename parsing and lifetime fitting from 
% batchfit_temporalv4 (not floating exponential time center anymore),
% power normalization (with ND and BP145B2 beamsplitter curve), and
% unique file extensions for output files (still all text files)
%%% v4 - added image fitting from fit_temporalv8, cleaned up power
% normalization, fixed DFGWN parsing for filenameconv 2 & 3
%%% v5 - added parsing experimental parameters from supplemental .txts (for
% .tif data), consolidated inspecting single files with batch processing,
% added lifetime amplitude comparison
%%% v6 - added trimfirstT for multi-FOV BLIM (piezostage stabilization)
%%% v7 - changed trimfirst/lastT to allow user-defined number of points,
% added photobleach correction (not good)
%%% v8 - added option to adjust time delay for probe tuning
%%% v9 - improved time delay compensation with interpolation (and 
% preserving full data), added normalization by initial reference scan and
% PMT gain, converted reloading previous data to a config option

% Initialize
clear; clc; close all;

% Configuration options - check before running!
reloadprevious = 0; % 0 = run new analysis, 1 = reload previous analysis
testrun = 0; % 0 = process all, 1 = no-save test run with a few files,
% 2 = examine a single file
testfilesperfolder = 1; % number of files to test per folder (default 1)
targetfolders = [2 3 4]; % indices of folders to process, [] = dialog
filetype = 2; % 1 = .txt, 2 = .tif (solution), 3 = .tif (image)
filenameconv = 3; % 0 = other
% 1 = [pumpWL]-[pumppower]-[signalWL]-[IRpower]-[ND]-[PMTgain]-
% [chopper]-[PMTBW]-[etc],
% 2 = [FOV]_[etc]_[size]_[idler]_[DFG]_[power]_[channel]
% 3 = [IRWN]_[probeWL]_[etc]_[size]_[idler]_[DFG]_[power]_[channel]
t0pos = 225.3; % specify t0 position (mm), [] = autofind
tuningrate = 0.03; % PicoEmerald signal t0 tuning rate (ps/nm)
minprobe = 750; maxprobe = 880; % min and max probeWL for aligning tD
pairedDCAC = 1; % 1 = both DC and AC files present and paired, 0 = only AC
writeprocyn = 1; % 1 = write batch processed files, 0 = not
writefigsyn = 1; % 1 = write figure files, 0 = not
DFGyn = 1; % 1 = IR from DFG, 0 = IR from idler
powernormyn = 4; % 0 = no normalization, 1 = normalize by IR power,
% 2 = normalize by probe and IR powers, 3 = 2 + photobleach correction 
% 4 = 2 + reference scan normalization, 5 = 4 + PMT gain correction
tempmod = 1; % 1 = temporal modulation in place, 0 = no beamsplitters
ND = 1; % ND filter strength in probe path (script will update if possible)
PMT = 5; % PMT gain (script will update if possible)
prbpowerset = 250; % probe power in mW (script will update if possible)
IRpowerset = 70; % IR power in mW (script will update if possible)
normIRpower = 50; % IR power on-sample to normalize to (default is 50)
normprobepower = 1.5; % probe power on-sample to normalize to (default 1.5)
trimlastT = 1; % how many points to remove from end of Tlist (default 0)
trimfirstT = 0; % points to remove from start of Tlist (default 0)
basefittype = 1; % 0 = no baseline fit, 1 = linear baseline fit,
% 2 = exponential baseline fit
cutlow = 0.03; % lower baseline fit cutoff, (default 10th %ile)
cuthigh = 0.90; % upper baseline fit cutoff, (default 90th %ile)
ltfittype = 2; % 0 = no lifetime fitting, 1 = Gaussian*monoexp,
% 2 = Gaussian*biexp, 3 = beating Gaussian*exp, 4 = beating Gaussian*biexp
setpulsewidth = 2.3; % define pulse width (ps) in fit, [] = float
ltmin = [];  % define minimum lifetime (ps) in fit, [] = default (0.1)
ltmax = []; % define maximum lifetime (ps) in fit, [] = default (100)
annotatefigure = 1; % 1 = annotate figure (unfinished, leave as 1)

D = pwd; % get current directory
S = dir(fullfile(D,'*')); % search current directory
N = setdiff({S([S.isdir]).name},{'.','..'}); % subfolders of D

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
indrIR = 1; % "indr" = row index
indrprobe = indrIR+1; % 2
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
indrDCpeak = indralignlength+1; % 19
indrDCpeaksd = indrDCpeak+1; % 20
indrintegDC = indrDCpeaksd+1; % 21
indrDCintegsd = indrintegDC+1; % 22
indrDCsnr = indrDCintegsd+1; % 23
indrDCsbr = indrDCsnr+1; % 24
indrDCnoise = indrDCsbr+1; % 25
indrDCbleach = indrDCnoise+1; % 26 - minimum value of maxTlength

if reloadprevious == 0 % run new analysis
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
                if ltfittype > 0 % check if raw data work for convolution fitting
                    deltax = zeros(length(x)-1,1);
                    for i=1:length(x)-1
                        deltax(i) = round(1e2*(x(i+1)-x(i)))/1e2;
                    end
                    xspacing = unique(deltax);
                    if length(xspacing) == 1 && mod(length(x),2) == 1 % data ok
                    else % add extra length for interpolation
                        currentTlength = 2*currentTlength-1;
                    end
                end
                if currentTlength > maxTlength
                    maxTlength = currentTlength;
                end
            end
        end
    end
    
    if 10*maxTlength-1 < indrDCbleach % make sure master array is large enough
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
    master = NaN(maxTlength,indcDCbase-3*(1-pairedDCAC),length(folders),maxfilesinfolder); % add 3 columns for DC raw, corrected, & baseline fit
    
    if powernormyn > 0 % prep for power normalization
        IRpower = []; prbpower = [];
        if tempmod == 1 % pre-load beamsplitter response curve if needed
            BSdata = importdata('/Users/pkocheril/Documents/Caltech/Wei Lab/Spreadsheets/BP145B2/BP145B2.csv');
        end
        if powernormyn > 2 % normalization across samples
            normarray = NaN(maxTlength,indcDCbase-3,2,numel(N));
            normstop = 5-powernormyn; % 2 if powernormyn = 3, 1 if powernormyn = 4 -> 5-powernormyn
            for ii=1:normstop % loop over normalization folder(s)
                if filetype == 1
                    normfilelist = dir(fullfile(D,N{end-(2-ii)},'*.txt'));
                    norminfofilesraw = [];
                else
                    normfilelist = dir(fullfile(D,N{end-(2-ii)},'*CH2*.tif'));
                    norminfofilesraw = dir(fullfile(D,N{end-(2-ii)},'*IRsweep*.txt'));
                end
                normfiles = {normfilelist(~[normfilelist.isdir]).name}; % norm files
                norminfofiles = {norminfofilesraw(~[norminfofilesraw.isdir]).name}; % info files
                if numel(normfiles) == numel(norminfofiles)
                    norminfofound = 1;
                else
                    norminfofound = 0;
                end
                for jj=1:numel(normfiles)
                    F = fullfile(D,N{end-(2-ii)},normfiles{jj});
                    
                    % Filename parsing
                    folderinfo = split(F,"/"); % split path into folders
                    if filenameconv == 1 % parse .txts
                        fileinfo = split(folderinfo(end),"-"); % split filename at "-"s
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
                        if DFGyn == 1
                            IRWN = 1e7/DFGWL;
                        else
                            IRWN = 1e7/idlerWL;
                        end
                        NDstrfull = char(fileinfo(5)); % extract ND
                        NDstr = NDstrfull(end); % remove "ND"
                        ND = sscanf(NDstr,'%f');
                    else % parse default .tifs
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
                        prbWL = 1; % <-- unknown from filename
                        idlerWL = 1E7/idlerWN;
                        signalWL = 1/((1/1031.2)-(1/idlerWL));
                        IRpower = IRpowermeter*300; % assuming 300 mW scale
                        if DFGyn == 1
                            IRWN = DFGWN;
                        else
                            IRWN = idlerWN;
                        end
                    end
                    if norminfofound == 1 % parse info from IRsweep.txt files
                        infofilename = fullfile(D,N{end-(2-ii)},norminfofiles{jj});
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
                        prbWL = round(infofile(13));
                        prbpower = infofile(14);
                        ND = infofile(15);
                        idlerWN = round(infofile(16));
                        DFGWN = round(infofile(17));
                        IRpower = infofile(18);
                        if DFGyn == 1
                            IRWN = DFGWN;
                        else
                            IRWN = idlerWN;
                        end
                        dwelltime = infofile(19);
                    end
            
                    % Load data - AC only for normalization
                    if filetype == 1 % load data from .txt
                        data = importdata(F); % load data
                        x = data(:,1); % save delay position as x
                        if width(data) >= 3
                            signal = data(:,3); % save channel 2 as signal
                        else
                            signal = data(:,2);
                        end
                        x = x(trimfirstT+1:end-trimlastT);
                        signal = signal(trimfirstT+1:end-trimlastT);
                        % No standard deviation data present - write as zeros
                        sigsds = zeros(height(signal),width(signal));
                    else % load Tlist and .tif
                        Tlistname = fullfile(D,N{end-(2-ii)},'Tlist.txt'); % current Tlist
                        x = importdata(Tlistname); % load Tlist
                        data = double(tiffreadVolume(F));
                        imagesize = size(data);
                        x = x(trimfirstT+1:end-trimlastT);
                        data = data(:,:,trimfirstT+1:end-trimlastT);
                        if length(x) > length(data)
                            x = x(1:length(data));
                        end
                        if filetype == 2 % average & stdev of XY data in .tif
                            signal = squeeze(mean(data, [1 2]));
                            sigsds = squeeze(std(data, 0, [1 2]));
                            imagesize = size(signal);
                        else % image data - make arrays to write outputs into
                            r2array = zeros([imagesize(1) imagesize(2)]);
                            lt1array = r2array; lt2array = r2array;
                            ssresidarray = r2array; FWHMarray = r2array;
                        end
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
                    
                            % Power normalize
                            if isempty(IRpower) == 1 % if no IR power found
                                IRpower = IRpowerset; % set to power from config
                            end
                            if isempty(prbpower) == 1 % if no probe power found
                                prbpower = prbpowerset; % set to power from config
                            end
                            signal = signal*normIRpower/IRpower;
                            sigsds = sigsds*normIRpower/IRpower;
                            if powernormyn > 1 % IR and probe
                                if tempmod == 1 % correct for temporal modulation
                                    prbpower = 0.25*prbpower; % loss from pinhole, 250 mW -> 62.5 mW
                                    transmittance = BSdata.data(prbWL,1)/100; % picoEmerald is p-polarized (horizontal)
                                    reflectance = BSdata.data(prbWL,4)/100; % CSV in %
                                    prbpower = transmittance*reflectance*prbpower; % 62.5 mW -> 15.6 mW
                                end
                                prbpower = prbpower/(10^(ND)); % ND1: 15.6 mW -> 1.56 mW
                                signal = signal*(normprobepower/prbpower);
                                sigsds = sigsds*(normprobepower/prbpower);
                            end
                            
                            % Setup baseline fitting
                            ilow = ceil(length(t)*cutlow); % using ceiling to round up
                            ihigh = ceil(length(t)*cuthigh);
                            tbase = [t(2:ilow); t(ihigh:end)]; % trimmed t
                            sigbase = [signal(2:ilow); signal(ihigh:end)]; % trimmed signal
                            sbr = max(signal(ilow:ihigh))/mean(signal(ihigh:end));
                            
                            % Baseline fitting
                            basefit = fit(tbase,sigbase,'poly1'); % fit baseline
                            basecoef = coeffvalues(basefit); % 1 = slope, 2 = int
                            basecurve = polyval(basecoef,t);
                            
                            % Calculate peak height and integrated signal
                            corrsig = signal - basecurve;
                            [sigpeak,maxindex] = max(corrsig);
                            peaksd = sigsds(maxindex);
                            integsig = sum(corrsig);
                            integsd = mean(sigsds,'all');
                            corrsigbase = [corrsig(2:ilow); corrsig(ihigh:end)];
                            noise = range(corrsigbase);
                            basesd = std(corrsigbase);
                            snr = sigpeak/noise;
                            snrsweep = corrsig/noise;
                            intsnr = sum(snrsweep);
                            
                            % Make figure annotation
                            annot = {'IR = '+string(IRWN)+' cm-1','probe = '+...
                                string(prbWL)+' nm'};
                            annotdim = [0.35 0.6 0.2 0.2];
                            annot(length(annot)+1) = {'Peak height = '+...
                                string(sigpeak)};
                            outname = string(F(1:end-4));
                            fig = figure('visible','off');
                            hold on; errorbar(t,signal,sigsds,'o')
                            plot(tbase,sigbase,'go',t,basecurve,'-'); hold off;
                            xlabel('Time (ps)'); ylabel('AC signal (AU)');
                            legend('Data','Baseline data','Baseline');
                            annotation('textbox',annotdim,'String',annot);
                            if writefigsyn == 1
                                saveas(fig,outname+'_proc.png')
                            end
    
                            normarray(1:length(t),indct,ii,jj) = t(:);
                            normarray(1:length(corrsig),indccorrsig,ii,jj) = corrsig(:);
                            normarray(indrIR,indcvalue,ii,jj) = IRWN;
                            normarray(indrprobe,indcvalue,ii,jj) = prbWL;
                            normarray(indrsigpeak,indcvalue,ii,jj) = sigpeak;
                            normarray(indrpeaksd,indcvalue,ii,jj) = peaksd;
                            normarray(indrintegsig,indcvalue,ii,jj) = integsig;
                            normarray(indrintegsd,indcvalue,ii,jj) = integsd;
                            normarray(indrsnr,indcvalue,ii,jj) = snr;
                            normarray(indrsbr,indcvalue,ii,jj) = sbr;
                            normarray(indrnoise,indcvalue,ii,jj) = noise;
                            normarray(indrbasesd,indcvalue,ii,jj) = basesd;
                            normarray(indrintsnr,indcvalue,ii,jj) = intsnr;
                            normarray(indrtlist,indcvalue,ii,jj) = length(t);
                            normarray(1:length(signal),indcrawsig,ii,jj) = signal(:);
                            
                            if writeprocyn == 1 % write files
                                temp = normarray(:,:,ii,jj);
                                writematrix(temp,outname+'_norm.dat',...
                                    'FileType','text')
                            end
                        end
                    end
                    prbpower = []; % clear powers for next file
                    IRpower = []; % (otherwise, normalization can compound)
                end
            end
        end
    end
    
    % Loop through all files and analyze
    for ii = folders % ii = subfolder number
        if filetype == 1  % look for .txt files in subfolders
            T = dir(fullfile(D,N{ii},'*.txt'));
            infofilesraw = [];
        else % look for CH2 tifs
            T = dir(fullfile(D,N{ii},'*CH2*.tif'));
            infofilesraw = dir(fullfile(D,N{ii},'*IRsweep*.txt'));
        end
        C = {T(~[T.isdir]).name}; % all data files
        infofiles = {infofilesraw(~[infofilesraw.isdir]).name}; % info files
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
                    if DFGyn == 1
                        IRWN = 1e7/DFGWL;
                    else
                        IRWN = 1e7/idlerWL;
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
                    imgdim = char(fileinfo(end-4)); % '5X5'
                    prbWL = 1; % <-- unknown from filename
                    idlerWL = 1E7/idlerWN;
                    signalWL = 1/((1/1031.2)-(1/idlerWL));
                    IRpower = IRpowermeter*300; % assuming 300 mW scale
                    if DFGyn == 1
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
                    prbWL = round(infofile(13));
                    prbpower = infofile(14);
                    ND = infofile(15);
                    idlerWN = round(infofile(16));
                    DFGWN = round(infofile(17));
                    IRpower = infofile(18);
                    if DFGyn == 1
                        IRWN = DFGWN;
                    else
                        IRWN = idlerWN;
                    end
                    dwelltime = infofile(19);
                end
                if filenameconv == 3 % parse 2D sweep .tifs
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
                            if totalfilesC <= 70
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
                    imgdim = char(fileinfo(end-4)); % '5X5'
                    idlerWL = 1E7/idlerWN;
                    signalWL = 1/((1/1031.2)-(1/idlerWL));
                    IRpower = IRpowermeter*300; % assuming 300 mW scale
                end
    
                % Load data
                if filetype == 1 % load data from .txt
                    data = importdata(F); % load data
                    x = data(:,1); % save delay position as x
                    if pairedDCAC == 1
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
                    x = x(trimfirstT+1:end-trimlastT);
                    DC = DC(trimfirstT+1:end-trimlastT);
                    signal = signal(trimfirstT+1:end-trimlastT);
                    % No standard deviation data present - write as zeros
                    sigsds = zeros(height(signal),width(signal));
                    if pairedDCAC == 1
                        DCsds = zeros(height(signal),width(signal));
                    end
                else % load Tlist and .tif
                    Tlistname = fullfile(D,N{ii},'Tlist.txt'); % current Tlist
                    x = importdata(Tlistname); % load Tlist
                    data = double(tiffreadVolume(F));
                    imagesize = size(data);
                    if pairedDCAC == 1
                        F_DC = F;
                        F_DC(end-4) = "1";
                        DCdata = double(tiffreadVolume(F_DC));
                    else
                        DCdata = zeros(height(data),width(data),length(x));
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
                            if pairedDCAC == 1
                                DC = squeeze(DCdata(iii,jjj,:));
                                DCsds = zeros(sigsds);
                            else
                                DC = sigsds; DCsds = sigsds;
                            end
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

                        if powernormyn > 0 % prepare for power normalization
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

                            if powernormyn > 1 % IR and probe
                                if tempmod == 1 % correct for temporal modulation
                                    prbpower = 0.25*prbpower; % loss from pinhole, 250 mW -> 62.5 mW
                                    transmittance = BSdata.data(prbWL,1)/100; % picoEmerald is p-polarized (horizontal)
                                    reflectance = BSdata.data(prbWL,4)/100; % CSV in %
                                    prbpower = transmittance*reflectance*prbpower; % 62.5 mW -> 15.6 mW
                                end
                                prbpower = prbpower/(10^(ND)); % ND1: 15.6 mW -> 1.56 mW
                                signal = signal*(normprobepower/prbpower);
                                sigsds = sigsds*(normprobepower/prbpower);
                                if pairedDCAC == 1
                                    DC = DC*(normprobepower/prbpower);
                                    DCsds = DCsds*(normprobepower/prbpower);
                                end
                                if powernormyn == 3 % photobleach correction
                                    progscaling = (jj-startjj)/(totalfilesC-startjj);
                                    normpeak = normarray(indrsigpeak,indcvalue,1,ii);
                                    pbpeak = normarray(indrsigpeak,indcvalue,2,ii);
                                    photobleachloss = normpeak-pbpeak;
                                    pbcorrfactor = normpeak/(normpeak-(progscaling*photobleachloss));
                                    if pbcorrfactor <= 10
                                        signal = signal*pbcorrfactor;
                                        sigsds = sigsds*pbcorrfactor;
                                        if pairedDCAC == 1
                                            DC = DC*pbcorrfactor;
                                            DCsds = DCsds*pbcorrfactor;
                                        end
                                    else
                                        if pbcorrfactor <= 100
                                            warning('Large photobleaching - using log-scale correction')
                                            signal = signal*(1+log(pbcorrfactor));
                                            sigsds = sigsds*(1+log(pbcorrfactor));
                                            if pairedDCAC == 1
                                                DC = DC*(1+log(pbcorrfactor));
                                                DCsds = DCsds*(1+log(pbcorrfactor));
                                            end
                                        else
                                            warning('Massive photobleaching - correction failed')
                                        end
                                    end
                                end
                                if powernormyn == 4 % reference scan normalization
                                    normpeak = normarray(indrsigpeak,indcvalue,1,ii);
                                    normfactor = 0.5/normpeak; % normalizing to peak height of 0.5;
                                    signal = signal*normfactor;
                                    sigsds = sigsds*normfactor;
                                    if pairedDCAC == 1
                                        DC = DC*normfactor;
                                        DCsds = DCsds*normfactor;
                                    end
                                end
                                if powernormyn == 5 % reference scan normalization and PMT gain correction
                                    normpeak = normarray(indrsigpeak,indcvalue,1,ii);
                                    % PMT gain stuff
                                    gaincurve = importdata('/Users/pkocheril/Documents/Caltech/Wei Lab/Data/2023_04_28/PMT Gain Data/gain_calibration_curve.txt');
                                    gainfactor = gaincurve(PMT,2);
                                    normfactor = 0.5/(normpeak*gainfactor); % normalizing to peak height of 0.5, PMT gain 10
                                    signal = signal*normfactor;
                                    sigsds = sigsds*normfactor;
                                    if pairedDCAC == 1
                                        DC = DC*normfactor;
                                        DCsds = DCsds*normfactor;
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
                        if pairedDCAC == 1
                            DCbase = [DC(2:ilow); DC(ihigh:end)]; % trimmed DC
                            DCsbr = max(DC(ilow:ihigh))/mean(DC(ihigh:end));
                        end
                        
                        % Fitting
                        if basefittype == 0 % no baseline fitting, set basecurve to 0
                            basecurve = zeros(height(t),width(t));
                            bleachrate = 0;
                            if pairedDCAC == 1
                                DCbasecurve = zeros(height(t),width(t));
                                DCbleachrate = 0;
                                DCsbr = 0;
                            end
                        else
                            if basefittype == 2 % fit baseline to exponential w/ vertical offset
                                basefn = @(b) b(1)+b(2)*exp(-b(3)*(tbase-b(4)))-sigbase; % error function
                                baseopt=optimoptions(@lsqnonlin);
                                baseopt.MaxFunctionEvaluations = 1e6; % let the fit run longer
                                baseopt.MaxIterations = 1e6; % let the fit run longer
                                baseopt.FunctionTolerance = 1e-12; % make the fit more accurate
                                baseopt.OptimalityTolerance = 1e-12; % make the fit more accurate
                                baseopt.StepTolerance = 1e-12; % make the fit more precise
                                baseopt.Display = 'off'; % silence console output
                                basegs = [min(sigbase) 0.1*(max(sigbase)-min(sigbase)) 0 0];
                                basefit = lsqnonlin(basefn,basegs,[],[],baseopt);
                                basecurve = basefit(1)+basefit(2)*exp(-basefit(3)*(t-basefit(4)));
                                if basefit(2) > 0 && basefit(3) > 0 
                                    bleachrate = basefit(3); % ps-1
                                else % bleach rate doesn't have physical meaning
                                    bleachrate = 0;
                                end
                                if pairedDCAC == 1 % fit DC baseline
                                    DCbasefn = @(b) b(1)+b(2)*exp(-b(3)*(tbase-b(4)))-DCbase;
                                    DCbgs = [min(DCbase) 0.1*(max(DCbase)-min(DCbase)) 0 0];
                                    DCbasefit = lsqnonlin(DCbasefn,DCbgs,[],[],baseopt);
                                    DCbasecurve = DCbasefit(1)+DCbasefit(2)*exp(-DCbasefit(3)*(t-DCbasefit(4)));
                                    if DCbasefit(2) > 0 && DCbasefit(3) > 0
                                        DCbleachrate = DCbasefit(3); % ps-1
                                    else
                                        DCbleachrate = 0;
                                    end
                                end
                            else
                                if basefittype == 1 % linear baseline fit
                                    basefit = fit(tbase,sigbase,'poly1'); % fit baseline
                                    basecoef = coeffvalues(basefit); % 1 = slope, 2 = int
                                    basecurve = polyval(basecoef,t);
                                    if basecoef(1) <= 0 
                                        bleachrate = -basecoef(1); % ps^-1
                                    else
                                        bleachrate = 0;
                                    end
                                    if pairedDCAC == 1
                                        DCbasefit = fit(tbase,DCbase,'poly1');
                                        DCbasecoef = coeffvalues(DCbasefit);
                                        DCbasecurve = polyval(DCbasecoef,t);
                                        if basecoef(1) <= 0 
                                            DCbleachrate = -DCbasecoef(1); % ps^-1
                                        else
                                            DCbleachrate = 0;
                                        end
                                    else
                                        DCbleachrate = 0;
                                    end
                                end
                            end
                        end

                        % Calculate t spacing
                        deltat = zeros(length(t)-1,1); % check t spacing
                        for i=1:length(t)-1
                            deltat(i) = round(1e2*(t(i+1)-t(i)))/1e2; % round to nearest 0.01 ps (had issues without rounding)
                        end
                        tspacing = unique(deltat); % length 1 <-> evenly spaced t

                        if ltfittype > 0 % lifetime fitting
                            % Check t spacing (even) and length (odd)
                            if length(tspacing) == 1  && mod(length(t),2) == 1
                                tint = t; sigint = signal;
                            else % interpolate for convolution
                                tint = linspace(t(1),t(end),length(t)*2-1).';
                                sigint = spline(t, signal, tint);
                            end
                            tc = tint(1:2:end); % odd values of tint --> convolution
                            if powernormyn > 0
                                IRFamp = IRpower/1e3; % IRF amplitude (~ IR power in W)
                            else
                                IRFamp = 50/1e3; % default to 50 mW
                            end
                        
                            % Define error function(r) as fit minus signal
                            fun = @(r) conv(IRFamp*exp(-r(2)*(tc-r(1)).^2), ...
                                r(3)*exp(-r(4)*(tc))+r(5)*exp(-r(6)*(tc))+...
                                r(7)*exp(-r(8)*(tc)).*cos(r(9)*tc-r(10))+ ...
                                r(11)*exp(-r(12)*(tc)).*cos(r(13)*tc-r(14))+ ...
                                r(15)*exp(-r(16)*(tc)).*cos(r(17)*tc-r(18)))+ ...
                                r(19)*tint+r(20)+r(21)*exp(-r(22)*(tint-r(23)))-sigint;
                            
                            % Initial guesses for r - Gaussian terms
                            cent = 10; % r(1), ps -- center
                            gdec = 4*log(2)/(4^2); % r(2), ps^-2 -- Gaussian decay
                            
                            % Initial guesses for r - exponential terms
                            eamp1 = 0.05; % r(3) -- amplitude 1
                            edec1 = 1/2; % r(4), ps^-1 -- 1/lifetime 1
                            eamp2 = 0.005; % r(5)  -- amplitude 2
                            edec2 = 1/10; % r(6), ps^-1  -- 1/lifetime 2
                            
                            % Initial guesses for r - beating terms
                            bamp1 = 0.1; % r(7) -- beat amplitude 1
                            bdec1 = 0.1; % r(8), ps^-1  -- beat decay 1
                            fr1 = 0.5; % r(9), rad ps^-1  -- beat freq 1
                            phs1 = 0; % r(10), rad -- beat phase 1
                            bamp2 = 0.1; % r(11) -- beat amplitude 2
                            bdec2 = 0.1; % r(12), ps^-1 -- beat decay 2
                            fr2 = 0.5; % r(13), rad ps^-1 -- beat freq 2
                            phs2 = 0; % r(14), rad -- beat phase 2
                            bamp3 = 0.1; % r(15) -- beat amplitude 3
                            bdec3 = 0.1; % r(16), ps^-1 -- beat decay 3
                            fr3 = 0.5; % r(17), rad ps^-1 -- beat freq 3
                            phs3 = 0; % r(18), rad -- beat phase 3
                                            
                            % Fit setup
                            r0 = [cent,gdec,eamp1,edec1,eamp2,edec2,... % Gaussian and exp decays
                                bamp1,bdec1,fr1,phs1,... % first beating term
                                bamp2,bdec2,fr2,phs2... % second beating term
                                bamp3,bdec3,fr3,phs3,... % third beating term
                                0,0,0,0,0]; % baseline terms
                            lb = [-10,0.01,0,1/100,0,1/100,... % lower bounds
                                0,0,0,-180,... % default lifetime max = 100 ps
                                0,0,0,-180,...
                                0,0,0,-180,...
                                0,0,0,0,0];
                            ub = [20,9,9,1/0.1,9,1/0.1,... % upper bounds
                                9,9,9,180,... % default lifetime min = 0.1 ps
                                9,9,9,180,...
                                9,9,9,180,...
                                9,9,9,9,9];
                            
                            % Implement fit options specified in config
                            if basefittype == 1 % linear baseline
                                lb(21:23) = 0; ub(21:23) = 0; % mute exponential terms
                                lb(19) = basecoef(1); ub(19) = lb(19);
                                lb(20) = basecoef(2); ub(20) = lb(20);
                            end
                            if basefittype == 2 % exponential baseline
                                lb(19) = 0; ub(19) = 0; % mute linear term
                                lb(20) = basefit(1); ub(20) = lb(20);
                                lb(21) = basefit(2); ub(21) = lb(21); 
                                lb(22) = basefit(3); ub(22) = lb(22);
                                lb(23) = basefit(4); ub(23) = lb(23);
                            end
                            if ltfittype == 1 % Gaussian*exp
                                lb(5:18) = 0; ub(5:18) = 0;
                            end
                            if ltfittype == 2 % Gaussian*biexp
                                lb(7:18) = 0; ub(7:18) = 0;
                            end
                            if ltfittype == 3 % Gaussian*exp with beating
                                lb(5:6) = 0; ub(5:6) = 0;
                            end
                            if isempty(setpulsewidth) == 0 % pulse width was specified
                                pulsedecay = log(2)/setpulsewidth^2;
                                lb(2) = pulsedecay; ub(2) = pulsedecay;
                            end
                            if isempty(ltmin) == 0 % minimum lifetime was specified
                                ub(4) = 1/ltmin; ub(6) = ub(4); % min lt -> ub decay
                            end
                            if isempty(ltmax) == 0 % maximum lifetime was specified
                                lb(4) = 1/ltmax; lb(6) = lb(4); % max lt -> lb decay
                            end
                        
                            % Run lsqnonlin - minimizes error function by adjusting r
                            options=optimoptions(@lsqnonlin);
                            options.MaxFunctionEvaluations = 1e6; % let the fit run longer
                            options.MaxIterations = 1e6; % let the fit run longer
                            options.FunctionTolerance = 1e-12; % make the fit more accurate
                            options.Display = 'off'; % silence console output
                            fitval = lsqnonlin(fun,r0,lb,ub,options); % fit coeffs
                        
                            % Generate curve for fitted parameters
                            fitcurve = conv(IRFamp*exp(-fitval(2)*(tc-fitval(1)).^2),...
                                fitval(3)*exp(-fitval(4)*(tc))+ ...
                                fitval(5)*exp(-fitval(6)*(tc))+ ...
                                fitval(7)*exp(-fitval(8)*(tc)).* ...
                                cos(fitval(9)*tc-fitval(10))+ ...
                                fitval(11)*exp(-fitval(12)*(tc)).* ...
                                cos(fitval(13)*tc-fitval(14))+ ...
                                fitval(15)*exp(-fitval(16)*(tc)).* ...
                                cos(fitval(17)*tc-fitval(18)))+...
                                fitval(19)*tint+fitval(20)+fitval(21)*...
                                exp(-fitval(22)*(tint-fitval(23)));
                            
                            % Calculate residuals, ssresid, and IRF FWHM
                            resid = sigint-fitcurve;
                            ssresid = sum(resid.^2); % sum of squares of residuals
                            tss = sum((sigint-mean(sigint)).^2); % total sum of squares
                            r2 = 1-(ssresid/tss); % coefficient of determination (R^2)
                            if r2 < 0 % remove nonsense R^2 values
                                r2 = 0;
                            end
                            FWHM = sqrt(log(2)/fitval(2)); % pulse width, ps
                            
                            % Lifetime assignment - want lifetime1 < lifetime2
                            if abs(1/fitval(4)) < abs(1/fitval(6))
                                if fitval(4) ~= 0 % could be 0 from fit type options
                                    lifetime1 = 1/fitval(4);
                                else
                                    lifetime1 = 0;
                                end
                                if fitval(6) ~= 0
                                    lifetime2 = 1/fitval(6);
                                else
                                    lifetime2 = 0;
                                end
                                lt1lt2 = fitval(3)/fitval(5); % A1/A2
                            else % means 1/fitval(4) is longer
                                if fitval(4) ~= 0
                                    lifetime2 = 1/fitval(4);
                                else
                                    lifetime2 = 0;
                                end
                                if fitval(6) ~= 0
                                    lifetime1 = 1/fitval(6);
                                else
                                    lifetime1 = 0;
                                end
                                lt1lt2 = fitval(5)/fitval(3); % still A1/A2
                            end
                        else % no lifetime fitting
                            tint = t; sigint = signal; fitcurve = signal;
                            lifetime1 = 0; lifetime2 = 0; lt1lt2 = 0;
                            fitval = zeros(23,1); r2 = 0;
                        end
                        
                        if filetype == 3
                            ssresidarray(iii,jjj) = ssresid;
                            r2array(iii,jjj) = r2;
                            FWHMarray(iii,jjj) = FWHM;
                            lt1array(iii,jjj) = lifetime1;
                            lt2array(iii,jjj) = lifetime2;
                        end
                        
                        % Make rounded strings for figure annotation
                        sDCbleachrate = string(sprintf('%0.2g',DCbleachrate)); % round to 2 sig figs
                        lt1 = string(round(10*lifetime1)/10); % round to nearest 0.1 ps
                        lt2 = string(round(10*lifetime2)/10);
                        beatstring = string(fitval(9)*1e12/(pi*29979245800))+...
                            ', '+string(fitval(13)*1e12/(pi*29979245800))+...
                            ', '+string(fitval(17)*1e12/(pi*29979245800));
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
                            if pairedDCAC == 1
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
                                    ' ps, beats = '+beatstring+' cm^{-1}, r^2 = '+sr2};
                            end
                            if ltfittype == 4
                                annot(length(annot)+1) = {'τ_1 = '+lt1+' ps, τ_2 = '+...
                                    lt2+' ps, A_{1}/A_{2} = '+slt1lt2+', beats = '+...
                                    beatstring+' cm^{-1}, r^2 = '+sr2};
                            end
                            if ltfittype > 0 && r2 < 0.7
                                annot(end) = {'Poor fit, r^2 = ' + sr2};
                            end
                        end
                        
                        % Prepare output filename(s)
                        if filetype > 1 && pairedDCAC == 1
                            outname = string(F(1:end-8)) + "_DCAC";
                        else
                            outname = string(F(1:end-4));
                        end
                        
                        % Calculate peak height and integrated signal
                        corrsig = signal - basecurve;
                        [sigpeak,maxindex] = max(corrsig);
                        peaksd = sigsds(maxindex);
                        integsig = sum(corrsig);
                        integsd = mean(sigsds,'all');
                        corrsigbase = [corrsig(2:ilow); corrsig(ihigh:end)];
                        noise = range(corrsigbase);
                        basesd = std(corrsigbase);
                        snr = sigpeak/noise;
                        snrsweep = corrsig/noise;
                        intsnr = sum(snrsweep);

                        % Interpolate for tD compensation
                        alignlength = ceil(((t(end)-t(1))/(tuningrate/2))+1); % want talignspacing >= tuningrate/2
                        talign = linspace(t(1),t(end),alignlength).';
                        sigalign = spline(t, corrsig, talign);
                        talignspacing = talign(2)-talign(1); % evenly spaced

                        % tD compensation
                        fulltimeoffset = tuningrate*(maxprobe-minprobe); % ps
                        currentprobetimeoffset = tuningrate*(prbWL-minprobe); % ps
                        alignstart = floor(currentprobetimeoffset/talignspacing); % index
                        maxalignstart = floor(fulltimeoffset/talignspacing);
                        alignend = length(talign)-(maxalignstart-alignstart);
                        alignstart = alignstart+1;
                        talign = talign(alignstart:alignend);
                        sigalign = sigalign(alignstart:alignend);
                
                        if ltfittype > 0 % prepare figure legend
                            legend1 = {'Data','Baseline data','Baseline','Fit'};
                        else
                            legend1 = {'Data','Baseline data','Baseline'};
                        end
                        
                        if pairedDCAC == 1 % make figure
                            corrDC = DC - DCbasecurve;
                            [DCpeak,DCmaxindex] = max(corrDC);
                            DCpeaksd = DCsds(DCmaxindex);
                            integDC = sum(corrDC);
                            DCintegsd = mean(DCsds,'all');
                            corrDCbase = [corrDC(2:ilow); corrDC(ihigh:end)];
                            DCnoise = range(corrDCbase);
                            DCsnr = DCpeak/DCnoise;
                        
                            % Plot DC and AC side-by-side
                            fig = figure('visible',visibility); 
                            if ltfittype > 0
                                tiledlayout(2,2);
                            else
                                tiledlayout(1,2);
                            end
                            nexttile([1 1]); 
                            hold on; errorbar(t,DC,DCsds,'o');
                            plot(tbase,DCbase,'o',t,DCbasecurve,'-'); hold off;
                            xlabel('Time (ps)'); ylabel('DC signal (AU)');
                            legend('Data','Baseline data','Baseline');
                            nexttile([1 1]); 
                            if ltfittype > 0
                                plot(tint,resid); xlabel('Time (ps)'); 
                                ylabel('Residuals (AU)'); legend('Residuals');
                                nexttile([1 2]);
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
                                tiledlayout(3,1); nexttile([1 1]);
                                plot(tint,resid); xlabel('Time (ps)'); 
                                ylabel('Residuals (AU)'); legend('Residuals');
                                nexttile([2 1]);
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
                            master(indrIR,indcvalue,ii,jj) = IRWN;
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
                
                            if pairedDCAC == 1
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
            end
        end
    end
    
    if writeprocyn == 1 && filetype ~= 3 % write batch file
        writematrix(master,'batch_flattenedv9.md','FileType','text') % backup, not very useful
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
    for ii = 1:numel(N)
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
figvis = 'on';

% Master array is 4D (x,y,ii,jj); ii,jj are folder,file
% x = rows (individual t points); y = columns
master_size = importdata('master_size.yml'); % load dimensions

% To pull columns: squeeze(master(1:TL,[indexcolumn],[folder],1:NF))
% TL = master(indrtlist,indcalue,[folder],1)
% NF = master(indrnfiles,indcvalue,[folder],1)
% (1:NF optional - can also just use :)

% To pull specific values: squeeze(master(indr...,indcvalue,[folder],1:NF))

% Use 'rmmissing' to remove NaN values

% Contours of signal vs t and probe at each IR freq
time = squeeze(master(1:master(indrtlist,indcvalue,2,1),indct,2,1)).';
probe = squeeze(master(indrprobe,indcvalue,2,1:master(indrnfiles,indcvalue,2,1)));
csig = NaN(length(probe),length(time),3);
csig(:,:,1) = squeeze(master(1:master(indrtlist,indcvalue,2,1),indccorrsig,2,:)).'; % 1500
csig(:,:,2) = squeeze(master(1:master(indrtlist,indcvalue,3,1),indccorrsig,3,:)).'; % 1510
csig(:,:,3) = squeeze(master(1:master(indrtlist,indcvalue,4,1),indccorrsig,4,:)).'; % 1520

[t1,probe1] = meshgrid(time,probe);

logsig(:,:,:) = log(abs(csig(:,:,:)));

% Linear scale
for i=1:3
    figure;%('visible','off');
    IRfreq = 1490+(10*i);
    contourf(t1,probe1,csig(:,:,i)); cb=colorbar;
    xlabel('Time delay (ps)'); ylabel('λ_{probe} (nm)');
    cb.Label.String='Corrected signal at '+string(IRfreq)+' cm^{-1} (AU)';
    cb.Label.Rotation=270; cb.Label.VerticalAlignment = "bottom";
end

% % Log scale - not useful
% for i=1:3
%     figure('visible','off');
%     IRfreq = 1490+(10*i);
%     contourf(t1,probe1,logsig(:,:,i)); cb=colorbar;
%     xlabel('Adjusted time delay (ps)'); ylabel('λ_{probe} (nm)');
%     cb.Label.String='Corrected signal at '+string(IRfreq)+' cm^{-1} (AU)';
%     cb.Label.Rotation=270; cb.Label.VerticalAlignment = "bottom";
% end

% Time delay compensated contours
timealign = squeeze(master(1:master(indralignlength,indcvalue,2,1),indctalign,2,1)).';
sigalign = NaN(length(probe),length(timealign),3);
sigalign(:,:,1) = squeeze(master(1:master(indralignlength,indcvalue,2,1),indcsigalign,2,:)).'; % 1500
sigalign(:,:,2) = squeeze(master(1:master(indralignlength,indcvalue,3,1),indcsigalign,3,:)).'; % 1510
sigalign(:,:,3) = squeeze(master(1:master(indralignlength,indcvalue,4,1),indcsigalign,4,:)).'; % 1520

[t2,probe2] = meshgrid(timealign,probe);

for i=1:3
    figure;%('visible','off');
    IRfreq = 1490+(10*i);
    contourf(t2,probe2,sigalign(:,:,i)); cb=colorbar;
    xlabel('Compensated time delay (ps)'); ylabel('λ_{probe} (nm)');
    cb.Label.String='Corrected signal at '+string(IRfreq)+' cm^{-1} (AU)';
    cb.Label.Rotation=270; cb.Label.VerticalAlignment = "bottom";
end


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

