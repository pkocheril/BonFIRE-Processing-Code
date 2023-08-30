%% Fit temporal sweep data
%%% v2 - fit to Gaussian*biexponential with optional beating
%%% v3 - added ability to fit linear background and compared 
% trust-region-reflective vs levenberg-marquardt (TRR is better)
%%% v4 - added spline interpolation for nonlinearly spaced data and
% tried Gaussian IRF pre-fitting, but it doesn't work
%%% v5 - cleaned up v4 and fixed FWHM calculation, also verified that
% Gaussian amplitude term isn't needed and doesn't affect lifetimes
%%% v6 - added R^2 calculation
%%% v7 - updated based on sweep_processv3, improved handling of low-info 
% filenames
%%% v8 - added image fitting from image_fit_temporalv5, power loss from 
% pinhole
%%% v9 - updated based on sweep_proc_alignv3

% Load data
clear; close all; clc;

%F = '516nm_1500cm-1 (07)_5X5_Idler4063.7_DFG1570.0_P0.1294_CH2.tif';
F = '516nm_1500cm-1 (06)_5X5_Idler4068.7_DFG1560.1_P0.1127_CH2.tif';
Tlistname = 'Tlist.txt';


% Configuration options - check before running!
visibility = 'on';
writefigsyn = 0; % 1 = write figure files, 0 = not

filetype = 2; % 1 = .txt, 2 = .tif (solution), 3 = .tif (image)
filenameconv = 2; % 0 = other,
% 1 = [pumpWL]-[pumppower]-[signalWL]-[IRpower]-[ND]-[PMTgain]-[modfreq]-[PMTBW]-[etc],
% 2 = [etc]_[size]_[idler]_[DFG]_[power]_[channel],
% 3 = [IRWN]_[probeWL]_[etc]_[size]_[idler]_[DFG]_[power]_[channel],
% 4 = [probeWL]_[IRWN]_[etc]_[size]_[idler]_[DFG]_[power]_[channel]
t0pos = 194.9; % specify t0 position (mm), [] = autofind
pairedDCAC = 1; % 2 = SPCM + PMT CH2, 1 = PMT CH1 + CH2, 0 = only PMT CH2
DFGyn = []; % [] = auto-choose, 1 = IR from DFG, 0 = IR from idler
powernormyn = 0; % 0 = no normalization, 1 = normalize by IR power,
% 2 = normalize by probe and IR powers, 3 = 2 + photobleach correction(bad)
% 4 = 2 + reference scan normalization, 5 = 4 + PMT gain correction
tempmod = 0; % 1 = temporal modulation in place, 0 = no beamsplitters
ND = 2; % ND filter strength in probe path (script will update if possible)
PMT = 1; % PMT gain (script will update if possible)
prbpowerset = 250; % probe power in mW (script will update if possible)
IRpowerset = 70; % IR power in mW (script will update if possible)
normIRpower = 50; % IR power on-sample to normalize to (default is 50)
normprobepower = 1.5; % probe power on-sample to normalize to (default 1.5)
trimlastT = 4; % how many points to remove from end of Tlist (default 0)
trimfirstT = 0; % points to remove from start of Tlist (default 0)
basefittype = 2; % 0 = no baseline fit, 1 = linear fit, 2 = exp fit
cutlow = 0.05; % lower baseline fit cutoff, (default 10th %ile)
cuthigh = 0.95; % upper baseline fit cutoff, (default 90th %ile)
ltfittype = 1; % [] = auto-choose, 0 = no fitting, 1 = Gaussian*monoexp,
% 2 = Gaussian*biexp, 3 = beating Gaussian*exp, 4 = beating Gaussian*biexp,
% 5 = two Gauss*exp, 6 = Gauss*exp + Gauss*biexp, 7 = two Gauss*biexp,
% 8 = two beating Gauss*exp, 9 = beating Gauss*exp + beating Gauss*biexp,
% 10 = two beating Gauss*biexp, 11 = three Gauss*exp, 
% 12 = two Gauss*exp + Gauss*biexp, 13 = two Gauss*biexp + Gauss*exp, 
% 14 = three Gauss*biexp, 15 = three beating Gauss*biexp
numbeats = 0; % number of beating freqs allowed (0-3)
padsignal = 0; % 0 = do nothing, 1 = pad end for fitting
setpulsewidth = 2.4; % define pulse width (ps) in fit, [] = float
ltmin = [];  % define minimum lifetime (ps) in fit, [] = default (0.1)
ltmax = []; % define maximum lifetime (ps) in fit, [] = default (100)
tuningrate = 0.033; % PicoEmerald signal t0 tuning rate (ps/nm)
minprobe = []; maxprobe = []; % probe range for aligning tD, [] = auto-find
annotatefigure = 1; % 1 = annotate figure (unfinished, leave as 1)

if filetype == 3 % for image data
    writefigsyn = 0; % don't write individual figures
    if powernormyn > 2 % don't do photobleach correction
        powernormyn = 2;
    end
end

if powernormyn > 0 % prep for power normalization
    if tempmod == 1 % pre-load beamsplitter response curve if needed
        BSdata = importdata('/Users/pkocheril/Documents/Caltech/Wei Lab/Spreadsheets/BP145B2/BP145B2.csv');
    end
    if powernormyn > 2
        powernormyn = 2;
    end
end

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
    if isempty(DFGyn) == 1
        if idlerWL < 2600
            DFGyn = 1;
        else
            DFGyn = 0;
        end
    end
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
    prbWL = 515.6; % <-- unknown from filename
    idlerWL = 1E7/idlerWN;
    signalWL = 1/((1/1031.2)-(1/idlerWL));
    if IRpowermeter ~= 0
        IRpower = IRpowermeter*300; % assuming 300 mW scale
    end
    if isempty(DFGyn) == 1
        if idlerWL < 2600 % nm
            DFGyn = 1;
        else
            DFGyn = 0;
        end
    end
    if DFGyn == 1
        IRWN = DFGWN;
    else
        IRWN = idlerWN;
    end
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
    imgdim = char(fileinfo(end-4)); % '5X5'
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
    imgdim = char(fileinfo(end-4)); % '5X5'
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
    x = importdata(Tlistname); % load Tlist
    data = double(tiffreadVolume(F));
    imagesize = size(data);
    if pairedDCAC > 0
        if pairedDCAC == 2 % SPCM + CH2
            F_DC = F;
            F_DC(end-6:end+1) = "SPCM.tif";
            DCdata = double(tiffreadVolume(F_DC));
        else
            F_DC = F;
            F_DC(end-4) = "1";
            DCdata = double(tiffreadVolume(F_DC));
        end
    else % no paired DCAC
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

startiii = 1; startjjj = 1;

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
                if pairedDCAC > 0
                    DC = DC*(normprobepower/prbpower);
                    DCsds = DCsds*(normprobepower/prbpower);
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
        
        % Fitting
        if basefittype == 0 % no baseline fitting, set basecurve to 0
            basecurve = zeros(height(t),width(t));
            bleachrate = 0;
            if pairedDCAC > 0
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
                if pairedDCAC > 0 % fit DC baseline
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
                    if pairedDCAC > 0
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

        % Calculate peak height and integrated signal
        corrsig = signal - basecurve;
        % Determine if signal peak is negative or positive
        abscorrsig = abs(corrsig);
        flipcorrsig = -1*corrsig;
        corrabsdist = sum((corrsig-abscorrsig).^2);
        flipabsdist = sum((flipcorrsig-abscorrsig).^2);
        if corrabsdist < flipabsdist % peak is positive
            sigsign = 1; % corrsig is closer to abscorrsig
        else % peak is negative
            sigsign = -1; % flipcorrsig is closer to abscorrsig
            corrsig = basecurve-signal; % corrsig becomes positive
        end
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
        %figure; plot(t,corrsig);

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

        % Calculate t spacing
        deltat = zeros(length(t)-1,1); % check t spacing
        for i=1:length(t)-1
            deltat(i) = round(1e2*(t(i+1)-t(i)))/1e2; % round to nearest 0.01 ps (had issues without rounding)
        end
        tspacing = unique(deltat); % length 1 <-> evenly spaced t
        
        if isempty(ltfittype) == 1 % auto-choose lifetime fitting
            emptytype = 1;
            fittypes = importdata('/Users/pkocheril/Documents/Caltech/Wei Lab/Spreadsheets/Lifetime fit typing.csv');
            fittypes = fittypes.data;
            if prbWL > 699
                ltfittype = fittypes(round(IRWN),2);
                numbeats = fittypes(round(IRWN),3);
            else % frequency doubling
                ltfittype = fittypes(round(IRWN),4);
                numbeats = fittypes(round(IRWN),5);
            end
        else
            emptytype = 0;
        end

        if ltfittype > 0 % lifetime fitting
            % Check t spacing (even) and length (odd)
            if length(tspacing) == 1  && mod(length(t),2) == 1
                tint = t; sigint = corrsig;
            else % interpolate for convolution
                tint = linspace(t(1),t(end),length(t)*2-1).';
                sigint = spline(t, corrsig, tint);
            end

            tintspace = tint(2)-tint(1);
            originallength = length(tint);
            tlong = linspace(tint(1),tint(end)+tintspace*100,length(tint)+100).';
            siglong = sigint;
            siglong(length(sigint)+1:length(tlong)) = min(sigint(end-5:end));
            if padsignal == 1 % signal padding for more accurate fitting
                tint = tlong;
                sigint = siglong;
            end

            tc = tint(1:2:end); % odd values of tint --> convolution
            if powernormyn > 0
                IRFamp = IRpower/1e3; % IRF amplitude (~ IR power in W)
            else
                IRFamp = 50/1e3; % default to 50 mW
            end

            if ltfittype > 4
                % Define error function(r) as fit minus signal
                fun = @(r) conv(exp(-r(2)*(tc-r(1)).^2), ... % conv 1
                    heaviside(tc).*(r(3)*exp(-r(4)*(tc))+r(5)*exp(-r(6)*(tc))+...
                    r(7)*exp(-r(8)*(tc)).*cos(r(9)*tc-r(10))+ ...
                    r(11)*exp(-r(12)*(tc)).*cos(r(13)*tc-r(14))+ ...
                    r(15)*exp(-r(16)*(tc)).*cos(r(17)*tc-r(18))))+ ...
                    r(19)*tint+r(20)+r(21)*exp(-r(22)*(tint-r(23)))+ ... % baseline
                    conv(exp(-r(25)*(tc-r(24)).^2), ... % conv 2
                    heaviside(tc).*(r(26)*exp(-r(27)*(tc))+r(28)*exp(-r(29)*(tc))+...
                    r(30)*exp(-r(31)*(tc)).*cos(r(32)*tc-r(33))+ ...
                    r(34)*exp(-r(35)*(tc)).*cos(r(36)*tc-r(37))+ ...
                    r(38)*exp(-r(39)*(tc)).*cos(r(40)*tc-r(41))))+ ...
                    conv(exp(-r(43)*(tc-r(42)).^2), ... % conv 3
                    heaviside(tc).*(r(44)*exp(-r(45)*(tc))+r(46)*exp(-r(47)*(tc))+...
                    r(48)*exp(-r(49)*(tc)).*cos(r(50)*tc-r(51))+ ...
                    r(52)*exp(-r(53)*(tc)).*cos(r(54)*tc-r(55))+ ...
                    r(56)*exp(-r(57)*(tc)).*cos(r(58)*tc-r(59)))) - sigint;
                
                % Fit setup - initial guesses
                r0 = [1.3,0.05,0.02,1/2,0.02,1/10,... % first conv - center (ps), decay (ps^-2), amplitude, decay (1/ps), amp, decay (1/ps)
                    0.1,0.1,0.5,0,... % beat 1 - amplitude, decay (1/ps), freq (rad/ps), phase (rad)
                    0.1,0.1,0.5,0,... % beat 2
                    0.1,0.1,0.5,0,... % beat 3
                    0,0,0,0,0,... % baseline terms
                    16,0.07,0.03,1/2,0.003,1/10,... % second conv
                    0.1,0.1,0.5,0,...
                    0.1,0.1,0.5,0,...
                    0.1,0.1,0.5,0,...
                    60,0.2,0.002,1/5,0,1/10,... % third conv
                    0.1,0.1,0.5,0,...
                    0.1,0.1,0.5,0,...
                    0.1,0.1,0.5,0];
                lb = [-10,0.01,0,1/100,0,1/100,... % lower bounds
                    0,0,0,-180,... % default lifetime max = 100 ps
                    0,0,0,-180,... % max IRF width = 8.3 ps (0.01)
                    0,0,0,-180,...
                    0,0,0,0,0,...
                    10,0.01,0,1/100,0,1/100,... % second conv
                    0,0,0,-180,...
                    0,0,0,-180,...
                    0,0,0,-180,...
                    50,0.01,0,1/100,0,1/100,... % third conv
                    0,0,0,-180,...
                    0,0,0,-180,...
                    0,0,0,-180];
                ub = [20,1,9,1/0.1,9,1/0.1,... % upper bounds
                    9,9,9,180,... % default lifetime min = 0.1 ps
                    9,9,9,180,... % min IRF width = 0.83 ps (1)
                    9,9,9,180,...
                    9,9,9,9,9,...
                    50,1,9,1/0.1,9,1/0.1,... % second conv
                    9,9,9,180,...
                    9,9,9,180,...
                    9,9,9,180,...
                    70,1,9,1/0.1,0,1/0.1,... % third conv
                    9,9,9,180,...
                    9,9,9,180,...
                    9,9,9,180];
                
                % Pass baseline fit values to lifetime fit
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
                
                lb(19:23) = 0; ub(19:23) = 0;
                
                % Implement fit options from config
                if ltfittype == 1 % Gaussian*exp
                    lb(5:18) = 0; ub(5:18) = 0;
                end
                if ltfittype == 2 % Gaussian*biexp
                    lb(7:18) = 0; ub(7:18) = 0;
                end
                if ltfittype == 3 % Gaussian*exp with beating
                    lb(5:6) = 0; ub(5:6) = 0;
                end
                if ltfittype < 5 % single peak
                    lb(24:end) = 0; ub(24:end) = 0;
                end
                if ltfittype == 5 % two Gauss*exp
                    lb([5:18 28:41]) = 0; ub([5:18 28:41]) = 0;
                end
                if ltfittype == 6 % Gauss*exp + Gauss*biexp
                    lb([5:18 30:41]) = 0; ub([5:18 30:41]) = 0;
                end
                if ltfittype == 7 % two Gauss*biexp
                    lb([7:18 30:41]) = 0; ub([7:18 30:41]) = 0;
                end
                if ltfittype == 8 % two beating Gauss*exp
                    lb([5:6 28:29]) = 0; ub([5:6 28:29]) = 0;
                end
                if ltfittype == 9 % beating Gauss*exp + Gauss*biexp
                    lb(5:6) = 0; ub(5:6) = 0;
                end
                if ltfittype < 11 % two peaks
                    lb(42:end) = 0; ub(42:end) = 0;
                end
                if ltfittype == 11 % three Gauss*exp
                    lb([5:18 28:41 46:59]) = 0; ub([5:18 28:41 46:59]) = 0;
                end
                if ltfittype == 12 % two Gauss*exp + Gauss*biexp (middle)
                    lb([5:18 30:41 46:59]) = 0; ub([5:18 30:41 46:59]) = 0;
                end
                if ltfittype == 13 % two Gauss*biexp + Gauss*exp (middle)
                    lb([7:18 28:41 48:59]) = 0; ub([7:18 28:41 48:59]) = 0;
                end
                if ltfittype == 14 % three Gauss*biexp
                    lb([7:18 30:41 48:59]) = 0; ub([7:18 30:41 48:59]) = 0;
                end
                if numbeats > 0 % beating control
                    if numbeats == 1
                        lb([11:18 34:41 52:59]) = 0; ub([11:18 34:41 52:59]) = 0;
                    end
                    if numbeats == 2
                        lb([15:18 38:41 56:59]) = 0; ub([15:18 38:41 56:59]) = 0;
                    end
                else % no beating
                    lb([7:18 30:41 48:59]) = 0; ub([7:18 30:41 48:59]) = 0;
                end

                % Parameter constraints from config
                if isempty(setpulsewidth) == 0 % pulse width was specified
                    pulsedecay = log(2)/setpulsewidth^2;
                    lb(2) = pulsedecay; ub(2) = pulsedecay;
                end
                if isempty(ltmin) == 0 % minimum lifetime was specified
                    ub(4) = 1/ltmin; ub(6) = ub(4); % min lt -> ub decay
                    ub(27) = ub(4); ub(29) = ub(4);
                    ub(45) = ub(4); ub(47) = ub(4);
                end
                if isempty(ltmax) == 0 % maximum lifetime was specified
                    lb(4) = 1/ltmax; lb(6) = lb(4); % max lt -> lb decay
                    lb(27) = lb(4); lb(29) = lb(4);
                    lb(45) = lb(4); lb(47) = lb(4);
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
            
                % Generate curve for fitted parameters
                fitcurve = conv(exp(-fitval(2)*(tc-fitval(1)).^2),... % first conv
                    heaviside(tc).*(fitval(3)*exp(-fitval(4)*(tc))+ ...
                    fitval(5)*exp(-fitval(6)*(tc))+ ...
                    fitval(7)*exp(-fitval(8)*(tc)).* ...
                    cos(fitval(9)*tc-fitval(10))+ ...
                    fitval(11)*exp(-fitval(12)*(tc)).* ...
                    cos(fitval(13)*tc-fitval(14))+ ...
                    fitval(15)*exp(-fitval(16)*(tc)).* ...
                    cos(fitval(17)*tc-fitval(18))))+...
                    fitval(19)*tint+fitval(20)+fitval(21)*... % baseline terms
                    exp(-fitval(22)*(tint-fitval(23)))+...
                    conv(exp(-fitval(25)*(tc-fitval(24)).^2),... % second conv
                    heaviside(tc).*(fitval(26)*exp(-fitval(27)*(tc))+ ...
                    fitval(28)*exp(-fitval(29)*(tc))+ ...
                    fitval(30)*exp(-fitval(31)*(tc)).* ...
                    cos(fitval(32)*tc-fitval(33))+ ...
                    fitval(34)*exp(-fitval(35)*(tc)).* ...
                    cos(fitval(36)*tc-fitval(37))+ ...
                    fitval(38)*exp(-fitval(39)*(tc)).* ...
                    cos(fitval(40)*tc-fitval(41))))+...
                    conv(exp(-fitval(43)*(tc-fitval(42)).^2),... % third conv
                    heaviside(tc).*(fitval(44)*exp(-fitval(45)*(tc))+ ...
                    fitval(46)*exp(-fitval(47)*(tc))+ ...
                    fitval(48)*exp(-fitval(49)*(tc)).* ...
                    cos(fitval(50)*tc-fitval(51))+ ...
                    fitval(52)*exp(-fitval(53)*(tc)).* ...
                    cos(fitval(54)*tc-fitval(55))+ ...
                    fitval(56)*exp(-fitval(57)*(tc)).* ...
                    cos(fitval(58)*tc-fitval(59))));
                
                % Lifetime assignment - want τ3 < τ4, τ5 < τ6
                lifetime3 = min(abs(1/fitval(27)),abs(1/fitval(29)));
                lifetime4 = max(abs(1/fitval(27)),abs(1/fitval(29)));
                lifetime5 = min(abs(1/fitval(45)),abs(1/fitval(47)));
                lifetime6 = max(abs(1/fitval(45)),abs(1/fitval(47)));
                if lifetime3 == Inf % possible from fit type options
                    lifetime3 = 0;
                end
                if lifetime4 == Inf
                    lifetime4 = 0;
                end
                if lifetime5 == Inf
                    lifetime5 = 0;
                end
                if lifetime6 == Inf
                    lifetime6 = 0;
                end

                % Decay ratios
                if abs(1/fitval(27)) < abs(1/fitval(29))
                    lt3lt4 = fitval(26)/fitval(28); % A3/A4
                else % means 1/fitval(27) is longer
                    lt3lt4 = fitval(28)/fitval(26); % still A3/A4
                end
                if abs(1/fitval(45)) < abs(1/fitval(47))
                    lt5lt6 = fitval(44)/fitval(46); % A5/A6
                else % means 1/fitval(45) is longer
                    lt5lt6 = fitval(46)/fitval(44); % still A5/A6
                end
            else % single-peak fitting
                % Define error function(r) as fit minus signal
                fun = @(r) conv(exp(-r(2)*(tc-r(1)).^2), ...
                    heaviside(tc).*(r(3)*exp(-r(4)*(tc))+r(5)*exp(-r(6)*(tc))+...
                    r(7)*exp(-r(8)*(tc)).*cos(r(9)*tc-r(10))+ ...
                    r(11)*exp(-r(12)*(tc)).*cos(r(13)*tc-r(14))+ ...
                    r(15)*exp(-r(16)*(tc)).*cos(r(17)*tc-r(18))))+ ...
                    r(19)*tint+r(20)+r(21)*exp(-r(22)*(tint-r(23))) - sigint;
                
                % Fit setup - initial guesses
                r0 = [1.3,0.05,0.02,1/2,0.02,1/10,... % first conv - center (ps), decay (ps^-2), amplitude, decay (1/ps), amp, decay (1/ps)
                    0.1,0.1,0.5,0,... % beat 1 - amplitude, decay (1/ps), freq (rad/ps), phase (rad)
                    0.1,0.1,0.5,0,... % beat 2
                    0.1,0.1,0.5,0,... % beat 3
                    0,0,0,0,0];... % baseline terms
                lb = [-10,0.01,0,1/100,0,1/100,... % lower bounds
                    0,0,0,-180,... % default lifetime max = 100 ps
                    0,0,0,-180,... % max IRF width = 8.3 ps (0.01)
                    0,0,0,-180,...
                    0,0,0,0,0];
                ub = [20,1,9,1/0.1,9,1/0.1,... % upper bounds
                    9,9,9,180,... % default lifetime min = 0.1 ps
                    0,9,9,180,... % min IRF width = 0.83 ps (1)
                    0,9,9,180,...
                    9,9,9,9,9];
                
                % Pass baseline fit values to lifetime fit
                if basefittype == 1 % linear baseline
                    lb(21:23) = 0; ub(21:23) = 0; % mute exponential terms
                    lb(19) = basecoef(1); ub(19) = lb(19);
                    lb(20) = basecoef(2); ub(20) = lb(20);
                    r0(19) = basecoef(1); r0(20) = basecoef(2);
                end
                if basefittype == 2 % exponential baseline
                    lb(19) = 0; ub(19) = 0; % mute linear term
                    lb(20) = basefit(1); ub(20) = lb(20);
                    lb(21) = basefit(2); ub(21) = lb(21);
                    lb(22) = basefit(3); ub(22) = lb(22);
                    lb(23) = basefit(4); ub(23) = lb(23);
                    r0(20) = basefit(1); r0(21) = basefit(2);
                    r0(22) = basefit(3); r0(23) = basefit(4);
                end
                
                % Baseline removed
                lb(19:23) = 0; ub(19:23) = 0;

                % Implement fit options from config
                if ltfittype == 1 % Gaussian*exp
                    lb(5:18) = 0; ub(5:18) = 0;
                end
                if ltfittype == 2 % Gaussian*biexp
                    lb(7:18) = 0; ub(7:18) = 0;
                end
                if ltfittype == 3 % Gaussian*exp with beating
                    lb(5:6) = 0; ub(5:6) = 0;
                end
                lb(24:end) = 0; ub(24:end) = 0; % single peak
                if numbeats > 0 % beating control
                    if numbeats == 1
                        lb(11:18) = 0; ub(11:18) = 0;
                    end
                    if numbeats == 2
                        lb(15:18) = 0; ub(15:18) = 0;
                    end
                else % no beating
                    lb(7:18) = 0; ub(7:18) = 0;
                end

                % Parameter constraints from config
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
                options.OptimalityTolerance = 1e-12; % make the fit more accurate
                options.StepTolerance = 1e-12; % make the fit more precise
                options.Display = 'off'; % silence console output
                fitval = lsqnonlin(fun,r0,lb,ub,options); % fit coeffs
            
                % Generate curve for fitted parameters
                fitcurve = conv(exp(-fitval(2)*(tc-fitval(1)).^2),... % first conv
                    heaviside(tc).*(fitval(3)*exp(-fitval(4)*(tc))+ ...
                    fitval(5)*exp(-fitval(6)*(tc))+ ...
                    fitval(7)*exp(-fitval(8)*(tc)).* ...
                    cos(fitval(9)*tc-fitval(10))+ ...
                    fitval(11)*exp(-fitval(12)*(tc)).* ...
                    cos(fitval(13)*tc-fitval(14))+ ...
                    fitval(15)*exp(-fitval(16)*(tc)).* ...
                    cos(fitval(17)*tc-fitval(18))))+...
                    fitval(19)*tint+fitval(20)+fitval(21)*... % baseline terms
                    exp(-fitval(22)*(tint-fitval(23)));

                lifetime3 = 0; lifetime4 = 0; lt3lt4 = 0;
                lifetime5 = 0; lifetime6 = 0; lt5lt6 = 0;
            end
            
            % Lifetime assignment - want τ1 < τ2
            lifetime1 = min(abs(1/fitval(4)),abs(1/fitval(6)));
            lifetime2 = max(abs(1/fitval(4)),abs(1/fitval(6)));
            if lifetime1 == Inf % possible from fit type options
                lifetime1 = 0;
            end
            if lifetime2 == Inf
                lifetime2 = 0;
            end
            % Decay ratio
            if abs(1/fitval(4)) < abs(1/fitval(6))
                lt1lt2 = fitval(3)/fitval(5); % A1/A2
            else % means 1/fitval(4) is longer
                lt1lt2 = fitval(5)/fitval(3); % still A1/A2
            end

            % Return to original length and sign
            %figure; plot(tint,sigint,tint,fitcurve);
            tint = tint(1:originallength);
            sigint = sigint(1:originallength);
            fitcurve = fitcurve(1:originallength);
            sigint = basecurve+sigsign*sigint(1:originallength);
            fitcurve = basecurve+sigsign*fitcurve(1:originallength);
            %figure; plot(tint,sigint);

            % Calculate residuals, ssresid, and IRF FWHM
            resid = sigint-fitcurve;
            ssresid = sum(resid.^2); % sum of squares of residuals
            tss = sum((sigint-mean(sigint)).^2); % total sum of squares
            r2 = 1-(ssresid/tss); % coefficient of determination (R^2)
            if r2 < 0 % remove nonsense R^2 values
                r2 = 0;
            end
            FWHM = sqrt(log(2)/fitval(2)); % pulse width, ps
        else % no lifetime fitting
            tint = t; sigint = signal; fitcurve = signal;
            lifetime1 = 0; lifetime2 = 0; lt1lt2 = 0;
            lifetime3 = 0; lifetime4 = 0; lt3lt4 = 0;
            lifetime5 = 0; lifetime6 = 0; lt5lt6 = 0;
            fitval = zeros(59,1); r2 = 0; FWHM = 0;
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
        if numbeats == 3 % make beating annotation
            beatstring = string(fitval(9)*1e12/(pi*29979245800))+...
                ', '+string(fitval(13)*1e12/(pi*29979245800))+...
                ', '+string(fitval(17)*1e12/(pi*29979245800));
        else
            if numbeats == 2
                beatstring = string(fitval(9)*1e12/(pi*29979245800))+...
                    ', '+string(fitval(13)*1e12/(pi*29979245800));
            else
                if numbeats == 1
                    beatstring = string(fitval(9)*1e12/(pi*29979245800));
                else
                    beatstring = "None";
                end
            end
        end
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
                tiledlayout(2,2);
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
    end
end
if writefigsyn == 1 & filetype == 3 % write processed image tif
    Save_tiff(outname+'_proc.tif',lt1array,...
        lt2array,ssresidarray,FWHMarray,r2array)
end

conv1curve = basecurve+conv(exp(-fitval(2)*(tc-fitval(1)).^2),... % first conv
    heaviside(tc).*(fitval(3)*exp(-fitval(4)*(tc))+ ...
    fitval(5)*exp(-fitval(6)*(tc))+ ...
    fitval(7)*exp(-fitval(8)*(tc)).* ...
    cos(fitval(9)*tc-fitval(10))+ ...
    fitval(11)*exp(-fitval(12)*(tc)).* ...
    cos(fitval(13)*tc-fitval(14))+ ...
    fitval(15)*exp(-fitval(16)*(tc)).* ...
    cos(fitval(17)*tc-fitval(18))));
conv2curve = basecurve+conv(exp(-fitval(25)*(tc-fitval(24)).^2),... % second conv
    heaviside(tc).*(fitval(26)*exp(-fitval(27)*(tc))+ ...
    fitval(28)*exp(-fitval(29)*(tc))+ ...
    fitval(30)*exp(-fitval(31)*(tc)).* ...
    cos(fitval(32)*tc-fitval(33))+ ...
    fitval(34)*exp(-fitval(35)*(tc)).* ...
    cos(fitval(36)*tc-fitval(37))+ ...
    fitval(38)*exp(-fitval(39)*(tc)).* ...
    cos(fitval(40)*tc-fitval(41))));

legend2 = legend1;
legend2(end+1) = {'Curve 1'};
legend2(end+1) = {'Curve 2'};

figure;
hold on; errorbar(t,signal,sigsds,'o')
plot(tbase,sigbase,'go',t,basecurve,'-'); 
plot(tint,fitcurve,'-','LineWidth',2);
plot(tint,conv1curve,'--','LineWidth',2);
plot(tint,conv2curve,'r--','LineWidth',2);
hold off;
xlabel('Time (ps)'); ylabel('AC signal (AU)');
legend(legend2);

% % FFT of residuals
% dtp = abs(tint(1)-tint(2)); % time step, ps
% fsp = 1/dtp; % freq BW, ps-1
% ts = tint*1e-12; % time, s
% dt = abs(ts(1)-ts(2)); % time step, s
% fs = 1/dt; % freq BW, Hz
% n = length(resid);
% 
% y = fft(resid); % raw FFT
% yshift = fftshift(y); % centered FFT
% power = abs(yshift).^2/n; % power spectrum
% f = (-n/2:n/2-1)*(fs/n); % frequency vector in Hz
% 
% wn = f/29979245800; % freq in cm-1
% 
% startfit = 23; endfit = 27;
% fftfit = fit(wn(startfit:endfit).',power(startfit:endfit),'gauss1');
% fftcoeffs = coeffvalues(fftfit);
% beatwn = fftcoeffs(2);
% beatannot = 'Beat = '+string(beatwn)+' cm^{-1}';
% 
% figure;
% hold on;
% plot(wn,power);
% %plot(fftfit);
% hold off;
% xlabel('Frequency (cm-1)')
% xlim([0 max(wn)]) % for cm-1
% ylabel('FFT power spectrum')
% %xlabel('Frequency (Hz)')
% %xlim([0 3e11]) % for Hz
% %legend('FFT of residuals','Gaussian fit');
% %annotation('textbox',[0.5 0.5 0.2 0.2],'String',beatannot,...
% %    'FitBoxToText','on');
% 
% xvector = [1 2 4 1 2 4 1 2 4].';
% yvector = [3.5 6.4 13.4 3.5 6.8 13.9 2.9 6.4 13.4].';
% 
% beatfit = fitlm(xvector,yvector);
% 
% figure; plot(beatfit)



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
