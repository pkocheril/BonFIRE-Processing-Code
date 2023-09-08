%% Fit temporal sweep data (second generation)
%%% -- replaced by sweep_proc_alignv4

% Load data
cd '/Users/pkocheril/Documents/Caltech/Wei Lab/Data/MATLAB/'
clear; close all; clc;

F = '1590cm-1_750 (15)_5X5_Idler4053.8_DFG1589.8_P0.1981_CH2.tif';
Tlistname = 'Tlist_1590.txt';

% Configuration options - check before running!
visibility = 'on';
writefigsyn = 0; % 1 = write figure files, 0 = not

filetype = 2; % 1 = .txt, 2 = .tif (solution), 3 = .tif (image)
filenameconv = 2; % 0 = other,
% 1 = [probeWL]-[probepower]-[signalWL]-[IRpower]-[ND]-[PMTgain]-[modfreq]-[PMTBW]-[etc],
% 2 = [etc]_[size]_[idler]_[DFG]_[power]_[channel],
% 3 = [IRWN]_[probeWL]_[etc]_[size]_[idler]_[DFG]_[power]_[channel],
% 4 = [probeWL]_[IRWN]_[etc]_[size]_[idler]_[DFG]_[power]_[channel]
t0pos = []; % specify t0 position (mm), [] = autofind
pairedDCAC = 0; % 2 = SPCM + PMT CH2, 1 = PMT CH1 + CH2, 0 = only PMT CH2
powernormyn = 0; % 0 = no normalization, 1 = normalize by IR power,
% 2 = normalize by probe and IR powers
tempmod = 0; % 1 = temporal modulation in place, 0 = no beamsplitters
ND = 2; % ND filter strength in probe path (script will update if possible)
PMT = 1; % PMT gain (script will update if possible)
prbpowerset = 250; % probe power in mW (script will update if possible)
IRpowerset = 70; % IR power in mW (script will update if possible)
normIRpower = 50; % IR power on-sample to normalize to (default is 50)
normprobepower = 1.5; % probe power on-sample to normalize to (default 1.5)
trimlastT = 0; % how many points to remove from end of Tlist (default 0)
trimfirstT = 0; % points to remove from start of Tlist (default 0)
basefittype = 1; % 0 = no baseline fit, 1 = linear fit, 2 = exp, 3 = exp + line
floatbase = 0; % percentage to float baseline coeffs (10 -> +/- 10%)
cutlow = 0.05; % lower baseline fit cutoff, (default 10th %ile)
cuthigh = 0.95; % upper baseline fit cutoff, (default 90th %ile)
ltfittype = 2; % [] = auto-choose, 0 = no fitting, 1 = Gaussian*monoexp,
% 2 = Gaussian*biexp, 3 = Gaussian*stretchexp
setpulsewidth = []; % define pulse width (ps) in fit, [] = float
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
    imgdim = char(fileinfo(end-4)); % '5X5'
    prbWL = 515.6; % <-- unknown from filename
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
if filenameconv == 3 % parse probe sweep .tifs
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
    if x(1) == x(end) && trimlastT == 0
        trimlastT = 1;
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
    if x(1) == x(end) && trimlastT == 0
        trimlastT = 1;
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

        %signal = -signal+2;

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
            % Check baseline fitting
            %figure; plot(t,signal,'o',t,basecurve,'-');
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
                tint = t; sigint = signal; % data are good for fitting
            else % interpolate for convolution
                tint = linspace(t(1),t(end),length(t)*2-1).';
                sigint = spline(t, signal, tint);
            end

            if max(tint) < 30 % auto-pad signal if tint doesn't go past 30 ps
                padsignal = 1;
            else
                padsignal = 0;
            end

            % Signal padding for short-tailed data
            tintspace = tint(2)-tint(1);
            originallength = length(tint);
            tlong = linspace(tint(1),tint(end)+tintspace*100,length(tint)+100).';
            siglong = sigint;
            siglong(length(sigint)+1:length(tlong)) = min(sigint(end-5:end));
            if padsignal == 1
                tint = tlong;
                sigint = siglong;
            end

            tc = tint(1:2:end); % odd values of tint --> convolution input
            if powernormyn > 0
                IRFamp = IRpower/1e3; % IRF amplitude (~ IR power in W)
            else
                IRFamp = 50/1e3; % default to 50 mW
            end

            fun = @(r) r(1)+r(2)*exp(-r(3)*(tint-r(4)))+r(5)*tint+ ... % baseline
                conv(exp(-((tc-r(6))./(r(7)/(2*sqrt(log(2))))).^2), ...
                IRFamp*heaviside(tc).*(r(8)*exp(-(tc)./r(9))+r(10)*exp(-((tc)./(r(11))).^r(12)))) - sigint;
            
            % Initial guesses
            r0 = [basefit,0,4,... % basecoefs, IRF center (ps), IRF width (ps)
                0.02,2,0.02,10,1]; % amp 1, τ1 (ps), amp 2, τ2 (ps), β (stretch)

            if floatbase > 0 % float baseline coeffs
                newlbbase = lbbase; newubbase = ubbase;
                for i=1:length(lbbase)
                    if lbbase(i) == 0
                        newlbbase(i) = 0;
                    else
                        newlbbase(i) = min((1-floatbase/100)*basefit(i),(1+floatbase/100)*basefit(i));
                    end
                    if ubbase(i) == 0
                        newubbase(i) = 0;
                    else
                        newubbase(i) = max((1-floatbase/100)*basefit(i),(1+floatbase/100)*basefit(i));
                    end
                end
            else
                newlbbase = basefit; newubbase = basefit;
            end

            % Lower and upper bounds
            lb = [newlbbase,-30,2,... % basecoefs, IRF center (ps), IRF min (ps)
                -9,0.1,... % amp 1, τ1min (ps)
                -9,0.1,... % amp 2, τ2min (ps)
                1e-2]; % β min
            ub = [newubbase,30,9,... % basecoefs, IRF center (ps), IRF max (ps)
                9,100,... % amp 1, τ1max (ps)
                9,100,... % amp 2, τ2max (ps)
                1e2]; % β max

            % Implement fit type options
            if ltfittype == 1 % Gauss*monoexp
                lb(10) = 0; ub(10) = 0;
            end
            if ltfittype == 2 % Gauss*biexp
                lb(12) = 1; ub(12) = 1;
            end
            if ltfittype == 3 % Gauss*stretchexp
                lb(8) = 0; ub(8) = 0;
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
            
            IRFwidth = fitval(7); % ps
            FWHM = IRFwidth/sqrt(2); % ps
            lifetime1 = min(fitval(9),fitval(11)); % ps
            lifetime2 = max(fitval(9),fitval(11)); % ps
            % Lifetime amplitude ratio
            if abs(fitval(9)) < abs(fitval(11))
                lt1lt2 = fitval(8)/fitval(10); % A1/A2
            else % means fitval(9) is longer (τ2)
                lt1lt2 = fitval(10)/fitval(8); % still A1/A2
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

            % Check fit
            %figure; tiledlayout(3,1); nexttile([1 1]); plot(tint,resid); nexttile([2 1]); plot(t,signal,'o',tint,fitcurve);
            %figure; plot(wn,powerfft);
        else % no lifetime fitting
            tint = t; sigint = signal; fitcurve = signal;
            lifetime1 = 0; lifetime2 = 0; lt1lt2 = 0;
            fitval = zeros(12,1); r2 = 0; FWHM = 0;
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
    end
end
if writefigsyn == 1 & filetype == 3 % write processed image tif
    Save_tiff(outname+'_proc.tif',lt1array,...
        lt2array,ssresidarray,FWHMarray,r2array)
end

% % Red-white-blue gradient for contour map
% fineness = 100;
% map1 = [linspace(0,1,fineness) ones(1,fineness-1)];
% map3 = flip(map1);
% map2 = min(map1,map3);
% map = [map1.' map2.' map3.'];
% 
% figure; % contour plot for w
% contourf(conc,rad,w)
% xlabel('c (M)')
% ylabel('r (nm)')
% set(gca,'xscale','log')
% xticks([1E-3 1E-2 1E-1 1E0]);
% colormap(map);
% cb = colorbar;

% % Plotting multiple convolutions
% conv1curve = basecurve+conv(exp(-fitval(2)*(tc-fitval(1)).^2),... % first conv
%     heaviside(tc).*(fitval(3)*exp(-fitval(4)*(tc))+ ...
%     fitval(5)*exp(-fitval(6)*(tc))+ ...
%     fitval(7)*exp(-fitval(8)*(tc)).* ...
%     cos(fitval(9)*tc-fitval(10))+ ...
%     fitval(11)*exp(-fitval(12)*(tc)).* ...
%     cos(fitval(13)*tc-fitval(14))+ ...
%     fitval(15)*exp(-fitval(16)*(tc)).* ...
%     cos(fitval(17)*tc-fitval(18))));
% conv2curve = basecurve+conv(exp(-fitval(25)*(tc-fitval(24)).^2),... % second conv
%     heaviside(tc).*(fitval(26)*exp(-fitval(27)*(tc))+ ...
%     fitval(28)*exp(-fitval(29)*(tc))+ ...
%     fitval(30)*exp(-fitval(31)*(tc)).* ...
%     cos(fitval(32)*tc-fitval(33))+ ...
%     fitval(34)*exp(-fitval(35)*(tc)).* ...
%     cos(fitval(36)*tc-fitval(37))+ ...
%     fitval(38)*exp(-fitval(39)*(tc)).* ...
%     cos(fitval(40)*tc-fitval(41))));
% 
% legend2 = legend1;
% legend2(end+1) = {'Curve 1'};
% legend2(end+1) = {'Curve 2'};
% 
% figure;
% hold on; errorbar(t,signal,sigsds,'o')
% plot(tbase,sigbase,'go',t,basecurve,'-'); 
% plot(tint,fitcurve,'-','LineWidth',2);
% plot(tint,conv1curve,'--','LineWidth',2);
% plot(tint,conv2curve,'r--','LineWidth',2);
% hold off;
% xlabel('Time (ps)'); ylabel('AC signal (AU)');
% legend(legend2);



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
