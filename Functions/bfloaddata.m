% Quick load BonFIRE data
function [T,SIGNAL,INFO,DATA,SIGSDS,DELAYPOS] = bfloaddata(varargin)
%%% This function loads BonFIRE data from a given input file, meant to be
% a lighter, easier alternative to "loaddata"

%%%%%%%% needs some work for if info file isn't there

    % Use input parser and varargin to enable calling as "bfloaddata()"
    defaultFile = 'nofile';
    p = inputParser;
    addOptional(p,'filename',defaultFile,@isfile);
    parse(p,varargin{:})

    CURRENTFILE = p.Results.filename;

    % If file doesn't exist or no file specified, use file browser GUI
    if ~exist(CURRENTFILE,'file') || isempty(CURRENTFILE) || ~isfile(CURRENTFILE)
        [filename,folder] = uigetfile('*.*');
        CURRENTFILE = fullfile(folder,filename);
    end

    % Filename parsing
    macseparator = count(CURRENTFILE,'/');
    winseparator = count(CURRENTFILE,'\');
    if macseparator > winseparator % MacOS/UNIX
        delimiter = "/";
    else % Windows
        delimiter = "\";
    end
    folderinfo = split(CURRENTFILE,delimiter); % split path into folders
    filename = char(folderinfo(end));

    % Figure out solution vs image data by filename
    if ~contains(CURRENTFILE,'FOV') % filename doesn't have 'FOV'
        currdatatype = 0; % solution data
    else % filename has 'FOV'
        currdatatype = 1; % image data
    end

    % Figure out file naming convention
    if contains(CURRENTFILE,'.txt') %count(currentfile,'.txt') > 0
        lookforhyphens = count(filename,'-');
        if numel(lookforhyphens) > 5
            currfilenameconv = 1; % my old .txts
        else
            currfilenameconv = 0; % even older .txts
        end
    else % .tif or .raw
        if strcmp(filename(5:10),'cm-1_7') && contains(filename,' (') % filename ####cm-1_750 (##)...
            currfilenameconv = 3; % manual probe sweep (before updated VI)
        else
            currfilenameconv = 2; % not a manual probe sweep
        end
    end

    % Parse info from filename
    if currfilenameconv == 1 % parse my old .txts
        fileinfo = split(filename,"-"); % split filename at "-"s
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
        filebase = strjoin(fileinfo(1:end-(numel(fileinfo)-1)),'-');
    else % image file (.raw or .tif)
        fileinfo = split(folderinfo(end),"_"); % split at "_"s
        % channel = char(fileinfo(end)); % 'CH1','CH2','SPCM','CMOS'
        if contains(filename,'DFG') && contains(filename,'Idler') % numel(fileinfo) > 5 % Haomin's VI naming
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
            filebase = strjoin(fileinfo(1:end-(numel(fileinfo)-1)),'_');
        end
    end

    % Guess matching info file name
    matchIRsweep = folderinfo; matchIRsweep(end) = {strjoin({filebase,'IRsweep.txt'},'_')}; matchIRsweep = strjoin(matchIRsweep,delimiter);
    matchXYZTsweep = folderinfo; matchXYZTsweep(end) = {strjoin({filebase,'XYZT.txt'},'_')}; matchXYZTsweep = strjoin(matchXYZTsweep,delimiter);
    matchZsweep = folderinfo; matchZsweep(end) = {strjoin({filebase,'Z.txt'},'_')}; matchZsweep = strjoin(matchZsweep,delimiter);
    matchTsweep = folderinfo; matchTsweep(end) = {strjoin({filebase,'T.txt'},'_')}; matchTsweep = strjoin(matchTsweep,delimiter);
    matchTXYZsweep = folderinfo; matchTXYZsweep(end) = {strjoin({filebase,'XYZT.txt'},'_')}; matchTXYZsweep = strjoin(matchTXYZsweep,delimiter);
    
    % Look for matching info file
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
                    if isfile(matchTXYZsweep)
                        infofilename = matchTXYZsweep;
                    else
                        % warning('No matching .txt infofile found')
                        infofilename = [];
                    end
                end
            end
        end
    end
    
    % If info file found
    if ~isempty(infofilename) % parse info from .txt files
        infotable = readtable(infofilename); infofile = table2array(infotable(3:end,2));
        if length(infofile) < 12 % import failed
        	infoarray = importdata(infofilename); infofile = infoarray.data;
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
                prbWL = infofile(20); % manual probe WL, nm
            end
            if round(infofile(21)) ~= 0
                prbpower = infofile(21); % manual probe power, mW
            end
            if length(infofile) > 21 % if using updated code
                pmtgain = infofile(22); % PMT gain
                pmtBW = 1e3*infofile(23); % PMT bandwidth, Hz (kHz in file)
                modfreq = 1e3*infofile(24); % IR modulation frequency, Hz (kHz in file)
                pinholeyn = infofile(25); % 0 = bypassed, 1 = in place
                prbdichroic = filterlookup(infofile(26));
                pmtcmosmirror = filterlookup(infofile(27));
                pmtfilter = filterlookup(infofile(28));
                spcmfilter = filterlookup(infofile(29));
                cmosfilter = filterlookup(infofile(30));
            end
            if length(infofile) > 30
                dutycycle = infofile(31);
            end
            if length(infofile) > 31
                probeFWHM = infofile(32); % nm
                levanteFWHM = infofile(33); % nm
            end
            if length(infofile) > 39 % new unified VI
                timeconstant = infofile(34); % ms
                wheelfilter = filterlookup(infofile(35));
                concentration = infofile(36);
                dlstime = infofile(37); % ms
                piezotime = infofile(38); % ms
                TCorder = infofile(39);
                levanteraw = table2array(infotable(end,2:end));
                % levantespectrum(:,1) = levanteraw(1:2:end-1).';
                % levantespectrum(:,2) = levanteraw(2:2:end).';
            end
        else % if not extended info
            pmtgain = 1; pmtBW = 1; modfreq = 0; pinholeyn = 0;
            prbdichroic = 0; pmtcmosmirror = 0; pmtfilter = 0;
            spcmfilter = 0; cmosfilter = 0; dutycycle = 0;
            probeFWHM = 0; levanteFWHM = 0; timeconstant = 0;
            wheelfilter = 0; concentration = 0; dlstime = 0;
            piezotime = 0; TCorder = 0; levanteraw = 0; 
            % levantespectrum = 0;
            % if verbose == 2
            %     verbose = 1; % can't print out experimental parameters
            % end
        end
    % else % no info file found
    %     xsteps = 1; ysteps = 1; % if isempty?
    end

    % Summarize filters
    bandpass = 'PMT'+string(pmtfilter)+'| SPCM'+string(spcmfilter)+'| CMOS'+string(cmosfilter)+'| Wheel'+string(wheelfilter);

    % Assign IRWN by idler wavelength
    idlerWL = 1e7/idlerWN; % signalWL = 1/((1/1031.2)-(1/idlerWL));
    if idlerWL < 2600 % nm
        IRWN = DFGWN;
    else
        IRWN = idlerWN;
    end
    
    % Calculate wavelength for manual probe sweep
    if currfilenameconv == 3
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

    % Write experimental parameters to current structure
    INFO.parameters.IRWN = IRWN; INFO.parameters.IRpower = IRpower; INFO.parameters.dwelltime = dwelltime;
    INFO.parameters.prbWL = prbWL; INFO.parameters.prbpower = prbpower; INFO.parameters.ND = ND;
    INFO.parameters.xsteps = xsteps; INFO.parameters.ysteps = ysteps; INFO.parameters.zsteps = zsteps; INFO.parameters.tsteps = tsteps;
    INFO.parameters.xinitial = xinitial; INFO.parameters.yinitial = yinitial; INFO.parameters.zinitial = zinitial; INFO.parameters.tinitial = tinitial;
    INFO.parameters.xfinal = xfinal; INFO.parameters.yfinal = yfinal; INFO.parameters.zfinal = zfinal; INFO.parameters.tfinal = tfinal;
    INFO.parameters.pmtgain = pmtgain; INFO.parameters.pmtBW = pmtBW; INFO.parameters.modfreq = modfreq; INFO.parameters.pinholeyn = pinholeyn;
    INFO.parameters.prbdichroic = prbdichroic; INFO.parameters.pmtcmosmirror = pmtcmosmirror; INFO.parameters.bandpass = bandpass;
    INFO.parameters.pmtfilter = pmtfilter; INFO.parameters.spcmfilter = spcmfilter; INFO.parameters.cmosfilter = cmosfilter;
    INFO.parameters.dutycycle = dutycycle; INFO.parameters.probeFWHM = probeFWHM; INFO.parameters.levanteFWHM = levanteFWHM;
    INFO.parameters.timeconstant = timeconstant; INFO.parameters.concentration = concentration; INFO.parameters.dlstime = dlstime;
    INFO.parameters.piezotime = piezotime; INFO.parameters.TCorder = TCorder; INFO.parameters.levanteraw = levanteraw;

    % Guess Tlist name (per-file)
    guessTlist = folderinfo; guessTlist(end) = {strjoin({filebase,'Tlist.txt'},'_')}; guessTlistname = strjoin(guessTlist,delimiter);

    % Guess file type from filename
    if contains(CURRENTFILE,'.txt') % solution .txt
        DATA = importdata(CURRENTFILE); % load data
        DELAYPOS = DATA(:,1); % save delay position
        % CH1DATA = DATA(:,2);
        SIGNAL = DATA(:,3);
        SIGSDS = zeros(length(SIGNAL));
    else % raw or tif
        % Load Tlist
        if isfile(guessTlistname)
            DELAYPOS = importdata(guessTlistname);
        else
            if isfile(fullfile(folder,'Tlist.txt'))
                DELAYPOS = importdata(fullfile(folder,'Tlist.txt'));
            else
                warning('No Tlist found - guessing Tlist')
                DELAYPOS = linspace(1,4,61);
            end
        end

        % Load data
        if contains(CURRENTFILE,'.raw') % .raw
            imgid = fopen(CURRENTFILE);
            imgvector = fread(imgid,'real*8'); % load raw file as a vector
            fclose('all');
            matrixarea = length(imgvector)/length(DELAYPOS);
            if matrixarea == xsteps*ysteps
                % data match expected size - do nothing
            else
                if mod(sqrt(matrixarea),1) == 0 % if data are a square
                    xsteps = sqrt(matrixarea); ysteps = sqrt(matrixarea);
                else
                    if currdatatype == 1 % if image data
                        if trueTlist == 1 % and real Tlist file found
                            warning('Mismatch between dimensions of .raw file and filename.')
                        else % generated Tlist
                            warning('Mismatch in dimensions of file - check Tlist specification and file dimensions.')
                        end
                    end
                end
            end
            tempdata = reshape(imgvector,[ysteps,xsteps,length(DELAYPOS)]);
            DATA = zeros(xsteps,ysteps,length(DELAYPOS));
            for i=1:length(DELAYPOS)
                DATA(:,:,i) = tempdata(:,:,i).'; % flip image to match Fiji
            end
        else % tif
            DATA = double(tiffreadVolume(CURRENTFILE));
        end

        % Guess data type from filename
        if contains(CURRENTFILE,'FOV') % image
            SIGNAL = DATA;
            SIGSDS = [];
        else % solution - average points together
            SIGNAL = squeeze(mean(DATA, [1 2]));
            SIGSDS = squeeze(std(DATA, 0, [1 2]));
        end
    end

    % Converting to time
    cmmps = 299792458*1E3*1E-12; % c in mm/ps
    T = zeros([length(DELAYPOS),1]);
    % Guess time zero by peak
    indexoffset = floor(length(DELAYPOS)/10); % using an offset to ignore sharp decay at start
    [~, t0index] = max(SIGNAL(indexoffset:end)); % estimating t0 by peak
    T(:) = (DELAYPOS(:)-DELAYPOS(t0index+indexoffset-1))*4/cmmps; % time vector in ps

    return
end
