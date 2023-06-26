%%% Batch processing of 2D sweep data
%%% initially written from batchfit_temporal_v4

% Initialize
clear; clc; close all;

% Configuration options
testrunyn = 0; % 1 = no-save test run with a few files, 0 = process all
targetfolders = [4 5 7]; % indices of folders to process, [] = all
t0pos = 225.1; % specify t0 position (mm), [] = autofind

writeindyn = 1; % 1 = write individual _proc.dat files, 0 = not
writebatchyn = 0; % 1 = write batch_processed file, 0 = not
writefigsyn = 1; % 1 = write figure files, 0 = not
fileexts = 2; % 1 = .txt, 2 = .tif with Tlist in folder
filenameconv = 3; % 0 = other
% 1 = [pumpWL]-[pumppower]-[signalWL]-[IRpower]-[ND]-[PMTgain]-
% [chopper]-[PMTBW]-[etc],
% 2 = [FOV]_[etc]_[size]_[idler]_[DFG]_[power]_[channel]
% 3 = [IRWN]_[probeWL]_[etc]_[size]_[idler]_[DFG]_[power]_[channel]

% Analysis options
cutlow = 0.15; % lower baseline fit cutoff, 15th %ile
cuthigh = 0.45; % upper baseline fit cutoff, 55th %ile
if testrunyn == 1 % setup test run conditions
    visibility = 'on'; % show figures
    writeindyn = 0; % don't write individual _proc.dat files
    writebatchyn = 0; % don't write batch_processed file
    writefigsyn = 0; % don't write figure files
else
    visibility = 'off'; % hide figures
end

D = pwd; % get current directory
S = dir(fullfile(D,'*')); % search current directory
N = setdiff({S([S.isdir]).name},{'.','..'}); % subfolders of D

if testrunyn == 1 % folder logic for test run
    if length(targetfolders) >= 3
        folders = targetfolders(1:3);
    else
        folders = targetfolders;
    end
else
    if isempty(targetfolders) == 0
        folders = targetfolders;
    else
        folders = 1:numel(N);
    end
end
for ii = folders % ii = subfolder number
    if fileexts == 1  % look for .txt files in subfolders
        T = dir(fullfile(D,N{ii},'*.txt'));
    end
    if fileexts == 2 % look for .tifs
        T = dir(fullfile(D,N{ii},'*.tif'));
        CH1T = dir(fullfile(D,N{ii},'*CH1.tif'));
        CH2T = dir(fullfile(D,N{ii},'*CH2.tif'));
        if length(CH1T) == length(CH2T)
            bothDCAC = 1;
        end
    end
    C = {T(~[T.isdir]).name}; % all data files
    if testrunyn == 1 % limit test run to 1 file per folder
        totalfilesC = 2;
    else
        totalfilesC = numel(C);
    end
    if isempty(T) == 0 % skip subfolders without data files
        for jj = 1:totalfilesC % jj = file number in subfolder ii
            F = fullfile(D,N{ii},C{jj}); % current file
            % Filename parsing
            folderinfo = split(F,"/"); % split path into folders
            if filenameconv == 3
                subfolderinfo = split(N{ii}," "); % split subfolder at " "s
                fileinfo = split(folderinfo(end),"_"); % split filename at "_"s
                IRWNstrfull = char(fileinfo(1)); % extract IRWN - "1598cm-1"
                IRWNstr = IRWNstrfull(1:end-4); % remove "cm-1" suffix - "1598"
                IRWN = sscanf(IRWNstr,'%f'); % convert to double - 1598
                prbstrfull = char(fileinfo(2)); % probe WL - "780 (1)"
                if length(prbstrfull) >= 6 % "780 (10)" --> 790 nm
                    prbsplit = split(prbstrfull," "); % "780" "(1)"
                    if length(prbsplit) > 1
                        prbWL1char = char(prbsplit(1)); % '780'
                        prbWL1 = sscanf(prbWL1char,'%f'); % 780
                        prbWL2char = char(prbsplit(2)); % '(10)' or 'finalcheck'
                        if length(prbWL2char) >= 5 % 'finalcheck'
                            prbWL2 = 0;
                        else % '(10)'
                            prbWL2str = prbWL2char(2:end-1); % '10'
                            prbWL2 = sscanf(prbWL2str,'%f'); % 10
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
            channel = sscanf(F(end-4),'%f'); % identify DC/AC data type

            % Load data
            if fileexts == 1 % load data from .txt
                data = importdata(F); % load data
                x = data(:,1); % save delay position as x
                DC = data(:,2); % save channel 1 as DC
                signal = data(:,3); % save channel 2 as signal
            end
            if fileexts == 2 % load data from Tlist.txt and .tif
                Tlistname = fullfile(D,N{ii},'Tlist.txt'); % current Tlist
                x = importdata(Tlistname); % load Tlist
                data = double(tiffreadVolume(F));
                imagesize = size(data);
                signal = squeeze(mean(data, [1 2]));
                x = x(1:end-1); % last delay position moved back to start
                signal = signal(1:end-1);
            end
            if isempty(t0pos) == 1 % converting to time
                cmmps = 299792458*1E3*1E-12; % c in mm/ps
                t = zeros([length(x),1]);
                indexoffset = floor(length(x)/6); % using an offset to ignore sharp decay at start
                [sigmax, t0index] = max(signal(indexoffset:end)); % estimating t0 by peak
                t(:) = (x(:)-x(t0index+indexoffset-1))*4/cmmps; % time vector in ps
            else
                cmmps = 299792458*1E3*1E-12; % c in mm/ps
                t = zeros([length(x),1]);
                t(:) = (x(:)-t0pos)*4/cmmps; % time vector in ps
            end

            % Setup baseline fitting cutoffs
            ilow = ceil(length(t)*cutlow); % using ceiling to round up
            ihigh = ceil(length(t)-length(t)*cuthigh);
            
            % Fit baseline to biexponential decay with vertical offset
            tbase = [t(2:ilow); t(ihigh:end)]; % trimmed t
            sigbase = [signal(2:ilow); signal(ihigh:end)]; % trimmed signal
            basefn = @(b) b(1)+b(2)*exp(-b(3)*(tbase-b(4)))-sigbase; % error function
            baseopt=optimoptions(@lsqnonlin);
            baseopt.MaxFunctionEvaluations = 1e6; % let the fit run longer
            baseopt.MaxIterations = 1e6; % let the fit run longer
            baseopt.FunctionTolerance = 1e-12; % make the fit more accurate
            baseopt.OptimalityTolerance = 1e-12;
            baseopt.StepTolerance = 1e-12;
            baseopt.Display = 'off'; % silence console output
            basegs = [min(sigbase) 0.1*(max(sigbase)-min(sigbase)) 0.2 0];
            basefit = lsqnonlin(basefn,basegs,[0 0 0 -20],[10 10 10 20],baseopt);
            basecurve = basefit(1)+basefit(2)*exp(-basefit(3)*(t-basefit(4)));
            bleachrate = basefit(3); % ps-1

            if channel == 1 % DC data processing
                % Make figure annotation
                annot = 'Bleach Rate = ' + string(bleachrate) + ' ps^{-1}';
                annotdim = [0.7 0.6 0.2 0.2]; % [posx posy sizex sizey]
                suffix = '_baseline.png';
                figname = string(F(1:end-4)) + suffix;
    
                % Plot exponential baseline fit
                fig1 = figure('visible',visibility); 
                plot(t,signal,'o',tbase,sigbase,'o',t,basecurve,'-');
                xlabel('Time (ps)'); ylabel('DC signal (AU)')
                legend('Data','Baseline data','Baseline')
                annotation('textbox',annotdim,'String',annot);
                if writefigsyn == 1
                    saveas(fig1,figname)
                end
                
                % Write values to an array
                if bothDCAC == 1
                    savecolDC = (jj+1)/2;
                else
                    savecolDC = jj;
                end
                IRWNs(ii,savecolDC,1) = IRWN;
                prbWLs(ii,savecolDC,1) = prbWL;
                DCbleachrates(ii,savecolDC,1) = bleachrate;
            end

            if channel == 2 % AC data processing
                % Make figure annotation
                annot = 'Bleach Rate = ' + string(bleachrate) + ' ps^{-1}';
                annotdim = [0.7 0.6 0.2 0.2]; % [posx posy sizex sizey]
                suffix = '_corrAC.png';
                figname = string(F(1:end-4)) + suffix;

                % Calculate peak height and integrated signal
                corrsig = signal - basecurve;
                sigpeak = max(corrsig);
                integsig = sum(corrsig);
                
                % Plot data
                fig2 = figure('visible',visibility);
                plot(t,signal,'o',tbase,sigbase,'go',t,basecurve,'-');
                xlabel('Time (ps)'); ylabel('AC signal (AU)');
                legend('Data','Baseline data','Baseline');
                annotation('textbox',annotdim,'String',annot);
                if writefigsyn == 1
                    saveas(fig2,figname)
                end

                % Write values to an array
                if bothDCAC == 1
                    savecolAC = jj/2;
                else
                    savecolAC = jj;
                end
                IRWNs(ii,savecolAC,1) = IRWN;
                prbWLs(ii,savecolAC,1) = prbWL;
                ACpeak(ii,savecolAC,1) = sigpeak;
                ACint(ii,savecolAC,1) = integsig;
                ACbleachrates(ii,savecolAC,1) = bleachrate;
                summary(ii,jj,1)
            end
        end
    end
end

IRWNs = IRWNs.';
prbWLs = prbWLs.';
DCbleachrates = DCbleachrates.';
ACbleachrates = ACbleachrates.';
ACpeak = ACpeak.';
ACint = ACint.';

%%
IR = [1598 1608];
probe = [780 781 782 730 780];

signal = [1 1 1;
    2 2 0;
    0.5 0.5 0.5];

signal2 = [[1598 780 1] [1598 781 1] [1598 782 1];
    [1598 730 2] [1598 780 2] [1598 780 0];
    [1608 780 0.5] [1608 781 0.5] [1608 782 0.5]];

% tstep,t/sig/data,ii,jj
signal3 = zeros(61,3,3,3);
for i=1:length(t) % probably a better way to do this
    signal3(i,1,1,1) = t(i);
    signal3(i,2,1,1) = corrsig(i);
end
% make a summary vector that's zeros for length(t)
% then write summary(1) = IRWN
% sumamry(2) = probe, etc.
signal3(3,3,:,:) = [1 1 1;
    2 2 0;
    0.5 0.5 0.5];

%% Post-batch analysis
% Combine into single vectors
IRfreqraw = [IRWNs(:,4).' IRWNs(:,5).' IRWNs(:,7).'];
probeWLraw = [prbWLs(:,4).' prbWLs(:,5).' prbWLs(:,7).'];
DCbleachraw = [DCbleachrates(:,4).' DCbleachrates(:,5).' DCbleachrates(:,7).'];
ACbleachraw = [ACbleachrates(:,4).' ACbleachrates(:,5).' ACbleachrates(:,7).'];
intACraw = [ACint(:,4).' ACint(:,5).' ACint(:,7).'];
peakACraw = [ACpeak(:,4).' ACpeak(:,5).' ACpeak(:,7).'];

% Remove zero entries - (1:133) (165:246)
IRfreq = [IRfreqraw(1:133) IRfreqraw(165:246)];
probeWL = [probeWLraw(1:133) probeWLraw(165:246)];
DCbleach = [DCbleachraw(1:133) DCbleachraw(165:246)];
ACbleach = [ACbleachraw(1:133) ACbleachraw(165:246)];
intAC = [intACraw(1:133) intACraw(165:246)];
peakAC = [peakACraw(1:133) peakACraw(165:246)];

figure; plot(probeWL,ACbleach,'o')
% Looking at the individual sweeps, can see the profile changing even as we
% sweep the probe... should reanalyze so I can plot the temporal profiles
% of everything together (as in, signal as a function of t and probeWL, one
% contour for each IRfreq)




%% Inspect individual files
ii = 3;
jj = 1;

F = fullfile(D,N{ii},C{jj}); % current file
% Load data
if fileexts == 1 % load data from .txt
    data = importdata(F); % load data
    x = data(:,1); % save delay position as x
    DC = data(:,2); % save channel 1 as DC
    signal = data(:,3); % save channel 2 as signal
end
if fileexts == 2 % load data from Tlist.txt and .tif
    Tlistname = fullfile(D,N{ii},'Tlist.txt'); % current Tlist
    x = importdata(Tlistname); % load Tlist
    data = double(tiffreadVolume(F));
    imagesize = size(data);
    signal = squeeze(mean(data, [1 2]));
    x = x(1:end-1); % last delay position moved back to start
    signal = signal(1:end-1);
end
if isempty(t0pos) == 1 % converting to time
    cmmps = 299792458*1E3*1E-12; % c in mm/ps
    t = zeros([length(x),1]);
    indexoffset = floor(length(x)/6); % using an offset to ignore sharp decay at start
    [sigmax, t0index] = max(signal(indexoffset:end)); % estimating t0 by peak
    t(:) = (x(:)-x(t0index+indexoffset-1))*4/cmmps; % time vector in ps
else
    cmmps = 299792458*1E3*1E-12; % c in mm/ps
    t = zeros([length(x),1]);
    t0index = find(x == t0pos);
    t(:) = (x(:)-x(t0index))*4/cmmps; % time vector in ps
end

% Setup baseline fitting cutoffs
ilow = ceil(length(t)*cutlow); % using ceiling to round up
ihigh = ceil(length(t)-length(t)*cuthigh);

% Fit baseline to biexponential decay with vertical offset
tbase = [t(1:ilow); t(ihigh:end)]; % trimmed t
sigbase = [signal(1:ilow); signal(ihigh:end)]; % trimmed signal
basefn = @(b) b(1)+b(2)*exp(-b(3)*(tbase-b(4)))-sigbase; % error function
baseopt=optimoptions(@lsqnonlin);
baseopt.MaxFunctionEvaluations = 1e6; % let the fit run longer
baseopt.MaxIterations = 1e6; % let the fit run longer
baseopt.FunctionTolerance = 1e-12; % make the fit more accurate
baseopt.OptimalityTolerance = 1e-12;
baseopt.StepTolerance = 1e-12;
%baseopt.Display = 'off'; % silence console output
basegs = [min(sigbase) 0.1*(max(sigbase)-min(sigbase)) 0.2 0];
basefit = lsqnonlin(basefn,basegs,[0 0 0 -20],[10 10 10 20],baseopt);
basecurve = basefit(1)+basefit(2)*exp(-basefit(3)*(t-basefit(4)));
bleachrate = basefit(3); % ps-1

% Make figure annotation
annot = 'Bleach Rate = ' + string(bleachrate) + ' ps^{-1}';
suffix = '_baseline.png';
annotdim = [0.6 0.5 0.2 0.1]; % [posx posy sizex sizey]

% Plot exponential baseline fit
fig1 = figure; plot(t,signal,'o',tbase,sigbase,'o',t,basecurve,'-');
xlabel('Time (ps)'); ylabel('DC signal (AU)')
legend('Data','Baseline data','Baseline')
annotation('textbox',annotdim,'String',annot);
figname = string(F(1:end-4)) + suffix;
if writefigsyn == 1
    saveas(fig1,figname)
end
