%%% Batch processing of 2D sweep data
%%% initially written from batchfit_temporal_v4
%%% v2 - added writing to master array (4D data hypercube)

% Initialize
clear; clc; close all;

% Configuration options - may replace with menu, listdlg, or inputdlg
testrunyn = 0; % 1 = no-save test run with a few files, 0 = process all
targetfolders = [4 5 7]; % indices of folders to process, [] = dialog
t0pos = 225.1; % specify t0 position (mm), [] = autofind

pairedDCAC = 1; % 1 = both DC and AC files present and paired, 0 = only AC
writeprocyn = 1; % 1 = write batch processed files, 0 = not
writefigsyn = 1; % 1 = write figure files, 0 = not
fitbaseyn = 1; % 1 = fit baselines and correct, 0 = not
fileexts = 2; % 1 = .txt, 2 = .tif with Tlist in folder
filenameconv = 3; % 0 = other
% 1 = [pumpWL]-[pumppower]-[signalWL]-[IRpower]-[ND]-[PMTgain]-
% [chopper]-[PMTBW]-[etc],
% 2 = [FOV]_[etc]_[size]_[idler]_[DFG]_[power]_[channel]
% 3 = [IRWN]_[probeWL]_[etc]_[size]_[idler]_[DFG]_[power]_[channel]
trimlastT = 1; % 1 = remove last delay position (if sent back to start)
% 0 = retain last delay position (default)

D = pwd; % get current directory
S = dir(fullfile(D,'*')); % search current directory
N = setdiff({S([S.isdir]).name},{'.','..'}); % subfolders of D

% Analysis options
cutlow = 0.15; % lower baseline fit cutoff, 15th %ile
cuthigh = 0.65; % upper baseline fit cutoff, 65th %ile
if testrunyn == 1 % setup test run conditions
    visibility = 'on'; % show figures
    writeprocyn = 0; % don't write batch processed files
    writefigsyn = 0; % don't write figure files
    if length(targetfolders) >= 3 % folder logic for test run
        folders = targetfolders(1:3);
    else
        folders = targetfolders;
    end
else % not a test run
    visibility = 'off'; % hide figures
    if isempty(targetfolders) == 0 % if target folders are specified
        folders = targetfolders;
    else
        folders = 1:numel(N); % all folders
    end
end

% Pre-loop to figure out dimensions of master array
maxTlength = 0; % counter variable
maxfilesinfolder = 0; % counter variable
for ii = folders % pre-loop to figure out dimensions of master array
    if fileexts == 1  % look for .txt files in subfolders
        T = dir(fullfile(D,N{ii},'*.txt'));
    end
    if fileexts == 2 % look for CH2.tifs
        T = dir(fullfile(D,N{ii},'*CH2*.tif'));
    end
    C = {T(~[T.isdir]).name}; % all data files
    if numel(C) > maxfilesinfolder
        maxfilesinfolder = numel(C);
    end
    if isempty(T) == 0 % skips subfolders without data files
        for jj = 1:numel(C) % jj = file number in subfolder ii
            F = fullfile(D,N{ii},C{jj}); % current file
            if fileexts == 1 % load data from .txt
                data = importdata(F); % load data
                x = data(:,1); % save delay position as x
                if trimlastT == 1
                    x = x(1:end-1);
                end
            end
            if fileexts == 2 % load data from Tlist.txt and .tif
                Tlistname = fullfile(D,N{ii},'Tlist.txt'); % current Tlist
                x = importdata(Tlistname); % load Tlist
                if trimlastT == 1
                    x = x(1:end-1);
                end
            end
            currentTlength = length(x);
            if currentTlength > maxTlength
                maxTlength = currentTlength;
            end
        end
    end
end

if isempty(targetfolders) == 1 % dialog to select folders
    [indx,~] = listdlg('PromptString',{'Select folders to process.'},...
    'SelectionMode','multiple','ListString',N);
    targetfolders = indx;
end

if maxTlength < 20 % make sure master array is large enough
    maxTlength = 20;
end

% Make master array to store all data
master = NaN(maxTlength,6,length(folders),maxfilesinfolder);
if pairedDCAC == 1 % add 3 columns for DC raw, corrected, & fit
    master = NaN(maxTlength,9,length(folders),maxfilesinfolder);
end

% Loop through all files and analyze

for ii = folders % ii = subfolder number
    if fileexts == 1  % look for .txt files in subfolders
        T = dir(fullfile(D,N{ii},'*.txt'));
    end
    if fileexts == 2 % look for .tifs
        T = dir(fullfile(D,N{ii},'*CH2*.tif'));
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
                if pairedDCAC == 1
                    DC = data(:,2); % save channel 1 as DC
                    signal = data(:,3); % save channel 2 as signal
                else
                    if width(data) >= 3
                        signal = data(:,3); % save channel 2 as signal
                    else
                        signal = data(:,2);
                    end
                end
                if trimlastT == 1
                    x = x(1:end-1);
                    DC = DC(1:end-1);
                    signal = signal(1:end-1);
                end
            end
            if fileexts == 2 % load data from Tlist.txt and .tif
                Tlistname = fullfile(D,N{ii},'Tlist.txt'); % current Tlist
                x = importdata(Tlistname); % load Tlist
                data = double(tiffreadVolume(F));
                imagesize = size(data);
                signal = squeeze(mean(data, [1 2]));
                sigsds = squeeze(std(data, 0, [1 2]));
                if pairedDCAC == 1
                    F_DC = F;
                    F_DC(end-4) = "1";
                    DCdata = double(tiffreadVolume(F_DC));
                    DC = squeeze(mean(DCdata, [1 2]));
                    DCsds = squeeze(std(DCdata, 0, [1 2]));
                end
                if trimlastT == 1
                    x = x(1:end-1);
                    DC = DC(1:end-1);
                    signal = signal(1:end-1);
                    sigsds = sigsds(1:end-1);
                    DCsds = DCsds(1:end-1);
                end
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
            ihigh = ceil(length(t)*cuthigh);
            
            if fitbaseyn == 1 % fit baseline to exponential w/ vertical offset
                tbase = [t(2:ilow); t(ihigh:end)]; % trimmed t
                sigbase = [signal(2:ilow); signal(ihigh:end)]; % trimmed signal
                sbr = max(signal(ilow:ihigh))/mean(signal(ihigh:end));
                basefn = @(b) b(1)+b(2)*exp(-b(3)*(tbase-b(4)))-sigbase; % error function
                baseopt=optimoptions(@lsqnonlin);
                baseopt.MaxFunctionEvaluations = 1e6; % let the fit run longer
                baseopt.MaxIterations = 1e6; % let the fit run longer
                baseopt.FunctionTolerance = 1e-12; % make the fit more accurate
                baseopt.OptimalityTolerance = 1e-12; % make the fit more accurate
                baseopt.StepTolerance = 1e-12; % make the fit more precise
                baseopt.Display = 'off'; % silence console output
                basegs = [min(sigbase) 0.1*(max(sigbase)-min(sigbase)) 0.01 0];
                basefit = lsqnonlin(basefn,basegs,[],[],baseopt);
                basecurve = basefit(1)+basefit(2)*exp(-basefit(3)*(t-basefit(4)));
                if basefit(2) > 0 && basefit(3) > 0 
                    bleachrate = basefit(3); % ps-1
                else % bleach rate doesn't have physical meaning
                    bleachrate = 0;
                end
                if pairedDCAC == 1 % fit DC baseline
                    DCbase = [DC(2:ilow); DC(ihigh:end)]; % trimmed DC
                    DCsbr = max(DC(ilow:ihigh))/mean(DC(ihigh:end));
                    DCbasefn = @(b) b(1)+b(2)*exp(-b(3)*(tbase-b(4)))-DCbase;
                    DCbgs = [min(DCbase) 0.1*(max(DCbase)-min(DCbase)) 0.01 0];
                    DCbasefit = lsqnonlin(DCbasefn,DCbgs,[],[],baseopt);
                    DCbasecurve = DCbasefit(1)+DCbasefit(2)*exp(-DCbasefit(3)*(t-DCbasefit(4)));
                    if DCbasefit(2) > 0 && DCbasefit(3) > 0
                        DCbleachrate = DCbasefit(3); % ps-1
                    else
                        DCbleachrate = 0;
                        DCsbr = 0;
                    end
                end
            else % no fitting, set basecurve to 0
                basecurve = zeros(height(t),width(t));
                bleachrate = 0;
                sbr = 0;
                if pairedDCAC == 1
                    DCbasecurve = zeros(height(t),width(t));
                    DCbleachrate = 0;
                    DCsbr = 0;
                end
            end

            % Make figure annotation
            annot = 'IR = '+string(IRWN)+' cm-1, probe = '+string(prbWL)+...
                ' nm, bleach = '+string(bleachrate)+' ps^{-1}';
            annotdim = [0.35 0.6 0.18 0.2]; % [posx posy sizex sizey]
            if pairedDCAC == 1
                annot = annot+'(DC = '+string(DCbleachrate)+' ps^{-1})';
            end
            if fileexts == 1
                outname = string(F(1:end-4));
            end
            if fileexts == 2
                if pairedDCAC == 1
                    outname = string(F(1:end-8)) + "_DCAC";
                else
                    outname = string(F(1:end-4));
                end
            end

            % Calculate peak height and integrated signal
            corrsig = signal - basecurve;
            [sigpeak,maxindex] = max(corrsig);
            peaksd = sigsds(maxindex);
            integsig = sum(corrsig);
            integsd = mean(sigsds,'all');
            corrsigbase = [corrsig(2:ilow); corrsig(ihigh:end)];
            noise = range(corrsigbase);
            snr = sigpeak/noise;
            snrsweep = corrsig/noise;
            intsnr = sum(snrsweep);

            if pairedDCAC == 1
                corrDC = DC - DCbasecurve;
                [DCpeak,DCmaxindex] = max(corrDC);
                DCpeaksd = DCsds(DCmaxindex);
                integDC = sum(corrDC);
                DCintegsd = mean(DCsds,'all');
                corrDCbase = [corrDC(2:ilow); corrDC(ihigh:end)];
                DCnoise = range(corrDCbase);
                DCsnr = DCpeak/DCnoise;

                % Plot DC and AC side by side
                fig = figure('visible',visibility); tiledlayout(1,2);
                nexttile([1 1]); 
                hold on; errorbar(t,DC,DCsds,'o');
                plot(tbase,DCbase,'o',t,DCbasecurve,'-'); hold off;
                xlabel('Time (ps)'); ylabel('DC signal (AU)');
                legend('Data','Baseline data','Baseline');
                nexttile([1 1]); 
                hold on; errorbar(t,signal,sigsds,'o');
                plot(tbase,sigbase,'go',t,basecurve,'-'); hold off;
                xlabel('Time (ps)'); ylabel('AC signal (AU)');
                legend('Data','Baseline data','Baseline');
                annotation('textbox',annotdim,'String',annot);
                if writefigsyn == 1
                    saveas(fig,outname+'_proc.png')
                end
            else % plot only signal
                fig = figure('visible',visibility);
                hold on; errorbar(t,signal,sigsds,'o')
                plot(tbase,sigbase,'go',t,basecurve,'-'); hold off
                xlabel('Time (ps)'); ylabel('AC signal (AU)');
                legend('Data','Baseline data','Baseline');
                annotation('textbox',annotdim,'String',annot);
                if writefigsyn == 1
                    saveas(fig,outname+'_proc.png')
                end
            end

            % Write out data to master array
            master(1:length(t),1,ii,jj) = t(:);
            master(1:length(corrsig),2,ii,jj) = corrsig(:);
            master(1,3,ii,jj) = IRWN;
            master(2,3,ii,jj) = prbWL;
            master(3,3,ii,jj) = sigpeak;
            master(4,3,ii,jj) = peaksd;
            master(5,3,ii,jj) = integsig;
            master(6,3,ii,jj) = integsd;
            master(7,3,ii,jj) = snr;
            master(8,3,ii,jj) = sbr;
            master(9,3,ii,jj) = noise;
            master(10,3,ii,jj) = intsnr;
            master(11,3,ii,jj) = bleachrate;
            master(12,3,ii,jj) = length(t);
            master(1:length(signal),4,ii,jj) = signal(:);
            master(1:length(basecurve),5,ii,jj) = basecurve(:);
            master(1:length(basecurve),6,ii,jj) = snrsweep(:);

            if pairedDCAC == 1
                master(1:length(DC),7,ii,jj) = DC(:);
                master(1:length(corrDC),8,ii,jj) = corrDC(:);
                master(1:length(DCbasecurve),9,ii,jj) = DCbasecurve(:);
                master(13,3,ii,jj) = DCpeak;
                master(14,3,ii,jj) = DCpeaksd;
                master(15,3,ii,jj) = integDC;
                master(16,3,ii,jj) = DCintegsd;
                master(17,3,ii,jj) = DCsnr;
                master(18,3,ii,jj) = DCsbr;
                master(19,3,ii,jj) = DCnoise;
                master(20,3,ii,jj) = DCbleachrate; % why maxTlength >= 20
            end
            
            if writeprocyn == 1 % write individual files
                temp = master(:,:,ii,jj);
                writematrix(temp,outname+'_proc.dat')
            end
        end
    end
end

if writeprocyn == 1 % write batch file
    writematrix(master,'batch_flattened.dat') % backup, not very useful
    writematrix(size(master),'master_size.dat')
end
if testrunyn == 0 % close background figures if not a test run
    close all;
end

%% Reload previously processed data
master_size = importdata('master_size.dat');
master = NaN(master_size(1),master_size(2),master_size(3),master_size(4));
D = pwd; % get current directory
S = dir(fullfile(D,'*')); % search current directory
N = setdiff({S([S.isdir]).name},{'.','..'}); % subfolders of D
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.dat'));
    C = {T(~[T.isdir]).name}; % all data files
    if isempty(T) == 0
        for jj = 1:numel(C)
            F = fullfile(D,N{ii},C{jj});
            master(:,:,ii,jj) = importdata(F);
        end
    end
end

%% Post-batch analysis
clc; close all;

% Master array is 4D (x,y,ii,jj); ii,jj are folder,file
% x = rows (individual t points); y = columns
master_size = importdata('master_size.dat'); % load dimensions
% Column indices
tindex = 1;
corrsigindex = 2;
% 3rd column is individual values
rawsigindex = 4;
sigbaseindex = 5;
snrsweepindex = 6;
rawDCindex = 7;
corrDCindex = 8;
DCbaseindex = 9;
% Row indices for individual values
IRindex = 1;
probeindex = 2;
sigpeakindex = 3;
peaksdindex = 4;
integsigindex = 5;
integsdindex = 6;
snrindex = 7;
sbrindex = 8;
noiseindex = 9;
intsnrindex = 10;
bleachindex = 11;
tlistindex = 12;
DCpeakindex = 13;
DCpeaksdindex = 14;
integDCindex = 15;
DCintegsdindex = 16;
DCsnrindex = 17;
DCsbrindex = 18;
DCnoiseindex = 19;
DCbleachindex = 20;

% Photobleaching at 1598cm-1,780 nm - need to specify Tlist length
t780_4_final = squeeze(master(1:master(tlistindex,3,4,1),1,4,1));
sig780_4_final = squeeze(master(1:master(tlistindex,3,4,1),2,4,1));
t780_4_first = squeeze(master(1:master(tlistindex,3,4,2),1,4,2));
sig780_4_first = squeeze(master(1:master(tlistindex,3,4,2),2,4,2));
t780_5_final = squeeze(master(1:master(tlistindex,3,5,1),1,5,51));
sig780_5_final = squeeze(master(1:master(tlistindex,3,5,1),2,5,51));

%bleach780 = figure; hold on; plot(t780_4_first,sig780_4_first);
%plot(t780_4_final,sig780_4_final);
%plot(t780_5_final,sig780_5_final); hold off;
%xlabel('Time (ps)'); ylabel('AC signal (AU)');
%legend('780nm start','780nm after red sweep','780nm after red and blue sweeps')
% suggests bleaching is worse for shorter wavelengths, as expected

% Photobleaching rates as a function of probe WL
%ACbleachprobe = [];
%ACbleach = [];
DCbleachprobe = [];
DCbleach = [];
for i=1:master_size(3) % pull nonzero, non-NaN values and associated probe
    for j=1:master_size(4)
        %if master(bleachindex,3,i,j) ~= 0
        %    ACbleachprobe = [ACbleachprobe master(probeindex,3,i,j)];
        %    ACbleach = [ACbleach master(bleachindex,3,i,j)];
        %end
        if master(DCbleachindex,3,i,j) ~= 0 && ~isnan(master(DCbleachindex,3,i,j)) == 1
            DCbleachprobe = [DCbleachprobe master(probeindex,3,i,j)];
            DCbleach = [DCbleach master(DCbleachindex,3,i,j)];
        end
    end
end

%2channelbleach = figure; tiledlayout(1,2); nexttile([1 1]);
%plot(DCbleachprobe,DCbleach,'o'); nexttile([1 1]);
%plot(ACbleachprobe,ACbleach,'o');

% AC bleaching trends aren't useful (makes sense because that's bleaching
% that oscillated at 10 kHz -- nonsense)

[sortbleach,sortindices] = sort(DCbleachprobe);
DCbleachprobesorted = DCbleachprobe(sortindices);
DCbleachsorted = DCbleach(sortindices);

DCbleachfig = figure('visible','off'); 
plot(DCbleachprobesorted,DCbleachsorted)
xlabel('λ_{probe} (nm)'); ylabel('Photobleaching rate (ps^{-1})')
annotation('line',[0.28 0.28],[0.12 0.92],'LineStyle','--')
annotation('textarrow',[0.5 0.29],[0.5 0.29],'String','747 nm')

% Bleaching is minimal for probe > 747 nm
% -- will sweep 750 and redder from now on

% 1598cm-1 red sweep
redprobe1598 = [];
redpeak1598 = [];
redsd1598 = [];
redint1598 = [];
redintsd1598 = [];
rednoise1598 = [];
redsnr1598 = [];
redsbr1598 = [];
% 1598cm-1 blue sweep
blueprobe1598 = [];
bluepeak1598 = [];
bluesd1598 = [];
blueint1598 = [];
blueintsd1598 = [];
bluenoise1598 = [];
bluesnr1598 = [];
bluesbr1598 = [];
% 1608cm-1 sweep
probe1608 = [];
peak1608 = [];
sd1608 = [];
int1608 = [];
intsd1608 = [];
noise1608 = [];
snr1608 = [];
sbr1608 = [];

for i=4 % 1598cm-1 red sweep
    for j=2:master_size(4) % skip photobleach check at start
        if ~isnan(master(probeindex,3,i,j)) == 1 && master(probeindex,3,i,j) ~= 0
            redprobe1598 = [redprobe1598 master(probeindex,3,i,j)];
            redpeak1598 = [redpeak1598 master(sigpeakindex,3,i,j)];
            redsd1598 = [redsd1598 master(peaksdindex,3,i,j)];
            redint1598 = [redint1598 master(integsigindex,3,i,j)];
            redintsd1598 = [redintsd1598 master(integsdindex,3,i,j)];
            rednoise1598 = [rednoise1598 master(noiseindex,3,i,j)];
            redsnr1598 = [redsnr1598 master(snrindex,3,i,j)];
            redsbr1598 = [redsbr1598 master(sbrindex,3,i,j)];
        end
    end
end
for i=5 % 1598cm-1 blue sweep
    for j=1:master_size(4)
        if ~isnan(master(probeindex,3,i,j)) == 1 && master(probeindex,3,i,j) ~= 0
            blueprobe1598 = [blueprobe1598 master(probeindex,3,i,j)];
            bluepeak1598 = [bluepeak1598 master(sigpeakindex,3,i,j)];
            bluesd1598 = [bluesd1598 master(peaksdindex,3,i,j)];
            blueint1598 = [blueint1598 master(integsigindex,3,i,j)];
            blueintsd1598 = [blueintsd1598 master(integsdindex,3,i,j)];
            bluenoise1598 = [bluenoise1598 master(noiseindex,3,i,j)];
            bluesnr1598 = [bluesnr1598 master(snrindex,3,i,j)];
            bluesbr1598 = [bluesbr1598 master(sbrindex,3,i,j)];
        end
    end
end
for i=7 % 1608cm-1 sweep
    for j=1:master_size(4)-1 % skip photobleach check at end
        if ~isnan(master(probeindex,3,i,j)) == 1 && master(probeindex,3,i,j) ~= 0
            probe1608 = [probe1608 master(probeindex,3,i,j)];
            peak1608 = [peak1608 master(sigpeakindex,3,i,j)];
            sd1608 = [sd1608 master(peaksdindex,3,i,j)];
            int1608 = [int1608 master(integsigindex,3,i,j)];
            intsd1608 = [intsd1608 master(integsdindex,3,i,j)];
            noise1608 = [noise1608 master(noiseindex,3,i,j)];
            snr1608 = [snr1608 master(snrindex,3,i,j)];
            sbr1608 = [sbr1608 master(sbrindex,3,i,j)];
        end
    end
end

% SNR and SBR plots
snrfig = figure;%('visible','off');
hold on; plot(redprobe1598,redsnr1598);
plot(blueprobe1598,bluesnr1598,'b'); plot(probe1608,snr1608); hold off;
xlabel('λ_{probe} (nm)'); ylabel('SNR')
legend('1598 cm^{-1}, sweep 1','1598 cm^{-1}, sweep 2','1608 cm^{-1}')

%sbrfig = figure;%('visible','off');
%hold on; plot(redprobe1598,redsbr1598);
%plot(blueprobe1598,bluesbr1598,'b'); plot(probe1608,sbr1608); hold off;
%xlabel('λ_{probe} (nm)'); ylabel('SBR')
%legend('1598 cm^{-1}, sweep 1','1598 cm^{-1}, sweep 2','1608 cm^{-1}')

% Normalizing to 780 nm for stitching
normfactor = redpeak1598(1)/bluepeak1598(end);
normbluepeak1598 = bluepeak1598*normfactor;
normbluesd1598 = bluesd1598*normfactor;
normbluenoise1598 = bluenoise1598*normfactor;
intnormfactor = redint1598(1)/blueint1598(end);
normblueint1598 = blueint1598*intnormfactor;
normblueintsd1598 = blueintsd1598*intnormfactor;

% Normalized peak stitch
normstitch = figure('visible','off'); 
hold on; errorbar(redprobe1598,redpeak1598,[])
errorbar(blueprobe1598,normbluepeak1598,[],'b');
errorbar(probe1608,peak1608,[]); %xlim([750 860]);
xlabel('λ_{probe} (nm)'); ylabel('Normalized peak signal (AU)')
legend('1598 cm^{-1}, sweep 1','1598 cm^{-1}, sweep 2','1608 cm^{-1}')

normstitcheb = figure('visible','off'); 
hold on; errorbar(redprobe1598,redpeak1598,redsd1598)
errorbar(blueprobe1598,normbluepeak1598,normbluesd1598,'b');
errorbar(probe1608,peak1608,sd1608); %xlim([750 860]);
xlabel('λ_{probe} (nm)'); ylabel('Normalized peak signal (AU)')
legend('1598 cm^{-1}, sweep 1','1598 cm^{-1}, sweep 2','1608 cm^{-1}')

% Normalized integrated stitch (less noisy than peak data)
intstitch = figure('visible','off'); 
hold on; errorbar(redprobe1598,redint1598,[])
errorbar(blueprobe1598,normblueint1598,[],'b');
errorbar(probe1608,int1608,[]); %xlim([750 860]);
xlabel('λ_{probe} (nm)'); ylabel('Normalized integrated signal (AU)')
legend('1598 cm^{-1}, sweep 1','1598 cm^{-1}, sweep 2','1608 cm^{-1}')

intstitcheb = figure('visible','off'); 
hold on; errorbar(redprobe1598,redint1598,redintsd1598)
errorbar(blueprobe1598,normblueint1598,normblueintsd1598,'b');
errorbar(probe1608,int1608,intsd1608); %xlim([750 860]);
xlabel('λ_{probe} (nm)'); ylabel('Normalized integrated signal (AU)')
legend('1598 cm^{-1}, sweep 1','1598 cm^{-1}, sweep 2','1608 cm^{-1}')

% less noisy, but that's just from the difference between a peak height and
% integral (the integral is bigger)

% Making plots with error bar area shading - cool, but maybe not too useful
redint1598above = redint1598+redintsd1598;
redint1598below = redint1598-redintsd1598;
redprobe1598shade = [redprobe1598, fliplr(redprobe1598)];
redint1598shade = [redint1598above, fliplr(redint1598below)];

blueint1598above = normblueint1598+normblueintsd1598;
blueint1598below = normblueint1598-normblueintsd1598;
blueprobe1598shade = [blueprobe1598, fliplr(blueprobe1598)];
blueint1598shade = [blueint1598above, fliplr(blueint1598below)];

int1608above = int1608+intsd1608;
int1608below = int1608-intsd1608;
probe1608shade = [probe1608, fliplr(probe1608)];
int1608shade = [int1608above, fliplr(int1608below)];

shade = figure('visible','off'); 
hold on; fill(redprobe1598shade,redint1598shade,'b',...
    'FaceAlpha',0.5,'FaceColor',[0 0.4470 0.7410],'EdgeAlpha',0);
fill(blueprobe1598shade,blueint1598shade,'b',...
    'FaceAlpha',0.5,'FaceColor','b','EdgeAlpha',0);
fill(probe1608shade,int1608shade,'b',...
    'FaceAlpha',0.5,'FaceColor',[0.9290 0.6940 0.1250],'EdgeAlpha',0);
plot(redprobe1598,redint1598); plot(blueprobe1598,normblueint1598,'b');
plot(probe1608,int1608); hold off;

% also can see different shapes between integrated and peak plots
% --> lineshapes vary

% Contour of corrsig as a function of probeWL and t at 1598 cm-1
T1598 = master(:,tindex,4,1);
probe1598 = [squeeze(master(probeindex,3,4,:)).'...
    squeeze(master(probeindex,3,5,:)).'];

% Rows are t, columns are probeWL
signormfactor = max(master(:,corrsigindex,4,1))/max(master(:,corrsigindex,5,50));
corrsig1598temp = [squeeze(master(:,corrsigindex,4,:))...
    signormfactor*squeeze(master(:,corrsigindex,5,:))];

% Trim NaN values
probe1598 = probe1598(1:133);
clear corrsig1598;
corrsig1598(:,:) = corrsig1598temp(:,1:133);

% Sort data by increasing probe WL
[sortprobe, sortinds] = sort(probe1598);
probe1598sort = probe1598(sortinds);
corrsig1598sort = corrsig1598(:,sortinds);
logsig1598 = log(abs(corrsig1598sort));

% Make meshgrid
[probe,time] = meshgrid(probe1598sort,T1598);

% 3D plots
%timecont = figure; contourf(probe,time,corrsig1598sort);
%xlabel('λ_{probe} (nm)'); ylabel('Time (ps)')
surf1598 = figure('visible','off'); 
surf(probe,time,corrsig1598sort);
xlabel('λ_{probe} (nm)'); ylabel('Time (ps)')

% Log-scale contour - helps to show contrast
cont1598 = figure;%('visible','off'); 
contourf(probe,time,logsig1598); cb=colorbar;
xlabel('λ_{probe} (nm)'); ylabel('Time (ps)');
cb.Label.String = 'log|Corrected signal| (AU)';
cb.Label.Rotation = 270; cb.Label.VerticalAlignment = "bottom";

% Contour for 1608 cm-1 sweep
T1608 = master(1:master(tlistindex,3,7,1),tindex,7,1);
probe1608 = squeeze(master(probeindex,3,7,:)).'; % already sorted
corrsig1608 = squeeze(master(1:master(tlistindex,3,7,1),corrsigindex,7,:));

% Trim photobleach check (at end)
probe1608 = probe1608(1:end-1);
corrsig1608trim = corrsig1608(:,1:end-1);
logsig1608 = log(abs(corrsig1608trim));

[probe2,time2] = meshgrid(probe1608,T1608);
surf1608 = figure('visible','off'); 
surf(probe2,time2,corrsig1608trim);
xlabel('λ_{probe} (nm)'); ylabel('Time (ps)')

cont1608 = figure;%('visible','off'); 
contourf(probe2,time2,logsig1608); cb=colorbar;
xlabel('λ_{probe} (nm)'); ylabel('Time (ps)');
cb.Label.String = 'log|Corrected signal| (AU)';
cb.Label.Rotation = 270; cb.Label.VerticalAlignment = "bottom";

% 2DVE contour plot - integrated signal as a function of IR and probe
IR = [1598; 1608];
probe1608freq(:) = 1e7./probe1608(:);
check = squeeze(master(probeindex,3,4,:)); % remove second
check2 = squeeze(master(probeindex,3,7,:)); % remove last

sweep1598raw = squeeze(master(integsigindex,3,4,:)).';
sweep1608raw = squeeze(master(integsigindex,3,7,:)).';
sweep1598 = sweep1598raw;
sweep1608 = sweep1608raw;
sweep1598(2) = [];
sweep1608(end) = [];

intsig2D = [sweep1598; sweep1608];

[w1,w2] = meshgrid(probe1608,IR);
cont2D = figure;%('visible','off');
contourf(w1,w2,intsig2D); cb=colorbar;
xlabel('λ_{probe} (nm)'); ylabel('ω_{IR} (cm^{-1})');
cb.Label.String = 'Integrated signal (AU)';
cb.Label.Rotation = 270; cb.Label.VerticalAlignment = "bottom";

[w1a,w2a] = meshgrid(probe1608freq,IR);
cont2Da = figure;%('visible','off');
contourf(w1a,w2a,intsig2D); cb=colorbar;
xlabel('ω_{probe} (cm^{-1})'); ylabel('ω_{IR} (cm^{-1})');
cb.Label.String = 'Integrated signal (AU)';
cb.Label.Rotation = 270; cb.Label.VerticalAlignment = "bottom";

% 2DVE SNR contour plot
IR = [1598; 1608];
probe1608freq(:) = 1e7./probe1608(:);

snrsweep1598raw = squeeze(master(intsnrindex,3,4,:)).';
snrsweep1608raw = squeeze(master(intsnrindex,3,7,:)).';
snrsweep1598 = snrsweep1598raw;
snrsweep1608 = snrsweep1608raw;
snrsweep1598(2) = [];
snrsweep1608(end) = [];

snr2D = [snrsweep1598; snrsweep1608];

snrcont2D = figure;%('visible','off');
contourf(w1,w2,snr2D); cb=colorbar;
xlabel('λ_{probe} (nm)'); ylabel('ω_{IR} (cm^{-1})');
cb.Label.String = 'Integrated signal (AU)';
cb.Label.Rotation = 270; cb.Label.VerticalAlignment = "bottom";

snrcont2Da = figure;%('visible','off');
contourf(w1a,w2a,snr2D); cb=colorbar;
xlabel('ω_{probe} (cm^{-1})'); ylabel('ω_{IR} (cm^{-1})');
cb.Label.String = 'Integrated signal (AU)';
cb.Label.Rotation = 270; cb.Label.VerticalAlignment = "bottom";


%% Inspect individual files
D = pwd; % get current directory
S = dir(fullfile(D,'*')); % search current directory
N = setdiff({S([S.isdir]).name},{'.','..'}); % subfolders of D
for ii = folders % loop to scan files
    if fileexts == 1  % look for .txt files in subfolders
        T = dir(fullfile(D,N{ii},'*.txt'));
    end
    if fileexts == 2 % look for CH2.tifs
        T = dir(fullfile(D,N{ii},'*CH2*.tif'));
    end
    C = {T(~[T.isdir]).name}; % all data files
    if numel(C) > maxfilesinfolder
        maxfilesinfolder = numel(C);
    end
    if isempty(T) == 0 % skips subfolders without data files
        for jj = 1:numel(C) % jj = file number in subfolder ii
            F = fullfile(D,N{ii},C{jj}); % current file
            if fileexts == 1 % load data from .txt
                data = importdata(F); % load data
                x = data(:,1); % save delay position as x
                if trimlastT == 1
                    x = x(1:end-1);
                end
            end
            if fileexts == 2 % load data from Tlist.txt and .tif
                Tlistname = fullfile(D,N{ii},'Tlist.txt'); % current Tlist
                x = importdata(Tlistname); % load Tlist
                if trimlastT == 1
                    x = x(1:end-1);
                end
            end
            currentTlength = length(x);
            if currentTlength > maxTlength
                maxTlength = currentTlength;
            end
        end
    end
end

ii = 3;
jj = 1;

F = fullfile(D,N{ii},C{jj}); % current file

%% Planning
% IR = [1598 1608];
% probe = [780 781 782 730 780];
% 
% signal = [1 1 1;
%     2 2 0;
%     0.5 0.5 0.5];
% 
% signal2 = [[1598 780 1] [1598 781 1] [1598 782 1];
%     [1598 730 2] [1598 780 2] [1598 780 0];
%     [1608 780 0.5] [1608 781 0.5] [1608 782 0.5]];
% 
% % tstep,t/sig/data,ii,jj
% signal3 = zeros(61,3,3,3);
% for i=1:length(t) % probably a better way to do this
%     signal3(i,1,1,1) = t(i);
%     signal3(i,2,1,1) = corrsig(i);
% end
% % make a summary vector that's zeros for length(t)
% % then write summary(1) = IRWN
% % sumamry(2) = probe, etc.
% signal3(3,3,:,:) = [1 1 1;
%     2 2 0;
%     0.5 0.5 0.5];
% 
% 
% 
% % Verifying that the batch processing worked
% figure; hold on;
% plot(nonzeros(squeeze(master(2,3,4,:))),nonzeros(squeeze(master(3,3,4,:))));
% plot(nonzeros(squeeze(master(2,3,5,:))),nonzeros(squeeze(master(3,3,5,:))));
% hold off
% 
% figure; hold on;
% plot(nonzeros(squeeze(master(2,3,4,:))),nonzeros(squeeze(master(3,3,4,:))));
% plot(nonzeros(squeeze(master(2,3,7,:))),nonzeros(squeeze(master(3,3,7,:))));
% hold off

%% Extra code

% Un-normalized
%rawstitch = figure; hold on; plot(redprobe1598,redpeak1598,'o')
%plot(blueprobe1598,bluepeak1598,'bo');
%plot(probe1608,peak1608,'o')
%xlabel('λ_{probe} (nm)'); ylabel('Peak signal (AU)')
%legend('1598 cm^{-1}, sweep 1','1598 cm^{-1}, sweep 2','1608 cm^{-1}')

%rawstitcheb = figure; hold on; errorbar(redprobe1598,redpeak1598,redsd1598,'o')
%errorbar(blueprobe1598,bluepeak1598,bluesd1598,'bo');
%errorbar(probe1608,peak1608,sd1608,'o')
%xlabel('λ_{probe} (nm)'); ylabel('Peak signal (AU)')
%legend('1598 cm^{-1}, sweep 1','1598 cm^{-1}, sweep 2','1608 cm^{-1}')

% Error bars from baseline noise - not much different
%normstitchebnoise = figure; hold on; errorbar(redprobe1598,redpeak1598,rednoise1598,'o')
%errorbar(blueprobe1598,normbluepeak1598,normbluenoise1598,'bo');
%errorbar(probe1608,peak1608,noise1608,'o')
%xlabel('λ_{probe} (nm)'); ylabel('Normalized peak signal (AU)')
%legend('1598 cm^{-1}, sweep 1','1598 cm^{-1}, sweep 2','1608 cm^{-1}')

% Raw contour - not useful (info already conveyed in sig vs probe)
%rawsig1598temp = [squeeze(master(:,rawsigindex,4,:))...
%    squeeze(master(:,rawsigindex,5,:))];
%rawsig1598(:,:) = rawsig1598temp(:,1:133);
%rawsig1598sort = rawsig1598(:,sortinds);
%rawtimeplot = figure; surf(probe,time,rawsig1598sort);
%xlabel('λ_{probe} (nm)'); ylabel('Time delay (ps)')

