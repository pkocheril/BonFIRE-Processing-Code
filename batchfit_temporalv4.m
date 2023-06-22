%%% Batch processing and fitting of .txt files and .tif solution files
%%% initially written from fit_temporal_v6
%%% v2 - improved general filename handling
%%% v3 - added lifetime sorting logic from image_fit_temporal_v4
%%% v4 - improved file writeout (batch and individual) -- in progress
%%% -- replaced by sweep_process

% Initialize
clear; clc; close all;

% Configuration options
testrunyn = 1; % 1 = no-save test run with a few files, 0 = process all
writeindyn = 0; % 1 = write individual _proc.dat files, 0 = not
writebatchyn = 0; % 1 = write batch_processed file, 0 = not
writefigsyn = 1; % 1 = write figure files, 0 = not
fileexts = 2; % 1 = .txt, 2 = .tif with Tlist in folder
filenameconv = 2; % 0 = other
% 1 = [pumpWL]-[pumppower]-[signalWL]-[IRpower]-[ND]-[PMTgain]-
% [chopper]-[PMTBW]-[etc],
% 2 = [FOV]_[etc]_[size]_[idler]_[DFG]_[power]_[channel]

% Analysis options
DFGyn = 0; % 1 = IR from DFG, 0 = IR from idler
powernormyn = 0; % 1 = do power normalization, 0 = no normalization
fitAC = 1; % 1 = fit AC data, 0 = fit DC data
cutlow = 0.05; % lower baseline fit cutoff, 5th %ile
cuthigh = 0.05; % upper baseline fit cutoff, 95th %ile
fittype = 2; % 0 = only baseline correction, 1 = Gaussian*monoexp,
% 2 = Gaussian*biexp, 3 = beating Gaussian*exp, 4 = beating Gaussian*biexp
setpulsewidth = 0; % define fixed pulse FWHM in ps for fitting, otherwise
% 0 = float pulse width in fit

% Tables to write data to
Table1 = []; % for writing to individual _proc.dat files
Table2 = []; % for combining all data into one table for further analysis

% Logic to clean up config
if filenameconv == 1 && powernormyn == 1 % when to normalize
    powernormyn = 1; % only power-normalizing for my .txt files
else
    powernormyn = 0;
end
if testrunyn == 1 % setup test run conditions
    visibility = 'on'; % show figures
    writeindyn = 0; % don't write individual _proc.dat files
    writebatchyn = 0; % don't write batch_processed file
    writefigsyn = 0; % don't write figure files
else
    visibility = 'off'; % hide figures
end

% Get files in current directory and subfolders
D = pwd; % get current directory
S = dir(fullfile(D,'*')); % search current directory
N = setdiff({S([S.isdir]).name},{'.','..'}); % subfolders of D

if testrunyn == 1 % limit test run to 3 folders
    totalfilesN = 3;
else
    totalfilesN = numel(N);
end
for ii = 1:totalfilesN % ii = subfolder number
    if fileexts == 1  % look for .txt files in subfolders
        T = dir(fullfile(D,N{ii},'*.txt'));
    end
    if fileexts == 2 % look for Tlists and .tifs
        Tlistlist = dir(fullfile(D,N{ii},'Tlist.txt')); % look for Tlists
        Tlistfiles = {Tlistlist(~[Tlistlist.isdir]).name}; % all Tlists
        if fitAC == 0
            T = dir(fullfile(D,N{ii},'*CH1.tif')); % look for DC .tifs
        else
            T = dir(fullfile(D,N{ii},'*CH2.tif')); % look for AC .tifs
        end
    end
    C = {T(~[T.isdir]).name}; % all data files
    if testrunyn == 1 % limit test run to 1 file per folder
        totalfilesC = 1;
    else
        totalfilesC = numel(C);
    end
    if isempty(T) == 0 % skip subfolders without data files
        for jj = 1:totalfilesC % jj = file number in subfolder ii
            F = fullfile(D,N{ii},C{jj}); % current file
    
            % Load data
            if fileexts == 1 % load data from .txt
                data = importdata(F); % load data
                x = data(:,1); % save delay position as x
                DC = data(:,2); % save channel 1 as DC
                AC = data(:,3); % save channel 2 as AC
            end
            if fileexts == 2 % load data from Tlist.txt and .tif
                Tlistname = fullfile(D,N{ii},'Tlist.txt'); % current Tlist
                x = importdata(Tlistname); % load Tlist
                data = double(tiffreadVolume(F));
                imagesize = size(data);
                if fitAC == 0
                    DC = squeeze(mean(data, [1 2]));
                else
                    AC = squeeze(mean(data, [1 2]));
                end
            end
            if fitAC == 0 % assign data to 'signal'
                signal = DC;
            else
                signal = AC;
            end
            sigamp = max(signal)-min(signal); % used for Gaussian amplitude
    
            % Converting to time
            cmmps = 299792458*1E3*1E-12; % c in mm/ps
            t = zeros([length(x),1]);
            [sigmax, t0index] = max(signal); % estimating t0 by peak
            t(:) = (x(:)-x(t0index))*4/cmmps; % time vector in ps

            % Interpolate for convolution
            samplingRateIncrease = 2;
            tint = linspace(t(1),t(end),length(t)*samplingRateIncrease-1).';
            sigint = spline(t, signal, tint);
            tc = tint(1:2:end); % pulls odd values of tint for convolution

            % Parse filename
            folderinfo = split(F,"/"); % split path into folders
            if filenameconv == 1
                fileinfo = split(folderinfo(end),"-"); % split filename at "-"s
                subfolderinfo = split(N{ii}," "); % split subfolder by spaces
                pumpWLstrfull = char(fileinfo(1)); % extract pump WL - "760m"
                pumpWLstr = pumpWLstrfull(1:end-1); % remove "m" suffix - "760"
                pumpWL = sscanf(pumpWLstr,'%f'); % convert to double
                pumppowerstrfull = char(fileinfo(2)); % extract pump power
                pumppowerstr = pumppowerstrfull(1:end-2); % remove "mW"
                pumppower = sscanf(pumppowerstr,'%f');
                signalWLstrfull = char(fileinfo(3)); % extract signal WL
                signalWLstr = signalWLstrfull(1:end-2); % remove "nm"
                signalWL = sscanf(signalWLstr,'%f');
                IRpowerstrfull = char(fileinfo(4)); % extract IR power
                IRpowerstr = IRpowerstrfull(1:end-2); % remove "mW"
                IRpower = sscanf(IRpowerstr,'%f');
                DFGWL = 1/((2/signalWL) - (1/1031.2)); % calculate DFG WL (nm)
                idlerWL = 1/((1/1031.2)-(1/signalWL));
                if DFGyn == 1
                    IRWN = 1e7/DFGWL;
                else
                    IRWN = 1e7/idlerWL;
                end
    
                % Create bookkeeping string array
                if length(tint) >= 10+length(fileinfo)+length(subfolderinfo)
                    bookkeeping = strings(length(tint),1);
                    bookkeeping(1) = 'subfolder = ' + string(N{ii});
                    bookkeeping(2) = 'file = ' + string(folderinfo(end));
                    bookkeeping(3) = 'pumpWL = ' + string(pumpWL) + ' nm';
                    bookkeeping(4) = 'pumppower = ' + string(pumppower) + ' mW';
                    bookkeeping(5) = 'signalWL = ' + string(signalWL) + ' nm';
                    bookkeeping(6) = 'IRpower = ' + string(IRpower) + ' mW';
                    bookkeeping(7) = 'DFGWN = ' + string(IRWN) + ' cm-1';
                    fileinfostart = 5;
                    for i = fileinfostart:length(fileinfo)
                        bookkeeping(i+9) = string(fileinfo(i));
                    end
                    bkl = length(fileinfo);
                    for i = 1:length(subfolderinfo)
                        bookkeeping(i+9+bkl) = string(subfolderinfo(i));
                    end
                end
                if length(tint) >= 10 % write numeric values to a vector
                    values = zeros(length(tint),1);
                    values(3) = pumpWL;
                    values(4) = pumppower;
                    values(5) = signalWL;
                    values(6) = IRpower;
                    values(7) = IRWN;
                end
            else
                bookkeeping = strings(length(tint),1);
                values = [];
            end

            if filenameconv == 2
                fileinfo = split(folderinfo(end),"_"); % split at "_"s
                channel = char(fileinfo(end)); % 'CH2'
                IRpowermeterstrfull = char(fileinfo(end-1)); % 'P0.4409'
                IRpowermeterstr = IRpowermeterstrfull(2:end); % '0.4409'
                IRpowermeter = sscanf(IRpowermeterstr,'%f'); % 0.4409
                DFGWLstrfull = char(fileinfo(end-2)); % 'DFG3797.9'
                DFGWLstr = DFGWLstrfull(4:end); % '3797.9'
                DFGWL = sscanf(DFGWLstr,'%f'); % 3797.9
                idlerWNstrfull = char(fileinfo(end-3)); % 'Idler2949.8'
                idlerWNstr = idlerWNstrfull(6:end); % '2949.8'
                idlerWN = sscanf(idlerWNstr,'%f'); % 2949.8
                imgdim = char(fileinfo(end-4)); % '5X5'
                pumpWL = 0; pumppower = 0; % <-- unknown from filename
                idlerWL = 1E7/idlerWN;
                signalWL = 1/((1/1031.2)-(1/idlerWL));
                IRpower = IRpowermeter*300; % assuming 300 mW scale
                if DFGyn == 1
                    IRWN = 1e7/DFGWL;
                else
                    IRWN = idlerWN;
                end
                if length(tint) >= 15
                    bookkeeping(1) = 'subfolder = ' + string(N{ii});
                    bookkeeping(2) = 'file = ' + string(folderinfo(end));
                    bookkeeping(3) = 'pumpWL = ' + string(pumpWL) + ' nm';
                    bookkeeping(4) = 'pumppower = ' + string(pumppower) + ' mW';
                    bookkeeping(5) = 'signalWL = ' + string(signalWL) + ' nm';
                    bookkeeping(6) = 'IRpower = ' + string(IRpower) + ' mW';
                    bookkeeping(7) = 'IRWN = ' + string(IRWN) + ' cm-1';
                end
                if length(tint) >= 10 % write numeric values to a vector
                    values = zeros(length(tint),1);
                    values(3) = pumpWL;
                    values(4) = pumppower;
                    values(5) = signalWL;
                    values(6) = IRpower;
                    values(7) = IRWN;
                end
            end
            
            if powernormyn == 1 % normalize to 800 mW probe and 70 mW IR
                pumpnorm = 800/pumppower; % normalization factor for pump
                IRnorm = 70/IRpower; % normalization factor for IR
                signal = signal*pumpnorm*IRnorm;
            end

            % Setup baseline fitting cutoffs - fit first 10% and last 5%
            tlow = t(ceil(length(t)*cutlow)); % using ceiling to round up
            thigh = t(end-ceil(length(t)*cuthigh));
            
            % Setup excluded points for baseline fitting
            excludeDupes = zeros([size(t),1]);
            for i = 1:length(t)
                if t(i) > tlow && t(i) < thigh
                    excludeDupes(i) = i;
                end
            end
            
            % Isolate data positions to exclude (unique, nonzero indices)
            exclude1 = nonzeros(unique(excludeDupes));
            
            % Fit baseline
            base = fit(t,signal,'poly1','Exclude',exclude1); % fit baseline
            basecoef = coeffvalues(base); % save coeff values
            basefitval = polyval(basecoef,t);
            corrsig = signal - basefitval;

            % Calculate peak height and integrated signal
            sigpeak = max(corrsig);
            integsig = sum(corrsig);
            
            % Define error function as fit (parameter vector r) minus signal
            fun = @(r) conv(sigamp*exp(-r(2)*(tc-r(1)).^2), ...
                r(3)*exp(-r(4)*(tc-r(1)))+r(5)*exp(-r(6)*(tc-r(1)))+ ...
                r(7)*exp(-r(8)*(tc-r(1))).*cos(r(9)*tc-r(10))+ ...
                r(11)*exp(-r(12)*(tc-r(1))).*cos(r(13)*tc-r(14))+ ...
                r(15)*exp(-r(16)*(tc-r(1))).*cos(r(17)*tc-r(18)))+ ...
                r(19)*tint+r(20)-sigint;
            
            % Initial guesses for r - Gaussian terms
            cent = 0; % r(1), ps -- center
            gdec = 4*log(2)/(4^2); % r(2), ps^-2 -- Gaussian decay
            
            % Initial guesses for r - exponential terms
            eamp1 = 1e-5; % r(3) -- amplitude 1
            edec1 = 1/2; % r(4), ps^-1 -- 1/lifetime 1
            eamp2 = 1e-3; % r(5)  -- amplitude 2
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
            
            % Initial guesses for r - baseline terms
            slope = basecoef(1); % r(19), slope from baseline fit
            intercept = basecoef(2); % r(20), intercept from baseline fit
            
            % Fit setup
            r0 = [cent,gdec,eamp1,edec1,eamp2,edec2,... % Gaussian and exp decays
                bamp1,bdec1,fr1,phs1,... % first beating term
                bamp2,bdec2,fr2,phs2... % second beating term
                bamp3,bdec3,fr3,phs3,... % third beating term
                slope,intercept]; % linear baseline
            lb = [-10,0,0,0.01,0,0.01,... % lower bounds
                0,0,0,-180,...
                0,0,0,-180,...
                0,0,0,-180,...
                slope,intercept];
            ub = [20,9,9,9,9,9,... % upper bounds
                9,9,9,180,...
                9,9,9,180,...
                9,9,9,180,...
                slope,intercept];
            
            if fittype == 0 % baseline fit only
                lb(1:18) = 0; ub(1:18) = 0;
            end
            if fittype == 1 % Gaussian*exp
                lb(5:18) = 0; ub(5:18) = 0;
            end
            if fittype == 2 % Gaussian*biexp
                lb(7:18) = 0; ub(7:18) = 0;
            end
            if fittype == 3 % Gaussian*exp with beating
                lb(5:6) = 0; ub(5:6) = 0;
            end
            if setpulsewidth ~= 0
                pulsedecay = log(2)/setpulsewidth^2;
                lb(2) = pulsedecay; ub(2) = pulsedecay;
            end
    
            % Run lsqnonlin - minimizes error function by adjusting r
            options=optimoptions(@lsqnonlin);
            options.MaxFunctionEvaluations = 1e6; % let the fit run longer
            options.MaxIterations = 1e6; % let the fit run longer
            options.FunctionTolerance = 1e-12; % make the fit more accurate
            options.Display = 'off'; % silence console output
            fitval = lsqnonlin(fun,r0,lb,ub,options); % fit coeffs
    
            % Generate curve for fitted parameters
            fitcurve = conv(sigamp*exp(-fitval(2)*(tc-fitval(1)).^2), ...
                fitval(3)*exp(-fitval(4)*(tc-fitval(1)))+ ...
                fitval(5)*exp(-fitval(6)*(tc-fitval(1)))+ ...
                fitval(7)*exp(-fitval(8)*(tc-fitval(1))).* ...
                cos(fitval(9)*tc-fitval(10))+ ...
                fitval(11)*exp(-fitval(12)*(tc-fitval(1))).* ...
                cos(fitval(13)*tc-fitval(14))+ ...
                fitval(15)*exp(-fitval(16)*(tc-fitval(1))).* ...
                cos(fitval(17)*tc-fitval(18)))+...
                fitval(19)*tint+fitval(20);
            
            % Calculate residuals, ssresid, and IRF FWHM
            resid = sigint-fitcurve;
            ssresid = sum(resid.^2); % sum of squares of residuals
            tss = sum((sigint-mean(sigint)).^2); % total sum of squares
            r2 = 1-(ssresid/tss); % coefficient of determination (R^2)
            FWHM = sqrt(log(2)/fitval(2)); % pulse width, ps
            
            % Lifetime calculation - want lifetime1 < lifetime2
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
            else        % means 1/fitval(4) is longer
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
            end
            
            if isempty(values) == 0
                bookkeeping(8) = 'τ1 = ' + string(lifetime1) + ' ps';
                bookkeeping(9) = 'τ2 = ' + string(lifetime2) + ' ps';
                bookkeeping(10) = 'R^2 = ' + string(r2);
                bookkeeping(11) = 'FWHM = ' + string(FWHM) + ' ps';
                bookkeeping(12) = 'Sig peak = ' + string(sigpeak);
                bookkeeping(13) = 'Integ sig = ' + string(integsig);
                values(8) = lifetime1;
                values(9) = lifetime2;
                values(10) = r2;
                values(11) = FWHM;
                values(12) = sigpeak;
                values(13) = integsig;
            else
                bookkeeping(1) = 'subfolder = ' + string(N{ii});
                bookkeeping(2) = 'file = ' + string(folderinfo(end));
                bookkeeping(3) = 'τ1 = ' + string(lifetime1) + ' ps';
                bookkeeping(4) = 'τ2 = ' + string(lifetime2) + ' ps';
                bookkeeping(5) = 'R^2 = ' + string(r2);
                bookkeeping(6) = 'FWHM = ' + string(FWHM) + ' ps';
                bookkeeping(7) = 'Sig peak = ' + string(sigpeak);
                bookkeeping(8) = 'Integ sig = ' + string(integsig);
                values = zeros(length(tint),1);
                values(3) = lifetime1;
                values(4) = lifetime2;
                values(5) = r2;
                values(6) = FWHM;
                values(7) = sigpeak;
                values(8) = integsig;
            end

            % Make figure annotation
            if fittype == 1 % single exponential decay
                annot = 'τ = ' + string(lifetime1) + ' ps';
                suffix = '_fit1e.png';
                annotdim = [0.15 0.65 0.15 0.05]; % [posx posy sizex sizey]
            end
            if fittype == 2 % biexponential decay
                annot = 'τ_{1} = ' + string(lifetime1) + ' ps, τ_{2} = '...
                    + string(lifetime2) + ' ps';
                suffix = '_fit2e.png';
                annotdim = [0.15 0.65 0.15 0.1];
            end

            % Plot fits
            fig1 = figure('visible',visibility); 
            tiledlayout(2,3); nexttile([1 1]); plot(base,t,signal,exclude1);
            xlabel('Time (ps)'); ylabel('Signal (AU)'); 
            nexttile([1 2]); plot(tint, resid)
            legend('Residuals'); xlabel('Time (ps)'); ylabel('Residuals (AU)');
            nexttile([1 3]); plot(t,signal,'o',tint,fitcurve,'-');
            xlabel('Time (ps)'); ylabel('Signal (AU)');
            legend('Data','Fit');
            annotation('textbox',annotdim,'String',annot);
            figname = string(F(1:end-4)) + suffix;
            if writefigsyn == 1
                saveas(fig1,figname)
            end
    
            % Save individual data (raw and fit) to a table
            Table1 = table(bookkeeping,values,tint,sigint,fitcurve,resid);
            procfile = string(F) + '_proc.dat';
            if writeindyn == 1
                writetable(Table1,procfile);
            end
    
            % Setup batch table labels (must be unique)
            filename = char(folderinfo(end));
            if length(filename) >= 41
                tablename = filename(end-40:end-4); % need char to be able to trim by position
            else
                tablename = filename(1:end-4);
            end
            bookname = "Info_" + string(ii) + "_" + string(jj) + "_" + ...
                string(tablename); % need quotes to make a string
            valuename = "Value_" + string(ii) + "_" + string(jj) + "_" + ...
                string(tablename);
            tname = "t_" + string(ii) + "_" + string(jj) + "_" + ...
                string(tablename); % need strings for table variable names
            signame = "sig_" + string(ii) + "_" + string(jj) + "_" + ...
                string(tablename);
            fitname = "fit_" + string(ii) + "_" + string(jj) + "_" + ...
                string(tablename);
            oldnames = ["bookkeeping","values","tint","sigint","fitcurve"]; 
            newnames = [bookname, valuename, tname, signame, fitname];
    
            if writebatchyn == 1 % combine processed data into one table
                TF = isempty(Table2);
                if TF == 1 % for the first dataset
                    Table2 = strings(length(tint),5*numdatafiles);
                    Table2 = table(Table2,bookkeeping,values,tint,sigint,fitcurve);
                    Table2(:,1) = []; % remove empty first column from table
                    % Rename variables to be unique by filename
                    Table2 = renamevars(Table2,oldnames,newnames);
                else % for successive datasets
                    %Table2 = table(Table2,bookkeeping,x,corrnormDC,corrnormAC);
                    datatable = table(bookkeeping,values,tint,sigint,fitcurve);
                    Table2 = [Table2 datatable]; % append new data to previous
                    Table2 = renamevars(Table2,oldnames,newnames);
                end
            end
        end
    end
end

if writebatchyn == 1 % write out all data to one txt file
    writetable(Table2,'batch_fitv4.txt');
end

%% Load batch-processed data

%%% Figure out for if I need to clear the workspace and don't want to wait
%%% for the analysis to run again
%Table2 = readtable('batch_processed.txt','VariableNamingRule','preserve');
% Close but not quite right


%% Case-by-case analysis
% Create infotable to figure out what files went where

[Table2size(1), Table2size(2)] = size(Table2);
infotable = strings([Table2size(1), Table2size(2)/5]);
for i = 1:width(infotable)
    j = 5*i-4;
    infotable(:,i) = Table2(:,j);
end

% Create array with only numeric data from Table2
dataarray = zeros(height(Table2),width(Table2));
dataarray(:,:) = table2array(Table2(:,:));

%x = dataarray(:,3); % in case I load in Table2 from batch_processed data

%%% File annotation
%	1	;	
%	6	;	
%	11	;	
%	16	;	
%	21	;	
%	26	;	
%	31	;	
%	36	;	
%	41	;	
%	46	;	
%	51	;	
%	56	;	
%	61	;	
%	66	;	
%	71	;	
%	76	;	
%	81	;	
%	86	;	
%	91	;	
%	96	;	
%	101	;	
%	106	;	
%	111	;	
%	116	;	
%	121	;	
%	126	;	
%	131	;	
%	136	;	
%	141	;	
%	146	;	
%	151	;	

% Example parsing
%pmtgain = [1 2 4 6 8 10 10 20 40 60 80 100];
%pmtgainlabel(:) = string("Gain " + pmtgain(:));
%pmtarrayAC(:,1) = dataarray(:,190); % PMT1
%pmtarrayAC(:,2) = dataarray(:,210); % PMT2

% Example finding time 0
%cmmps = 299792458*1E3*1E-12; % c in mm/ps
%[pmtarray1peak, t0index] = max(pmtarray(:,1));
%t(:) = (x(:)-x(t0index))*4/cmmps; % time vector in ps

% Example plot
%figure;
%plot(x,pmtarrayAC(:,1),x,pmtarrayAC(:,2),'-o') % the "-" means you get the trace and markers together
%xlabel('Position (mm)')
%ylabel('Corrected AC signal (AU)')
%legend(pmtgainlabel)

% Example peak vs integration data
%pmtmax = max(pmtarrayAC);
%pmtint = sum(pmtarrayAC);

%figure;
%plot(pmtgain,pmtint,'-o') % the "-" means you get the trace and markers together
%xlabel('PMT Gain')
%ylabel('Integrated AC signal (AU)')

% Example 3D-plot
%[gain,pos] = meshgrid(pmtgain,x);
%figure;
%surf(gain,pos,pmtarrayAC)
%xlabel('PMT Gain')
%ylabel('Position (mm)')
%zlabel('Corrected AC signal (AU)')

% Example error bar plot
%figure;
%errorbar(nbconc,meanints,stdevints)
%xlabel('[NB] (mM)')
%ylabel('Mean AC integrated signal (AU)')

%% Processing images
% Example loading measurements from Fiji
%unsortimage = readtable('./Images/BonFIRE/BFimages.csv','VariableNamingRule','preserve');
%imagetable = sortrows(unsortimage,2); % sort by second column (should be filename)

% Isolate RawIntDen (sum of pixel values) in a vector
%imagedata(:,1) = imagetable(:,8); % still a table at this point
%rawintden(:,:) = table2array(imagedata(:,:));

% Separate rawintden into individual columns in an array
%BFimages = zeros(3,9);
%BFimages(:,1) = rawintden(1:3); % 0.5 mM
%BFimages(:,9) = rawintden(4:6); % 10 mM
%BFimages(:,2) = rawintden(7:9); % 1 mM
%BFimages(:,3) = rawintden(10:12); % 3 mM
%BFimages(:,4) = rawintden(13:15); % 5 mM
%BFimages(:,5) = rawintden(16:18); % 6 mM
%BFimages(:,6) = rawintden(19:21); % 7 mM
%BFimages(:,7) = rawintden(22:24); % 8 mM
%BFimages(:,8) = rawintden(25:27); % 9 mM


%% Appendix
