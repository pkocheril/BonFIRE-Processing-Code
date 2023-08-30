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

% Load data
clear; close all;
F = 'FOV1_BF_50X41_Idler2228.1_DFG5241.3_P0.1368_CH2.tif';
Tlistname = 'Tlist.txt';
outputfile = 'Fit_' + string(F(1:end-4)) + '.png';
tiffname = 'Fit_' + string(F(1:end-4)) + '.tif';

IRWN = 2228;
prbWL = 825;

% Configuration options - check before running!
testrun = 0; % 0 = process all, 1 = no-save test run with a few files
filetype = 3; % 1 = .txt, 2 = .tif (solution), 3 = .tif (image)
t0pos = []; % specify t0 position (mm), [] = autofind
pairedDCAC = 0; % 1 = both DC and AC files present and paired, 0 = only AC
DFGyn = 0; % 1 = IR from DFG, 0 = IR from idler
powernormyn = 0; % 2 = normalize by probe and IR powers, 
% 1 = normalize by IR power, 0 = no normalization
tempmod = 1; % 1 = temporal modulation in place, 0 = no beamsplitters
ND = 1; % ND filter strength in probe path (script will update if possible)
prbpowerset = 250; % probe power in mW (script will update if possible)
IRpowerset = 70; % IR power in mW (script will update if possible)
normIRpower = 50; % IR power on-sample to normalize to (default is 50)
normprobepower = 1.5; % probe power on-sample to normalize to (default 1.5)
trimlastT = 0; % 1 = remove last delay position (if sent back to start)
% 0 = retain last delay position (default)
basefittype = 1; % 2 = exponential baseline fit,
% 1 = linear baseline fit, 0 = no baseline fit
cutlow = 0.10; % lower baseline fit cutoff, (default 10th %ile)
cuthigh = 0.90; % upper baseline fit cutoff, (default 90th %ile)
ltfittype = 1; % 0 = no lifetime fitting, 1 = Gaussian*monoexp,
% 2 = Gaussian*biexp, 3 = beating Gaussian*exp, 4 = beating Gaussian*biexp
setpulsewidth = []; % define pulse width (ps) in fit, [] = float
ltmin = [];  % define minimum lifetime (ps) in fit, [] = default (0.1)
ltmax = []; % define maximum lifetime (ps) in fit, [] = default (100)

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
    if trimlastT == 1
        x = x(1:end-1);
        DC = DC(1:end-1);
        signal = signal(1:end-1);
    end
    % No standard deviation data present - write as zeros
    sigsds = zeros(height(signal),width(signal));
    if pairedDCAC == 1
        DCsds = zeros(height(signal),width(signal));
    end
else % load Tlist and .tif
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
    if trimlastT == 1
        x = x(1:end-1);
        data = data(:,:,1:end-1);
        DCdata = DCdata(:,:,1:end-1);
    end
    if filetype == 2 % average & stdev of XY data in .tif
        signal = squeeze(mean(data, [1 2]));
        sigsds = squeeze(std(data, 0, [1 2]));
        DC = squeeze(mean(DCdata, [1 2]));
        DCsds = squeeze(std(DCdata, 0, [1 2]));
        imagesize = size(signal);
    else  % make arrays to write outputs into
        r2array = zeros([imagesize(1) imagesize(2)]);
        lt1array = r2array; lt2array = r2array;
        ssresidarray = r2array; FWHMarray = r2array;
    end
end

if filetype ~= 3
    imagesize(1) = 1; imagesize(2) = 1;
end

if testrun > 0
    startiii = round(imagesize(1)/2);
    startjjj = round(imagesize(2)/2);
    imagesize(1) = (3-1)+startiii; imagesize(2) = (2-1)+startjjj;
    savedatayn = 0;
else
    savedatayn = 1;
    startiii = 1; startjjj = 1;
end

for iii=startiii:imagesize(1)
    for jjj=startjjj:imagesize(2)
        if filetype == 3 % image data
            signal = squeeze(data(iii,jjj,:)); % pulls sweep data as vector
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
        end
        if powernormyn == 1 % power normalization
            signal = signal*normIRpower/IRpower;
            sigsds = sigsds*normIRpower/IRpower;
            if pairedDCAC == 1
                DC = DC*70/IRpower;
                DCsds = DCsds*70/IRpower;
            end
        else
            if powernormyn == 2 % IR (70 mW) and probe (1.5 mW)
                if tempmod == 1 % correct for temporal modulation
                    prbpower = 0.25*prbpower; % loss from pinhole, 250 mW -> 62.5 mW
                    transmittance = BSdata.data(prbWL,1)/100; % picoEmerald is p-polarized (horizontal)
                    reflectance = BSdata.data(prbWL,4)/100; % CSV in %
                    prbpower = transmittance*reflectance*prbpower; % 62.5 mW -> 15.6 mW
                end
                prbpower = prbpower/(10^(ND)); % ND1: 15.6 mW -> 1.56 mW
                signal = signal*(normIRpower/IRpower)*(normprobepower/prbpower);
                sigsds = sigsds*(normIRpower/IRpower)*(normprobepower/prbpower);
                if pairedDCAC == 1
                    DC = DC*(normIRpower/IRpower)*(normprobepower/prbpower);
                    DCsds = DCsds*(normIRpower/IRpower)*(normprobepower/prbpower);
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
        
        if ltfittype > 0 % lifetime fitting
            deltat = zeros(length(t)-1,1); % check t spacing and length
            for i=1:length(t)-1
                deltat(i) = round(1e3*(t(i+1)-t(i)))/1e3;
            end % * had issues without rounding
            tspacing = unique(deltat); % length 1 <-> evenly spaced t
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
                IRFamp = 50/1e3;
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
            end
        else % no lifetime fitting
            tint = t; sigint = signal; fitcurve = signal;
            lifetime1 = 0; lifetime2 = 0;
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
        sDCbleachrate = string(round(100*DCbleachrate)/100);
        lt1 = string(round(10*lifetime1)/10); 
        lt2 = string(round(10*lifetime2)/10);
        beatstring = string(fitval(9)*1e12/(pi*29979245800))+', '+...
            string(fitval(13)*1e12/(pi*29979245800))+', '+...
            string(fitval(17)*1e12/(pi*29979245800));
        sr2 = string(r2);
        
        % Make figure annotation
        annot = {'IR = '+string(IRWN)+' cm-1','probe = '+...
            string(prbWL)+' nm'};
        if ltfittype > 0
            annotdim = [0.47 0.35 0.2 0.2]; % [posx posy sizex sizey]
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
                lt2+' ps, r^2 = '+sr2};
        end
        if ltfittype == 3
            annot(length(annot)+1) = {'τ_{mono} = '+lt1+...
                ' ps, beats = '+beatstring+' cm^{-1}, r^2 = '+sr2};
        end
        if ltfittype == 4
            annot(length(annot)+1) = {'τ_1 = '+lt1+' ps, τ_2 = '+...
                lt2+' ps, beats = '+beatstring+' cm^{-1}, r^2 = '+sr2};
        end
        if ltfittype > 0 && r2 < 0.7
            annot(end) = {'Poor fit, r^2 = ' + sr2};
        end
        
        % Prepare output filename(s)
        if filetype == 1
            outname = string(F(1:end-4));
        else % fileexts > 1
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
        basesd = std(corrsigbase);
        snr = sigpeak/noise;
        snrsweep = corrsig/noise;
        intsnr = sum(snrsweep);

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
            fig = figure; 
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
        else % plot only signal (AC)
            fig = figure;
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
        end
    end
end

if savedatayn == 1 && filetype == 3 % save data
    Save_tiff(tiffname,lt1array,lt2array,ssresidarray,FWHMarray,r2array)
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
