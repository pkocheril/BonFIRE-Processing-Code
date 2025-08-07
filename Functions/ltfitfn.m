% Lifetime fitting
function [TINT,SIGINT,FITVECTOR,TFIT,FITCURVE,FITVAL,LIFETIME1,LIFETIME2,LT1LT2,FWHM,R2,RESID,SRR,FFTX,FFTY,CURRLTFITTYPE] = ...
    ltfitfn(T,SIGNAL,CURRLTFITTYPE,PEAKSIGN,SIGPEAK,SNR,CURRPULSEWIDTH,LTMIN,LTMAX,BASEFIT,CURRBASEFITTYPE,FLOATBASE,IRWNum,PROBE,TROUBLESHOOT)
%%% This function is performs lifetime fitting (convolution of Gaussian
% with vibrational decay) on an input signal vector for a given time
% vector. It can guess what type of fit to use based on frequency.

    if isempty(CURRLTFITTYPE) % auto-choose lifetime fitting
        if IRWNum > 1702 && IRWNum < 1800 % carbonyl
            CURRLTFITTYPE = 1; % monoexponential
        else % not carbonyl
            if IRWNum < 2400 && IRWNum > 2000 % triple bond
                if PROBE < 520 && PROBE > 390 % 515.6 nm or 400 nm
                    CURRLTFITTYPE = 2; % biexponential (Coumarin 337 nitrile or small cyanocoumarin)
                else
                    CURRLTFITTYPE = 1; % monoexponential (normal probe)
                end
            else % double bond or CH
                if (IRWNum == 1300 || IRWNum == 1500 || IRWNum == 1590) && PROBE < 770 && PROBE > 690
                    CURRLTFITTYPE = 4; % stretched/compressed biexponential
                else
                    CURRLTFITTYPE = 2; % biexponential
                end
            end
        end
    end
    % Lifetime fitting
    if CURRLTFITTYPE > 0
        % Prepare for lifetime fitting
        deltat = zeros(length(T)-1,1); % calculate t spacing
        for i=1:length(T)-1
            deltat(i) = round(1e2*(T(i+1)-T(i)))/1e2; % round to nearest 0.01 ps (had issues without rounding)
        end
        tspacing = unique(deltat); % length 1 <-> evenly spaced t
        % Check t spacing (even) and length (odd)
        if isscalar(tspacing) && mod(length(T),2) == 1
            TINT = T; SIGINT = SIGNAL; % data are good for fitting
        else % interpolate for convolution
            TINT = linspace(T(1),T(end),length(T)*2-1).';
            SIGINT = spline(T, SIGNAL, TINT);
            if TROUBLESHOOT == 1
                figure; plot(TINT,SIGINT,'o',T,SIGNAL,'o'); title('Interpolation')
            end
        end
        if max(TINT) < 100 % auto-pad signal if tint doesn't go past 100 ps
            padsignal = 1;
        else
            padsignal = 0;
        end
        % Prepare for signal padding
        tintspace = TINT(2)-TINT(1);
        originallength = length(TINT); % should be an odd number
        if mod(originallength-1,4) == 0
            padlength = round(originallength/2)+1;
        else
            padlength = round(originallength/2);
        end
        tlong = linspace(TINT(1),TINT(end)+tintspace*padlength,length(TINT)+padlength).';
        baselong = BASEFIT(1)+BASEFIT(2)*exp(-BASEFIT(3)*(tlong-BASEFIT(4)))+...
            BASEFIT(5)*tlong;
        if CURRBASEFITTYPE > 0
            siglong = baselong;
            siglong(1:length(SIGINT)) = SIGINT;
        else
            if PEAKSIGN == 1
                siglong(length(SIGINT)+1:length(tlong)) = min(SIGINT(end-5:end));
            else
                siglong(length(SIGINT)+1:length(tlong)) = max(SIGINT(end-5:end));
            end
        end
        % Signal padding for short-tailed data
        if padsignal == 1
            TINT = tlong;
            SIGINT = siglong;
            if TROUBLESHOOT == 1
                figure; plot(TINT,SIGINT,'-o',T,SIGNAL,'o'); title('Padding');
            end
        end
        % Prepare convolution input (same spacing, half length, symmetrically around zero)
        tc = linspace(0.5*min(TINT),0.5*max(TINT),0.5*(length(TINT)+1)).';
        fun = @(r) r(1)+r(2)*exp(-r(3)*(TINT-r(4)))+r(5)*TINT+ ... % baseline
            conv(exp(-((tc-r(6))./(r(7)/(2*sqrt(log(2))))).^2), ...
            heaviside(tc).*(r(8)*exp(-(tc)./r(9))+r(10)*exp(-((tc)./(r(11))).^r(12)))) - SIGINT;
        % Initial guesses
        r0 = [BASEFIT,0,2.6,... % basecoefs, IRF center (ps), IRF width (ps)
            4.5,1,1.5,10,1]; % amp 1, τ1 (ps), amp 2, τ2 (ps), β (stretch)
        % Float baseline coeffs
        newlbbase = BASEFIT; newubbase = BASEFIT;
        for i=1:length(BASEFIT)
            if abs(BASEFIT(i)) < 1e-6 % set very small values to zero
                newlbbase(i) = 0;
                newubbase(i) = 0;
            else
                newlbbase(i) = min((1-FLOATBASE)*BASEFIT(i),(1+FLOATBASE)*BASEFIT(i));
                newubbase(i) = max((1-FLOATBASE)*BASEFIT(i),(1+FLOATBASE)*BASEFIT(i));
            end
        end
        % Lower and upper bounds
        lb = [newlbbase,-10,1*sqrt(2),... % basecoefs, IRF center (ps), IRF min (ps)
            0,0.1,... % amp 1, τ1min (ps)
            0,0.1,... % amp 2, τ2min (ps)
            0.1]; % β min
        ub = [newubbase,20,3*sqrt(2),... % basecoefs, IRF center (ps), IRF max (ps)
            10*SIGPEAK,100,... % amp 1, τ1max (ps)
            10*SIGPEAK,100,... % amp 2, τ2max (ps)
            1.8]; % β max
        % Implement fit type options
        if CURRLTFITTYPE == 1 % Gauss*monoexp
            lb(10) = 0; ub(10) = 0; % first exp only
        end
        if CURRLTFITTYPE == 2 % Gauss*biexp
            lb(12) = 1; ub(12) = 1; % no stretching
        end
        if CURRLTFITTYPE == 3 % Gauss*stretchexp
            lb(8) = 0; ub(8) = 0; % second exp only
        end
        % Parameter constraints from config
        if ~isempty(CURRPULSEWIDTH) % pulse width was specified
            lb(7) = CURRPULSEWIDTH*sqrt(2); ub(7) = CURRPULSEWIDTH*sqrt(2);
        end
        if ~isempty(LTMIN) % minimum lifetime was specified
            lb(9) = LTMIN; lb(11) = LTMIN;
        end
        if ~isempty(LTMAX) % maximum lifetime was specified
            ub(9) = LTMAX; ub(11) = LTMAX;
        end
        % Run lsqnonlin - minimizes error function by adjusting r
        options=optimoptions(@lsqnonlin);
        if SNR > 10 % let run longer for high-SNR
            options.MaxFunctionEvaluations = 1e6; % let the fit run longer
            options.MaxIterations = 1e6; % let the fit run longer
            options.FunctionTolerance = abs(SIGPEAK)*1e-10; % make the fit more accurate
            options.OptimalityTolerance = abs(SIGPEAK)*1e-10; % make the fit more accurate
            options.StepTolerance = abs(SIGPEAK)*1e-10; % make the fit more precise
        else
            options.MaxFunctionEvaluations = 1e5; % let the fit run longer
            options.MaxIterations = 1e5; % let the fit run longer
            options.FunctionTolerance = abs(SIGPEAK)*1e-8; % make the fit more accurate
            options.OptimalityTolerance = abs(SIGPEAK)*1e-8; % make the fit more accurate
            options.StepTolerance = abs(SIGPEAK)*1e-8; % make the fit more precise
        end
        if TROUBLESHOOT == 0
            options.Display = 'off'; % silence console output
        end
        FITVAL = lsqnonlin(fun,r0,lb,ub,options); % fit coeffs
        FITVECTOR = FITVAL(1)+FITVAL(2)*exp(-FITVAL(3)*(TINT-FITVAL(4)))+FITVAL(5)*TINT+ ... % baseline
            conv(exp(-((tc-FITVAL(6))./(FITVAL(7)/(2*sqrt(log(2))))).^2), ...
            heaviside(tc).*(FITVAL(8)*exp(-(tc)./FITVAL(9))+FITVAL(10)*exp(-((tc)./(FITVAL(11))).^FITVAL(12))));
        if TROUBLESHOOT == 1
            figure; tiledlayout(2,1); nexttile([1 1]); 
            plot(TINT,SIGINT,'o',TINT,FITVECTOR,'-'); 
            nexttile([1 1]); plot(TINT,SIGINT-FITVECTOR); title('Fitting without crop')
        end
        % Trim padded signal to original length
        TINT = TINT(1:originallength);
        SIGINT = SIGINT(1:originallength);
        FITVECTOR = FITVECTOR(1:originallength);
        % Calculate pulse width from fit
        IRFwidth = FITVAL(7); % ps
        FWHM = IRFwidth/sqrt(2); % ps
        % Biexponential lifetime assignments
        LIFETIME1 = min(FITVAL(9),FITVAL(11)); % ps
        LIFETIME2 = max(FITVAL(9),FITVAL(11)); % ps
        % Lifetime amplitude ratio
        if abs(FITVAL(9)) < abs(FITVAL(11)) % if fitval(9) is shorter (τ1)
            LT1LT2 = FITVAL(8)/FITVAL(10); % A1/A2
        else % means fitval(9) is longer (τ2)
            LT1LT2 = FITVAL(10)/FITVAL(8); % still A1/A2
        end
        if CURRLTFITTYPE == 1 % Gauss*monoexp
            LIFETIME1 = FITVAL(9); LIFETIME2 = 0;
        end

        if CURRLTFITTYPE == 3 % Gauss*stretchexp
            LIFETIME1 = FITVAL(11); LIFETIME2 = 0;
            LT1LT2 = FITVAL(12); % stretching factor
        end
        if CURRLTFITTYPE == 4 % Gauss*strbiexp
            LIFETIME1 = FITVAL(9); LIFETIME2 = FITVAL(11);
            LT1LT2 = FITVAL(8)/FITVAL(10);
        end
        % Calculate residuals and ssresid
        RESID = SIGINT-FITVECTOR;
        RESID = spline(TINT,RESID,T);
        SRR = SIGPEAK./std(RESID(2:end));
        ssresid = sum(RESID.^2); % sum of squares of residuals
        tss = sum((SIGINT-mean(SIGINT)).^2); % total sum of squares
        R2 = 1-(ssresid/tss); % coefficient of determination (R^2)
        if R2 < 0 % remove nonsense R^2 values
            R2 = 0;
        end
        % FFT of residuals
        ts = TINT*1e-12; % time, s
        dt = abs(ts(1)-ts(2)); % time step, s
        fs = 1/dt; % freq BW, Hz
        resfft = fft(RESID); % raw FFT
        resfftshift = fftshift(resfft); % centered FFT
        FFTY = abs(resfftshift).^2/length(RESID); % power spectrum
        freqHz = (-length(RESID)/2+1/2:length(RESID)/2-1/2)*(fs/length(RESID)); % FFT frequency vector in Hz
        FFTX = freqHz/29979245800; % FFT freq in cm-1
        if TROUBLESHOOT == 1
            % Check fit
            figure; tiledlayout(3,1); nexttile([1 1]); plot(TINT,RESID); nexttile([2 1]); plot(T,SIGNAL,'o',TINT,FITVECTOR); title('Fitting')
            % Check FFT
            figure; plot(FFTX,FFTY); title('FFT of residuals');
        end
    else % no lifetime fitting
        TINT = T; SIGINT = SIGNAL; FITVECTOR = SIGNAL;
        RESID = SIGNAL; FFTX = T; FFTY = SIGNAL;
        LIFETIME1 = 0; LIFETIME2 = 0; LT1LT2 = 0;
        FITVAL = zeros(12,1); R2 = 0; FWHM = 0;
        originallength = length(T); SRR = 0;
    end
    % Brute-forced smoother fit curve
    TFIT = linspace(TINT(1),TINT(originallength),originallength*10-9).'; % smoother t vector
    FITCURVE = spline(TINT,FITVECTOR,TFIT);
    return
end
