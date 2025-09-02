% Basic lifetime fitting
function [TFIT,FITCURVE,FITVAL,RESID] = basicltfit(T,SIGNAL,CURRLTFITTYPE)
%%% This function performs lifetime fitting but without the extra "smarts"
% of the previous function.
    
    [sigmax,maxind] = max(SIGNAL);
    % Ensure column vectors
    if length(T) == width(T)
        T = T.';
    end
    if length(SIGNAL) == width(SIGNAL)
        SIGNAL = SIGNAL.';
    end
    % Make sure lengths match
    if length(T) ~= length(SIGNAL)
        if length(T) < length(SIGNAL) % SIGNAL is longer
            SIGNAL = SIGNAL(1:length(T));
        else % T is longer
            T = T(1:length(SIGNAL));
        end
    end
    
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
    siglong = SIGINT(end).*ones(length(tlong),1);
    siglong(1:length(SIGINT)) = SIGINT;
    
    % Signal padding for short-tailed data
    TINT = tlong;
    SIGINT = siglong;

    % Prepare convolution input (same spacing, half length, symmetrically around zero)
    tc = linspace(0.5*min(TINT),0.5*max(TINT),0.5*(length(TINT)+1)).';
    fun = @(r) r(1)+r(2)*exp(-r(3)*(TINT-r(4)))+r(5)*TINT+ ... % baseline
        conv(exp(-((tc-r(6))./(r(7)/(2*sqrt(log(2))))).^2), ...
        heaviside(tc).*(r(8)*exp(-(tc)./r(9))+r(10)*exp(-((tc)./(r(11))).^r(12)))) - SIGINT;
    % Initial guesses
    r0 = [0,0,0,0,0,... % basecoefs
        T(maxind),2.6,... % IRF center (ps), IRF width (ps)
        1,1,0.02,10,1]; % amp 1, τ1 (ps), amp 2, τ2 (ps), β (stretch)
    % Lower and upper bounds
    lb = [-Inf,-Inf,0,-Inf,-Inf,... % basecoefs
        min(T),0.1,... % IRF center (ps), IRF min (ps)
        0,0.1,... % amp 1, τ1min (ps)
        0,0.1,... % amp 2, τ2min (ps)
        1e-2]; % β min
    ub = [Inf,Inf,0,Inf,Inf,... % basecoefs, 
        max(T),3*sqrt(2),... % IRF center (ps), IRF max (ps)
        10*sigmax,100,... % amp 1, τ1max (ps)
        10*sigmax,100,... % amp 2, τ2max (ps)
        1e2]; % β max
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
    
    % Run lsqnonlin - minimizes error function by adjusting r
    options=optimoptions(@lsqnonlin);
    options.MaxFunctionEvaluations = 1e6; % let the fit run longer
    options.MaxIterations = 1e6; % let the fit run longer
    options.FunctionTolerance = abs(max(SIGINT)-min(SIGINT))*1e-12; % make the fit more accurate
    options.OptimalityTolerance = abs(max(SIGINT)-min(SIGINT))*1e-12; % make the fit more accurate
    options.StepTolerance = abs(max(SIGINT)-min(SIGINT))*1e-12; % make the fit more precise
    options.Display = 'off'; % silence console output
    FITVAL = lsqnonlin(fun,r0,lb,ub,options); % fit coeffs
    FITVECTOR = FITVAL(1)+FITVAL(2)*exp(-FITVAL(3)*(TINT-FITVAL(4)))+FITVAL(5)*TINT+ ... % baseline
        conv(exp(-((tc-FITVAL(6))./(FITVAL(7)/(2*sqrt(log(2))))).^2), ...
        heaviside(tc).*(FITVAL(8)*exp(-(tc)./FITVAL(9))+FITVAL(10)*exp(-((tc)./(FITVAL(11))).^FITVAL(12))));

    % Trim padded vectors to original length
    TINT = TINT(1:originallength);
    SIGINT = SIGINT(1:originallength);
    FITVECTOR = FITVECTOR(1:originallength);

    RESID = SIGINT-FITVECTOR;
    RESID = spline(TINT,RESID,T);

    % Brute-forced smoother fit curve
    TFIT = linspace(TINT(1),TINT(end),length(TINT)*10-9).'; % smoother t vector
    FITCURVE = spline(TINT,FITVECTOR,TFIT);
    return
end
