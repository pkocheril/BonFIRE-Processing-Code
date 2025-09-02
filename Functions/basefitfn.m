% Baseline fitting
function [BASECURVE,BLEACHRATE,BASEFIT,LBBASE,UBBASE,SBR,CORRSIG,CORRSIGBASE,TBASE,SIGBASE,CURRBASEFITTYPE] = ...
    basefitfn(T,SIGNAL,CURRBASEFITTYPE,CUTLOW,CUTHIGH,TROUBLESHOOT)
%%% This function performs baseline fitting and can guess which function
% to use based on the variation of the data.

%%% v2 - fixed SBR calculation (previously was peak value instead of peak -
% baseline)

    % Setup baseline fitting
    ilow = ceil(length(T)*CUTLOW); % using ceiling to round up
    ihigh = ceil(length(T)*CUTHIGH);
    if isempty(CURRBASEFITTYPE)
        frontmean = mean(SIGNAL(1:ilow));
        tailmean = mean(SIGNAL(ihigh:end));
        fronttaildiff = 2*(frontmean-tailmean)/(frontmean+tailmean);
        if fronttaildiff > 0.03
            CURRBASEFITTYPE = 3;
        else
            CURRBASEFITTYPE = 1;
        end
    end
    if CURRBASEFITTYPE >= 2 % if exponential baseline
        startbase = 1; % start at first point
    else % if linear baseline
        startbase = 2; % start at second point
    end
    TBASE = [T(startbase:ilow); T(ihigh:end)]; % trimmed t
    SIGBASE = [SIGNAL(startbase:ilow); SIGNAL(ihigh:end)]; % trimmed signal
    SBR = (max(SIGNAL(ilow:ihigh))-mean(SIGNAL(ihigh:end)))/mean(SIGNAL(ihigh:end));
    % Baseline fitting
    if CURRBASEFITTYPE == 0 % no baseline fitting, set basecurve to 0
        BASECURVE = zeros(height(T),width(T));
        BLEACHRATE = 0;
        BASEFIT = zeros(1,5); LBBASE = BASEFIT; UBBASE = BASEFIT;
    else
        basefn = @(b) b(1)+b(2)*exp(-b(3)*(TBASE-b(4)))+b(5)*TBASE-SIGBASE; % error function
        LBBASE = [-Inf -Inf 0 -Inf -Inf];
        UBBASE = [Inf Inf Inf Inf Inf];
        if CURRBASEFITTYPE == 1
            LBBASE(2:4) = 0; UBBASE(2:4) = 0;
        end
        if CURRBASEFITTYPE == 2
            LBBASE(5) = 0; UBBASE(5) = 0;
        end
        baseopt=optimoptions(@lsqnonlin);
        baseopt.MaxFunctionEvaluations = 1e6; % let the fit run longer
        baseopt.MaxIterations = 1e6; % let the fit run longer
        baseopt.FunctionTolerance = abs(max(SIGNAL))*1e-12; % make the fit more accurate
        baseopt.OptimalityTolerance = abs(max(SIGNAL))*1e-12; % make the fit more accurate
        baseopt.StepTolerance = abs(max(SIGNAL))*1e-12; % make the fit more precise
        if TROUBLESHOOT == 0
            baseopt.Display = 'off'; % silence console output
        end
        basegs = [min(SIGNAL) 0.1*(max(SIGNAL)-min(SIGNAL)) 0.1 0 0];
        BASEFIT = lsqnonlin(basefn,basegs,LBBASE,UBBASE,baseopt);
        BASECURVE = BASEFIT(1)+BASEFIT(2)*exp(-BASEFIT(3)*(T-BASEFIT(4)))+BASEFIT(5)*T;
        if BASEFIT(2) > 0 && BASEFIT(3) > 0 
            BLEACHRATE = BASEFIT(3); % ps-1
        else % bleach rate doesn't have physical meaning
            BLEACHRATE = 0;
        end
        if TROUBLESHOOT == 1 % Check baseline fitting
            figure; plot(T,SIGNAL,'o',T,BASECURVE,'-'); title('Baseline fitting');
        end
    end
    CORRSIG = SIGNAL-BASECURVE;
    CORRSIGBASE = [CORRSIG(2:ilow); CORRSIG(ihigh:end)];
    return
end
