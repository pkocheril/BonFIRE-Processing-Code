% Fitting peaks to Gaussians
function [GX,GY,SOLN,FWHMS] = pkfit(XVAL,YVAL,nGauss,nPoly,peakorient)
%%% This function fits multi-peak data to a configurable number of
% Gaussians (1-4) with a variable-order polynomial baseline (1-4).

    % Ensure column vectors
    if length(XVAL) == width(XVAL)
        XVAL = XVAL.';
    end
    if length(YVAL) == width(YVAL)
        YVAL = YVAL.';
    end
    % Make sure lengths match
    if length(XVAL) ~= length(YVAL)
        if length(XVAL) < length(YVAL) % YVAL is longer
            YVAL = YVAL(1:length(XVAL));
        else % XVAL is longer
            XVAL = XVAL(1:length(YVAL));
        end
    end

    if isempty(nGauss)
        nGauss = 1; % guess single peak if not specified
    end
    if isempty(nPoly)
        nPoly = 1; % guess linear baseline if not specified
    end

    % Multi-Gaussian fit with polynomial baseline
    GX = linspace(min(XVAL),max(XVAL),length(XVAL).*10); GX = GX.';
    fn = @(x) x(1)+x(2)*XVAL+x(3).*XVAL.^2+x(4).*XVAL.^3+x(5).*XVAL.^4+... % polynomial baseline
        x(6).*exp(-x(7).*(XVAL-x(8)).^2)+... % Gauss1
        x(9).*exp(-x(10).*(XVAL-x(11)).^2)+... % Gauss2
        x(12).*exp(-x(13).*(XVAL-x(14)).^2)+... % Gauss3
        x(15).*exp(-x(16).*(XVAL-x(17)).^2)+... % Gauss4
        -1.*YVAL;
    opt=optimoptions(@lsqnonlin);
    opt.MaxFunctionEvaluations = 1e6; % let the fit run longer
    opt.MaxIterations = 1e6; % let the fit run longer
    opt.FunctionTolerance = (max(YVAL)-min(YVAL)).*1e-12; % make the fit more accurate
    opt.OptimalityTolerance = (max(YVAL)-min(YVAL)).*1e-12; % make the fit more accurate
    opt.StepTolerance = (max(YVAL)-min(YVAL)).*1e-12; % make the fit more precise
    opt.Display = 'off'; % silence console output

    if max(XVAL) < 5000 && min(XVAL) > 825 % IR sweep
        gwidth = 1.5e-2; % guess ~13 cm-1
        gmax = 1e-1;
        gmin = 1e-3;
    else % IR sweep
        gwidth = 1e-5; % guess ~500 cm-1
        gmax = Inf; gmin = -Inf;
    end

    % Figuring out initial guesses
    [YMAX,maxind] = max(YVAL);
    TEMPRES1 = YVAL-exp(-gwidth.*(XVAL-XVAL(maxind)).^2);
    [YMAX1,maxind1] = max(TEMPRES1);
    TEMPRES2 = TEMPRES1-exp(-gwidth.*(XVAL-XVAL(maxind1)).^2);
    [YMAX2,maxind2] = max(TEMPRES2);
    TEMPRES3 = TEMPRES2-exp(-gwidth.*(XVAL-XVAL(maxind2)).^2);
    [YMAX3,maxind3] = max(TEMPRES3);
    % Ylength = round(length(YVAL)./3);
    % [YMAX1,maxind1] = max(YVAL(1:Ylength));
    % [YMAX2,maxind2] = max(YVAL(Ylength:2*Ylength));
    % [YMAX3,maxind3] = max(YVAL(2*Ylength:end));

    x0 = [min(YVAL) (YVAL(end)-YVAL(1))./(XVAL(end)-XVAL(1)) 0 0 0 ...
        YMAX-min(YVAL) gwidth XVAL(maxind) ...
        YMAX1 gwidth XVAL(maxind1)...
        YMAX2 gwidth XVAL(maxind2)...
        YMAX3 gwidth XVAL(maxind3)];
    lb = [-Inf -Inf -Inf -Inf -Inf ...
        -Inf gmin min(XVAL) ...
        -Inf gmin min(XVAL) ...
        -Inf gmin min(XVAL) ...
        -Inf gmin min(XVAL)];
    ub = [Inf Inf Inf Inf Inf ...
        10*YMAX gmax max(XVAL) ...
        10*YMAX gmax max(XVAL) ...
        10*YMAX gmax max(XVAL) ...
        10*YMAX gmax max(XVAL)];

    % Set Gaussian orientation
    if strcmp(peakorient,'positive')
        lb(6:3:15) = 0; % no negative peaks
    end
    if strcmp(peakorient,'negative')
        ub(6:3:15) = 0; % no negative peaks
    end

    % Set number of Gaussians (can use 0 for polynomial fit)
    nGauss = round(nGauss);
    if nGauss < 4 && nGauss >= 0
        % 3 -> 15, 2 -> 12, 1 -> 9, 0 -> 6
        lb(((nGauss*3)+6):3:15) = 0; ub(((nGauss*3)+6):3:15) = 0;
    end

    % Set degree of polynomial (can use -1 for bare Gaussians)
    nPoly = round(nPoly);
    if nPoly < 4 && nPoly >= -1
        lb((nPoly+2):5) = 0; ub((nPoly+2):5) = 0; % 1 -> 3:5 zeroed, etc.
    end

    SOLN = lsqnonlin(fn,x0,lb,ub,opt);
    GY = SOLN(1)+SOLN(2)*GX+SOLN(3).*GX.^2+SOLN(4).*GX.^3+SOLN(5).*GX.^4+... % polynomial baseline
        SOLN(6).*exp(-SOLN(7).*(GX-SOLN(8)).^2)+... % Gauss1
        SOLN(9).*exp(-SOLN(10).*(GX-SOLN(11)).^2)+... % Gauss2
        SOLN(12).*exp(-SOLN(13).*(GX-SOLN(14)).^2)+... % Gauss3
        SOLN(15).*exp(-SOLN(16).*(GX-SOLN(17)).^2); % Gauss4

    stdev = sqrt(1./(2*SOLN(7))); FWHMS(1) = 2*sqrt(2*log(2))*stdev;
    stdev = sqrt(1./(2*SOLN(10))); FWHMS(2) = 2*sqrt(2*log(2))*stdev;
    stdev = sqrt(1./(2*SOLN(13))); FWHMS(3) = 2*sqrt(2*log(2))*stdev;
    stdev = sqrt(1./(2*SOLN(16))); FWHMS(4) = 2*sqrt(2*log(2))*stdev;
    return
end