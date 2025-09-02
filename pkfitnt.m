% Fitting cell-silent BonFIRE spectra with a Fano fit for NDR-TPA
function [NX,NY,SOLN,FWHM] = pkfitnt(XVAL,YVAL)
%%% This function fits y(x) to a Gaussian with a linear background and a
% Fano resonance.
    
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

    % Set up fitting function
    fn = @(x) x(1)+x(2)*XVAL+... % linear baseline
        x(3).*exp(-x(4).*(XVAL-x(5)).^2)+... % Gaussian (BonFIRE)
        x(6).*((x(7).*x(8)+XVAL-x(5)).^2)./(x(8).^2+(XVAL-x(5)).^2)+... % Fano (NDR-TPA)
        -1.*YVAL;
    
    % Set options
    opt=optimoptions(@lsqnonlin);
    opt.MaxFunctionEvaluations = 1e6; % let the fit run longer
    opt.MaxIterations = 1e6; % let the fit run longer
    opt.FunctionTolerance = abs(max(YVAL)-min(YVAL)).*1e-12; % make the fit more accurate
    opt.OptimalityTolerance = abs(max(YVAL)-min(YVAL)).*1e-12; % make the fit more accurate
    opt.StepTolerance = abs(max(YVAL)-min(YVAL)).*1e-12; % make the fit more precise
    opt.Display = 'off'; % silence console output
    
    [~,maxind] = max(YVAL);

    % Set initial guesses
    x0 = [-1.3241 7.7256e-04 ...
        0.6148 0.0096 XVAL(maxind) ...
        -1 .44 10];
    lb = [-Inf -Inf ...
        0  0.0025 -Inf ... % 0.0025 --> max FWHM = 33.3 cm-1
        -5 -1 -Inf];
    ub = [Inf Inf ...
        Inf 0.04 Inf ... % 0.04 --> min FWHM = 8.3 cm-1
        5 1 Inf];
    
    % Calculate fit
    SOLN = lsqnonlin(fn,x0,lb,ub,opt);
    FWHM = 2*sqrt(2*log(2))*(sqrt(1./(2*SOLN(4))));
    
    % Calculate function for plotting
    NX = min(XVAL):0.1:max(XVAL);
    NY = SOLN(1)+SOLN(2)*NX+SOLN(3).*exp(-SOLN(4).*(NX-SOLN(5)).^2)+...
        SOLN(6).*((SOLN(7).*SOLN(8)+NX-SOLN(5)).^2)./(SOLN(8).^2+(NX-SOLN(5)).^2);
    return
end