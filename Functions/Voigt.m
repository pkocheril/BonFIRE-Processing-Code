% Create a Voigt lineshape for plotting from Fityk
function [VOIGTY,VOIGTX] = Voigt(HEIGHT,CENTER,GWIDTH,SHAPE,X)
%%% This function creates a Voigt lineshape (x and y vectors) for an
% input height, center, Gaussian width, shape, and x vector.
% Designed to take Fityk outputs as inputs.
% X must be odd in length for VOIGTX to equal X.

    % Convert GWIDTH and SHAPE to GFHWM and LFWHM
    GFWHM = GWIDTH*2*sqrt(log(2)); % Fityk "gwidth" term: GFWHM = gwidth * 2*sqrt(log(2)) = gwidth * 1.6651
    LFWHM = 2*GWIDTH*SHAPE; % Fityk "shape" term: LFWHM = 2*shape*gwidth

    % Make sure X is a column vector
    if length(X) == width(X)
        X = X.';
    end
    % Create X vector for convolution input (center around 0)
    if mod(length(X),2) == 0 % if length(X) is even
        X(end+1) = X(end)+(X(end)-X(end-1)); % extend by 1 point
    end
    INPUTX = linspace(0.5*(min(X)-CENTER),0.5*(max(X)-CENTER),0.5*(length(X)+1)).';
    % Make Gaussian
    Gstdev = GFWHM./(2*sqrt(2*log(2)));
    Gdecay = -1./(2*Gstdev.^2);
    GAUSSIAN = exp(Gdecay.*(INPUTX).^2); % centered at 0
    % Make Lorentzian
    gamma = LFWHM./2;
    LORENTZIAN = (gamma./pi).*((INPUTX).^2+gamma.^2).^-1; % centered at 0
    % Make Voigt
    VOIGTY = conv(GAUSSIAN,LORENTZIAN);
    VOIGTY = HEIGHT.*VOIGTY./max(VOIGTY);
    VOIGTX = linspace(min(X)-CENTER,max(X)-CENTER,length(X)).';
    VOIGTX = VOIGTX+CENTER;
    return
end
