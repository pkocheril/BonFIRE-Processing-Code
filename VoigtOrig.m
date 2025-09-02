% Create a Voigt lineshape for plotting
function [VOIGTY,VOIGTX] = VoigtOrig(HEIGHT,CENTER,GFWHM,LFWHM,X)
%%% This function creates a Voigt lineshape (x and y vectors) for an
% input height, center, Gaussian FWHM, Lorentzian FWHM, and x vector.
% X must be odd in length for VOIGTX to equal X.

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