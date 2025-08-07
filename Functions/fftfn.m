% Compute an FFT for a given time vector and data vector
function [FFTXHz,FFTY,FFTXcm] = fftfn(TIMESEC,DATA)
%%% This function computes the Fourier transform of a time-domain dataset,
% given a time input in seconds.
% (outputs are in cm-1 or Hz)

    % Ensure column vectors
    TIMESEC = ensurecolumn(TIMESEC); DATA = ensurecolumn(DATA);

    % Make sure lengths match
    [TIMESEC,DATA] = lengthmatch(TIMESEC,DATA);

    % Make sure time vector is evenly spaced
    deltat = zeros(length(TIMESEC)-1,1);
    for i=1:length(TIMESEC)-1
        deltat(i) = round(1e2*(TIMESEC(i+1)-TIMESEC(i)))/1e2; % round to nearest 0.01
    end
    % Check spacing
    tspacing = unique(deltat); % length 1 <-> evenly spaced t
    % Check t spacing (even) and length (odd)
    if isscalar(tspacing) && mod(length(TIMESEC),2) == 1
        TINT = TIMESEC; DATAINT = DATA; % data are good for fitting
    else % interpolate
        TINT = linspace(TIMESEC(1),TIMESEC(end),length(TIMESEC)*2-1).';
        DATAINT = spline(TIMESEC, DATA, TINT);
    end

    % Calculate FFT
    dt = abs(TINT(1)-TINT(2)); % time step, s
    fs = 1/dt; % freq bandwidth, Hz
    resfft = fft(DATAINT); % raw FFT
    resfftshift = fftshift(resfft); % centered FFT
    FFTY = abs(resfftshift).^2/length(DATAINT); % power spectrum
    FFTXHz = (-length(DATAINT)/2+1/2:length(DATAINT)/2-1/2)*(fs/length(DATAINT)); % FFT frequency vector in Hz
    FFTXcm = FFTXHz/29979245800; % FFT freq in cm-1

    % Crop FFT output to positive frequencies
    FFTY = FFTY(FFTXcm >= 0); % y-values (amplitude component of FFT)
    FFTXHz = FFTXHz(FFTXcm >= 0); % output in Hz
    FFTXcm = FFTXcm(FFTXcm >= 0); % output in cm-1
    return
end
