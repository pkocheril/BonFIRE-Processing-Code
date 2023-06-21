%% Fit temporal sweep image data
%%% initially written from fit_temporal_v3
%%% v2 - added lifetime parsing logic
%%% v3 - updated based on fit_temporal_v5
%%% v4 - updated based on fit_temporal_v6

clear; clc; close all; warning off;
% Load data
x = importdata("Tlist.txt");
data = double(tiffreadVolume("FOV1_64X64_Idler2159.9_DFG5377.6_P0.0657_CH2.tif"));
imagesize = size(data); % x by y by z -- x by y images, z is temporal

% Image to save to
tiffname = 'biexp_lifetimes_FOV1.tif';

% Test-run on subset of image
testsubset = 1; % 1 = do a test run, 0 = full processing
subsize = 2; % size of subset for test run

if testsubset == 1
    imagesize(1) = subsize; imagesize(2) = subsize; % pick subset of image
    showfits = 1; % show individual temporal fits
    savedatayn = 0; % don't save data
else
    showfits = 0; % don't show fits
    savedatayn = 1; % save data to tiff
end

% Matrices to write data to
FWHM = zeros([imagesize(1) imagesize(2)]);
lifetime2 = FWHM; ssresid = FWHM; lifetime1 = FWHM; tss = FWHM; r2 = FWHM;

% Estimating t0 index by peak
cmmps = 299792458*1E3*1E-12; % c in mm/ps
[roughmax, t0index] = max(data(1,1,:));
t = zeros([length(x),1]);
t(:) = (x(:)-x(t0index))*4/cmmps; % time vector in ps

% Interpolate for convolution
samplingRateIncrease = 2;
tint = linspace(t(1), t(end), length(t)*samplingRateIncrease-1);
tc = tint(1:2:end); % pulls odd values of tint for convolution

% Setup baseline fitting cutoffs - fit first 10% and last 5%
cutlow = 0.1; % 10th %ile
cuthigh = 0.05; % 95th %ile
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

for ii=1:imagesize(1)
    for jj=1:imagesize(2)
        % Set up vectors
        signal = squeeze(data(ii,jj,:)); % pulls sweep data as vector
        sigamp = max(signal)-min(signal); % used for Gaussian amplitude
        sigint = spline(t, signal, tint); % interpolate signal

        % Fit baseline
        base = fit(t,signal,'poly1','Exclude',exclude1); % fit baseline
        basecoef = coeffvalues(base); % save coeff values
        
        % Define error function as fit (parameter vector r) minus signal
        fun = @(r) conv(sigamp*exp(-r(2)*(tc-r(1)).^2), ...
            r(3)*exp(-r(4)*(tc-r(1)))+r(5)*exp(-r(6)*(tc-r(1)))+ ...
            r(7)*exp(-r(8)*(tc-r(1))).*cos(r(9)*tc-r(10))+ ...
            r(11)*exp(-r(12)*(tc-r(1))).*cos(r(13)*tc-r(14))+ ...
            r(15)*exp(-r(16)*(tc-r(1))).*cos(r(17)*tc-r(18)))+ ...
            r(19)*tint+r(20)-sigint;
        
        % Initial guesses for r
        % Gaussian terms
        cent = 0; % r(1), ps -- center
        gdec = 4*log(2)/(4^2); % r(2), ps^-2 -- Gaussian decay
        
        % Exponential terms
        eamp1 = 1e-5; % r(3) -- amplitude 1
        edec1 = 1/2; % r(4), ps^-1 -- 1/lifetime 1
        eamp11 = 1e-3; % r(5)  -- amplitude 1.1
        edec11 = 1/10; % r(6), ps^-1  -- 1/lifetime 1.1
        
        % Beating terms - can use 0.1, 0.1, 0.5, 0 for beating initial guesses
        eamp2 = 0; % r(7) -- amplitude 2 (beat 1)
        edec2 = 0; % r(8), ps^-1  -- lifetime 2 (beat 1)
        fr1 = 0; % r(9), rad ps^-1  -- beat freq 1
        phs1 = 0; % r(10), rad -- beat phase 1
        eamp3 = 0; % r(11) -- amplitude 3 (beat 2)
        edec3 = 0; % r(12), ps^-1 -- 1/lifetime 3 (beat 2)
        fr2 = 0; % r(13), rad ps^-1 -- beat freq 2
        phs2 = 0; % r(14), rad -- beat phase 2
        eamp4 = 0; % r(15) -- amplitude 4 (beat 3)
        edec4 = 0; % r(16), ps^-1 -- 1/lifetime 4 (beat 3)
        fr3 = 0; % r(17), rad ps^-1 -- beat freq 3
        phs3 = 0; % r(18), rad -- beat phase 3
        
        % Baseline terms
        slope = basecoef(1); % r(19), slope from baseline fit
        intercept = basecoef(2); % r(20), intercept from baseline fit
        
        % Fit setup
        r0 = [cent,gdec,eamp1,edec1,eamp11,edec11,... % Gaussian and exp decays
            eamp2,edec2,fr1,phs1,... % first beating term
            eamp3,edec3,fr2,phs2... % second beating term
            eamp4,edec4,fr3,phs3,... % third beating term
            slope,intercept];
        lb = [-10,0,0,0,0,0,... % lower bounds
            0,0,0,-180,...
            0,0,0,-180,...
            0,0,0,-180,...
            slope,intercept];
        ub = [20,9,9,9,9,9,... % upper bounds - eamp11=0 for monoexp
            0,9,9,180,... % set first value to 0 to mute beating
            0,9,9,180,... % set first value to 0 to mute beating
            0,9,9,180,... % set first value to 0 to mute beating
            slope,intercept];
        
        % Option to fit or just plot initial guesses
        fityn = 1; % 1 = do fitting, 0 = plot initial guesses
        
        % Option to use previous fit parameters as initial guesses
        %r0 = fitval; % uncomment to use previous fit
        
        if fityn == 1
            % Run lsqnonlin - minimizes error function by adjusting r
            options=optimoptions(@lsqnonlin);
            options.MaxFunctionEvaluations = 1e6; % let the fit run longer
            options.MaxIterations = 1e6; % let the fit run longer
            options.FunctionTolerance = 1e-12; % make the fit more accurate
            options.Display = 'off'; % silence console output
            fitval = lsqnonlin(fun,r0,lb,ub,options); % fit coeffs
        else
            % Plot initial guesses
            fitval = r0;
        end
        
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
        
        % Calculate residuals, ssresid, and IRF FWHMs
        resid = sigint-fitcurve;
        ssresid(ii,jj) = sum(resid.^2); % sum of squares of residuals
        tss(ii,jj) = sum((sigint-mean(sigint)).^2); % total sum of squares
        r2(ii,jj) = 1-(ssresid(ii,jj)/tss(ii,jj)); % coefficient of determination (R^2)
        FWHM(ii,jj) = sqrt(log(2)/fitval(2)); % pulse width, ps

        % Lifetime assignment logic - want lifetime1 to be shorter
        if 1/fitval(4) < 1/fitval(6)
            if 1/fitval(4) < 1e2 % only keep reasonable values < 100 ps
                lifetime1(ii,jj) = 1/fitval(4); % ps
            else
                lifetime1(ii,jj) = 0; %
            end
            if 1/fitval(6) < 1e2
                lifetime2(ii,jj) = 1/fitval(6); % ps
            else
                lifetime2(ii,jj) = 0; % ps
            end
        else % means 1/fitval(6) is shorter --> swap lifetime1 and 2
            if 1/fitval(4) < 1e2
                lifetime2(ii,jj) = 1/fitval(4); % ps
            else
                lifetime2(ii,jj) = 0; %
            end
            if 1/fitval(6) < 1e2
                lifetime1(ii,jj) = 1/fitval(6); % ps
            else
                lifetime1(ii,jj) = 0; % ps
            end
        end

        % Lifetime assignment logic 2 - constrained lifetime -> lifetime2
        if lb(6) == ub(6) % if only one lifetime floated, put in lifetime1
            if 1/fitval(4) < 1e2
                lifetime1(ii,jj) = 1/fitval(4); % ps
            else
                lifetime1(ii,jj) = 0; %
            end
            if 1/fitval(6) < 1e2
                lifetime2(ii,jj) = 1/fitval(6); % ps
            else
                lifetime2(ii,jj) = 0; % ps
            end
        end
        if lb(4) == ub(4) % if only one lifetime floated, put in lifetime1
            if 1/fitval(4) < 1e2 % only keep reasonable values < 100 ps
                lifetime2(ii,jj) = 1/fitval(4); % ps
            else
                lifetime2(ii,jj) = 0; %
            end
            if 1/fitval(6) < 1e2 % only keep reasonable values < 100 ps
                lifetime1(ii,jj) = 1/fitval(6); % ps
            else
                lifetime1(ii,jj) = 0; % ps
            end
        end
        
        if showfits == 1
            % Plot fits
            figure; tiledlayout(3,2); nexttile([1 1]); plot(base,t,signal,exclude1);
            xlabel('Time (ps)'); ylabel('Signal (AU)'); nexttile([1 1]);
            plot(t,signal,'r*'); hold on; plot(tint,sigint,'bo'); hold off
            legend('Data','Interpolation'); xlabel('Time (ps)'); ylabel('Signal (AU)')
            nexttile([1 2]); plot(t, signal, 'o', tint, fitcurve, '-');
            xlabel('Time (ps)'); ylabel('Signal (AU)');
            legend('Data','Fit'); nexttile([1 2]); plot(tint, resid)
            legend('Residuals'); xlabel('Time (ps)'); ylabel('Residuals (AU)');
        end
    end
end

if savedatayn == 1
    % Save data
    Save_tiff(tiffname,lifetime1,lifetime2,ssresid,FWHM,r2)
else
    disp('Test run completed');
end

%% Save data manually
Save_tiff(tiffname,lifetime1,lifetime2,ssresid,FWHM,r2)

%% Inspecting problematic fits
jj = 17; % x-value in Fiji
ii = 5; % y-value in Fiji

% Set up vectors
signal = squeeze(data(ii,jj,:)); % pulls sweep data as vector
sigamp = max(signal)-min(signal); % used for Gaussian amplitude
sigint = spline(t, signal, tint); % interpolate signal

% Fit baseline
base = fit(t,signal,'poly1','Exclude',exclude1); % fit baseline
basecoef = coeffvalues(base); % save coeff values

% Define error function as fit (parameter vector r) minus signal
fun = @(r) conv(sigamp*exp(-r(2)*(tc-r(1)).^2), ...
    r(3)*exp(-r(4)*(tc-r(1)))+r(5)*exp(-r(6)*(tc-r(1)))+ ...
    r(7)*exp(-r(8)*(tc-r(1))).*cos(r(9)*tc-r(10))+ ...
    r(11)*exp(-r(12)*(tc-r(1))).*cos(r(13)*tc-r(14))+ ...
    r(15)*exp(-r(16)*(tc-r(1))).*cos(r(17)*tc-r(18)))+ ...
    r(19)*tint+r(20)-sigint;

% Initial guesses for r
% Gaussian terms
cent = 0; % r(1), ps -- center
gdec = 4*log(2)/(4^2); % r(2), ps^-2 -- Gaussian decay

% Exponential terms
eamp1 = 1e-5; % r(3) -- amplitude 1
edec1 = 1/2; % r(4), ps^-1 -- 1/lifetime 1
eamp11 = 1e-3; % r(5)  -- amplitude 1.1
edec11 = 1/10; % r(6), ps^-1  -- 1/lifetime 1.1

% Beating terms - can use 0.1, 0.1, 0.5, 0 for beating initial guesses
eamp2 = 0; % r(7) -- amplitude 2 (beat 1)
edec2 = 0; % r(8), ps^-1  -- lifetime 2 (beat 1)
fr1 = 0; % r(9), rad ps^-1  -- beat freq 1
phs1 = 0; % r(10), rad -- beat phase 1
eamp3 = 0; % r(11) -- amplitude 3 (beat 2)
edec3 = 0; % r(12), ps^-1 -- 1/lifetime 3 (beat 2)
fr2 = 0; % r(13), rad ps^-1 -- beat freq 2
phs2 = 0; % r(14), rad -- beat phase 2
eamp4 = 0; % r(15) -- amplitude 4 (beat 3)
edec4 = 0; % r(16), ps^-1 -- 1/lifetime 4 (beat 3)
fr3 = 0; % r(17), rad ps^-1 -- beat freq 3
phs3 = 0; % r(18), rad -- beat phase 3

% Baseline terms
slope = basecoef(1); % r(19), slope from baseline fit
intercept = basecoef(2); % r(20), intercept from baseline fit

% Fit setup
r0 = [cent,gdec,eamp1,edec1,eamp11,edec11,... % Gaussian and exp decays
    eamp2,edec2,fr1,phs1,... % first beating term
    eamp3,edec3,fr2,phs2... % second beating term
    eamp4,edec4,fr3,phs3,... % third beating term
    slope,intercept];
lb = [-10,0,0,0,0,0,... % lower bounds
    0,0,0,-180,...
    0,0,0,-180,...
    0,0,0,-180,...
    slope,intercept];
ub = [20,9,9,9,9,9,... % upper bounds - eamp11=0 for monoexp
    0,9,9,180,... % set first value to 0 to mute beating
    0,9,9,180,... % set first value to 0 to mute beating
    0,9,9,180,... % set first value to 0 to mute beating
    slope,intercept];

% Option to fit or just plot initial guesses
fityn = 1; % 1 = do fitting, 0 = plot initial guesses

% Option to use previous fit parameters as initial guesses
%r0 = fitval; % uncomment to use previous fit

if fityn == 1
    % Run lsqnonlin - minimizes error function by adjusting r
    options=optimoptions(@lsqnonlin);
    options.MaxFunctionEvaluations = 1e6; % let the fit run longer
    options.MaxIterations = 1e6; % let the fit run longer
    options.FunctionTolerance = 1e-12; % make the fit more accurate
    options.Display = 'off'; % silence console output
    fitval = lsqnonlin(fun,r0,lb,ub,options); % fit coeffs
else
    % Plot initial guesses
    fitval = r0;
end

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
        
% Calculate residuals, ssresid, and IRF FWHMs
resid = sigint-fitcurve;
ssresidind = sum(resid.^2);
tssind = sum((sigint-mean(sigint)).^2); % total sum of squares
r2ind = 1-(ssresidind/tssind); % coefficient of determination (R^2)
FWHMind = sqrt(log(2)/fitval(2)); % pulse width, ps

% Lifetime assignment logic - want lifetime1 to be shorter
if 1/fitval(4) < 1/fitval(6)
    if 1/fitval(4) < 1e2 % only keep reasonable values < 100 ps
        lifetime1ind = 1/fitval(4); % ps
    else
        lifetime1ind = 0; %
    end
    if 1/fitval(6) < 1e2
        lifetime2ind = 1/fitval(6); % ps
    else
        lifetime2ind = 0; % ps
    end
else % means 1/fitval(6) is shorter --> swap lifetime1 and 2
    if 1/fitval(4) < 1e2
        lifetime2ind = 1/fitval(4); % ps
    else
        lifetime2ind = 0; %
    end
    if 1/fitval(6) < 1e2
        lifetime1ind = 1/fitval(6); % ps
    else
        lifetime1ind = 0; % ps
    end
end

% Lifetime assignment logic 2 - constrained lifetime -> lifetime2
if lb(6) == ub(6) % if only one lifetime floated, put in lifetime1
    if 1/fitval(4) < 1e2
        lifetime1ind = 1/fitval(4); % ps
    else
        lifetime1ind = 0; %
    end
    if 1/fitval(6) < 1e2
        lifetime2ind = 1/fitval(6); % ps
    else
        lifetime2ind = 0; % ps
    end
end
if lb(4) == ub(4) % if only one lifetime floated, put in lifetime1
    if 1/fitval(4) < 1e2 % only keep reasonable values < 100 ps
        lifetime2ind = 1/fitval(4); % ps
    else
        lifetime2ind = 0; %
    end
    if 1/fitval(6) < 1e2 % only keep reasonable values < 100 ps
        lifetime1ind = 1/fitval(6); % ps
    else
        lifetime1ind = 0; % ps
    end
end

% Inspect fit
figure; tiledlayout(3,2); nexttile([1 1]); plot(base,t,signal,exclude1);
xlabel('Time (ps)'); ylabel('Signal (AU)'); nexttile([1 1]);
plot(t,signal,'r*'); hold on; plot(tint,sigint,'bo'); hold off
legend('Data','Interpolation'); xlabel('Time (ps)'); ylabel('Signal (AU)')
nexttile([1 2]); plot(t, signal, 'o', tint, fitcurve, '-');
xlabel('Time (ps)'); ylabel('Signal (AU)');
legend('Data','Fit'); nexttile([1 2]); plot(tint, resid)
legend('Residuals'); xlabel('Time (ps)'); ylabel('Residuals (AU)');


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


