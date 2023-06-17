%% Fit temporal sweep image data
%%% initially written from fit_temporal_v3

clear; clc; close all; warning off;
% Load data
x = importdata("Tlist.txt");
data = double(tiffreadVolume("FOV1_64X64_Idler2159.9_DFG5377.6_P0.0657_CH2.tif"));
% x by y by z -- x by y images, z is temporal
imagesize = size(data);

% Estimating t0 index by AC peak (rough)
cmmps = 299792458*1E3*1E-12; % c in mm/ps
[roughmax, t0index] = max(data(1,1,:));

% Matrices to write lifetimes to
lifetime1 = zeros([imagesize(1) imagesize(2)]);
%lifetime1 = zeros([3,3]); % for testing with 3x3 subset of image
lifetime2 = lifetime1;
ssresid = lifetime1;
FWHM = lifetime1;

for ii=1:imagesize(1)
    for jj=1:imagesize(2)
        % Set up vectors
        t = zeros([length(x),1]);
        t(:) = (x(:)-x(t0index))*4/cmmps; % time vector in ps
        signal = squeeze(data(ii,jj,:)); % pulls sweep data as vector
        
        % Baseline fitting - fit first 10% and last 5%
        cut = length(t)*0.05;
        tlow = t(2*ceil(cut)); % xlow as 10th %ile
        thigh = t(end-ceil(cut)); % xhigh as 95th %ile
        
        % Setup excluded points for baseline fitting
        excludeDupes = zeros(size(t));
        for i = 1:length(t)
            if t(i) > tlow && t(i) < thigh
                excludeDupes(i) = i;
            end
        end
        
        % Remove duplicates from exclude list
        excludeunique = unique(excludeDupes(:).'); % adding the .' transposes it into a row vector
        exclude1 = nonzeros(excludeunique); % isolate nonzero data positions to exclude
        
        % Fit baseline
        base = fit(t,signal,'poly1','Exclude',exclude1); % fit baseline
        basecoef = coeffvalues(base); % save coeff values
                
        % Make odd-length tconv for convolution
        if mod(length(t),2) == 0 % if t is even
            t(length(t)+1) = t(end) + abs(t(end)-t(end-1)); % add a value to t
            signal(length(signal)+1) = signal(end); % add a value to signal
        end
        tconv = t(1:2:end); % pulls odd values of t
        
        % Define error function as fit (parameter vector r) minus signal
        fun = @(r) conv(exp(-r(2)*(tconv-r(1)).^2), ...
            r(3)*exp(-r(4)*(tconv-r(1)))+r(5)*exp(-r(6)*(tconv-r(1)))+ ...
            r(7)*exp(-r(8)*(tconv-r(1))).*cos(r(9)*tconv-r(10))+ ...
            r(11)*exp(-r(12)*(tconv-r(1))).*cos(r(13)*tconv-r(14))+ ...
            r(15)*exp(-r(16)*(tconv-r(1))).*cos(r(17)*tconv-r(18)))+ ...
            r(19)*t+r(20)-signal;
        
        % Initial guesses for r
        cent = 0; % r(1), ps -- center
        gwid = 0.1; % r(2), ps^-2 -- Gaussian width
        eamp1 = 1e-3; % r(3) -- amplitude 1
        edec1 = 1/2; % r(4), ps^-1 -- 1/lifetime 1
        eamp11 = 1e-4; % r(5)  -- amplitude 1.1
        edec11 = 1/10; % r(6), ps^-1  -- 1/lifetime 1.1
        slope = basecoef(1); % slope from baseline fit
        intercept = basecoef(2); % intercept from baseline fit
        
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
        
        % Fit setup
        r0 = [cent,gwid,eamp1,edec1,eamp11,edec11,... % Gaussian and exp decays
            eamp2,edec2,fr1,phs1,... % first beating term
            eamp3,edec3,fr2,phs2... % second beating term
            eamp4,edec4,fr3,phs3,... % third beating term
            slope,intercept];
        lb = [-10,0,0,0,0,0,... % lower bounds
            0,0,0,-180,...
            0,0,0,-180,...
            0,0,0,-180,...
            basecoef(1),basecoef(2)];
        ub = [20,9,9,9,9,9,... % upper bounds
            0,9,9,180,... % set first value to 0 to mute beating
            0,9,9,180,... % set first value to 0 to mute beating
            0,9,9,180,... % set first value to 0 to mute beating
            basecoef(1),basecoef(2)];
        
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
        fitcurve = conv(exp(-fitval(2)*(tconv-fitval(1)).^2), ...
            fitval(3)*exp(-fitval(4)*(tconv-fitval(1)))+ ...
            fitval(5)*exp(-fitval(6)*(tconv-fitval(1)))+ ...
            fitval(7)*exp(-fitval(8)*(tconv-fitval(1))).* ...
            cos(fitval(9)*tconv-fitval(10))+ ...
            fitval(11)*exp(-fitval(12)*(tconv-fitval(1))).* ...
            cos(fitval(13)*tconv-fitval(14))+ ...
            fitval(15)*exp(-fitval(16)*(tconv-fitval(1))).* ...
            cos(fitval(17)*tconv-fitval(18)))+...
            fitval(19)*t+fitval(20);
        
        % Calculate residuals, lifetimes, and Gaussian FWHMs
        resid = signal-fitcurve;
        ssresid(ii,jj) = sum(resid.^2);
        if 1/fitval(4) < 1e3
            lifetime1(ii,jj) = 1/fitval(4); % ps
        else
            lifetime1(ii,jj) = 0; %
        end
        if 1/fitval(6) < 1e3
            lifetime2(ii,jj) = 1/fitval(6); % ps
        else
            lifetime2(ii,jj) = 0; % ps
        end
        FWHM(ii,jj) = 2*log(2)*sqrt(1/(2*fitval(2))); % ps
    end
end

% Save data
Save_tiff('lifetimes.tif',lifetime1,lifetime2,ssresid,FWHM,[])

%% Checking out problematic fits
ii = 49; % y-value in Fiji
jj = 7; % x-value in Fiji

t = zeros([length(x),1]);
t(:) = (x(:)-x(t0index))*4/cmmps; % time vector in ps
signal = squeeze(data(ii,jj,:)); % pulls sweep data as vector

% Baseline fitting - fit first 10% and last 5%
cut = length(t)*0.05;
tlow = t(2*ceil(cut)); % xlow as 10th %ile
thigh = t(end-ceil(cut)); % xhigh as 95th %ile

% Setup excluded points for baseline fitting
excludeDupes = zeros(size(t));
for i = 1:length(t)
    if t(i) > tlow && t(i) < thigh
        excludeDupes(i) = i;
    end
end

% Remove duplicates from exclude list
excludeunique = unique(excludeDupes(:).'); % adding the .' transposes it into a row vector
exclude1 = nonzeros(excludeunique); % isolate nonzero data positions to exclude

% Fit baseline
base = fit(t,signal,'poly1','Exclude',exclude1); % fit baseline
basecoef = coeffvalues(base); % save coeff values
        
% Make odd-length tconv for convolution
if mod(length(t),2) == 0 % if t is even
    t(length(t)+1) = t(end) + abs(t(end)-t(end-1)); % add a value to t
    signal(length(signal)+1) = signal(end); % add a value to signal
end
tconv = t(1:2:end); % pulls odd values of t

% Define error function as fit (parameter vector r) minus signal
fun = @(r) conv(exp(-r(2)*(tconv-r(1)).^2), ...
    r(3)*exp(-r(4)*(tconv-r(1)))+r(5)*exp(-r(6)*(tconv-r(1)))+ ...
    r(7)*exp(-r(8)*(tconv-r(1))).*cos(r(9)*tconv-r(10))+ ...
    r(11)*exp(-r(12)*(tconv-r(1))).*cos(r(13)*tconv-r(14))+ ...
    r(15)*exp(-r(16)*(tconv-r(1))).*cos(r(17)*tconv-r(18)))+ ...
    r(19)*t+r(20)-signal;

% Initial guesses for r
cent = 0; % r(1), ps -- center
gwid = 0.1; % r(2), ps^-2 -- Gaussian width
eamp1 = 1e-3; % r(3) -- amplitude 1
edec1 = 1/2; % r(4), ps^-1 -- 1/lifetime 1
eamp11 = 1e-4; % r(5)  -- amplitude 1.1
edec11 = 1/10; % r(6), ps^-1  -- 1/lifetime 1.1
slope = basecoef(1); % slope from baseline fit
intercept = basecoef(2); % intercept from baseline fit

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

% Fit setup
r0 = [cent,gwid,eamp1,edec1,eamp11,edec11,... % Gaussian and exp decays
    eamp2,edec2,fr1,phs1,... % first beating term
    eamp3,edec3,fr2,phs2... % second beating term
    eamp4,edec4,fr3,phs3,... % third beating term
    slope,intercept];
lb = [-10,0,0,0,0,0,... % lower bounds
    0,0,0,-180,...
    0,0,0,-180,...
    0,0,0,-180,...
    basecoef(1),basecoef(2)];
ub = [20,9,9,9,9,9,... % upper bounds
    0,9,9,180,... % set first value to 0 to mute beating
    0,9,9,180,... % set first value to 0 to mute beating
    0,9,9,180,... % set first value to 0 to mute beating
    basecoef(1),basecoef(2)];

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
fitcurve = conv(exp(-fitval(2)*(tconv-fitval(1)).^2), ...
    fitval(3)*exp(-fitval(4)*(tconv-fitval(1)))+ ...
    fitval(5)*exp(-fitval(6)*(tconv-fitval(1)))+ ...
    fitval(7)*exp(-fitval(8)*(tconv-fitval(1))).* ...
    cos(fitval(9)*tconv-fitval(10))+ ...
    fitval(11)*exp(-fitval(12)*(tconv-fitval(1))).* ...
    cos(fitval(13)*tconv-fitval(14))+ ...
    fitval(15)*exp(-fitval(16)*(tconv-fitval(1))).* ...
    cos(fitval(17)*tconv-fitval(18)))+...
    fitval(19)*t+fitval(20);

% Calculate residuals, lifetimes, and Gaussian FWHM
resid = signal-fitcurve;
lifetime1(ii,jj) = 1/fitval(4); % ps
lifetime2(ii,jj) = 1/fitval(6); % ps
FWHM = 2*log(2)*sqrt(1/(2*fitval(2))); % ps

% Plot data & fit with residuals below
figure; tiledlayout(2,2); nexttile([1 2]) % top panel
hold on; plot(t,signal,'o'); plot(t,fitcurve,'-'); hold off
legend('Data','Fit'); xlabel('Time (ps)'); ylabel('Signal (AU)')
nexttile([1 2]); plot(t,resid) % bottom panel        
legend('Residuals'); xlabel('Time (ps)'); ylabel('Residuals (AU)')


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


