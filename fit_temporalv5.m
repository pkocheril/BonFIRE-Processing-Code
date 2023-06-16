%% Fit temporal sweep data
%%% v2 - fit to Gaussian*biexponential with optional beating
%%% v3 - added ability to fit linear background and compared 
% trust-region-reflective vs levenberg-marquardt (TRR is better)
%%% v4 - added spline interpolation for nonlinearly spaced data and
% tried Gaussian IRF pre-fitting, but it doesn't work
%%% v5 - cleaned up v4 and fixed FWHM calculation, also verified that
% Gaussian amplitude term isn't needed and doesn't affect lifetimes

% Load data
clear; close all;
data = importdata('779m-706mW-1770.8nm-68mW-ND1-PMT1-10kHz-250kHz.txt');
x = data(:,1);
DC = data(:,2);
AC = data(:,3);

% Converting to time
cmmps = 299792458*1E3*1E-12; % c in mm/ps
t = zeros([length(x),1]);
[corrACmax, t0index] = max(AC); % estimating t0 by AC peak (rough)
t(:) = (x(:)-x(t0index))*4/cmmps; % time vector in ps

% Fitting AC data
signal = AC;
sigamp = max(signal)-min(signal); % used for Gaussian amplitude

% Interpolate for convolution
samplingRateIncrease = 2;
tint = linspace(t(1), t(end), length(t)*samplingRateIncrease-1);
sigint = spline(t, signal, tint);
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

% Calculate residuals, ssresid, IRF FWHM, and lifetimes
resid = sigint-fitcurve;
ssresid = sum(resid.^2);
FWHM = sqrt(log(2)/fitval(2)); % pulse width, ps
lifetime1 = 1/fitval(4);
lifetime2 = 1/fitval(6);

% Plot fits
figure; tiledlayout(3,2); nexttile([1 1]); plot(base,t,signal,exclude1);
xlabel('Time (ps)'); ylabel('Signal (AU)'); nexttile([1 1]);
plot(t,signal,'r*'); hold on; plot(tint,sigint,'bo'); hold off
legend('Data','Interpolation'); xlabel('Time (ps)'); ylabel('Signal (AU)')
nexttile([1 2]); plot(t, signal, 'o', tint, fitcurve, '-');
xlabel('Time (ps)'); ylabel('Signal (AU)');
legend('Data','Fit'); nexttile([1 2]); plot(tint, resid)
legend('Residuals'); xlabel('Time (ps)'); ylabel('Residuals (AU)');

