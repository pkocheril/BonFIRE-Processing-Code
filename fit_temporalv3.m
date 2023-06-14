%% Fit temporal sweep data
%%% v2 - fit to Gaussian*biexponential with optional beating
%%% v3 - added ability to fit linear background and compared 
% trust-region-reflective vs levenberg-marquardt (TRR is better)

% Load data
data = importdata('779m-706mW-1770.8nm-68mW-ND1-PMT1-10kHz-250kHz.txt');
x = data(:,1);
DC = data(:,2);
AC = data(:,3);

% Converting to time
cmmps = 299792458*1E3*1E-12; % c in mm/ps
t = zeros([length(x),1]);
[corrACmax, t0index] = max(AC); % estimating t0 by AC peak (rough)
t(:) = (x(:)-x(t0index))*4/cmmps; % time vector in ps

% Baseline fitting - fit first 10% and last 5%
cut = length(t)*0.05;
tlow = t(2*ceil(cut)); % xlow as 10th %ile
thigh = t(end-floor(cut)); % xhigh as 95th %ile

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

% Fit baselines with lines
baseDC = fit(t,DC,'poly1','Exclude',exclude1);
baseAC = fit(t,AC,'poly1','Exclude',exclude1);

% Save coefficient values from fit to vector
basecoefDC = coeffvalues(baseDC);
basecoefAC = coeffvalues(baseAC);

% Convolution fitting
signal = AC;

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
slope = basecoefAC(1); % slope from baseline fit
intercept = basecoefAC(2); % intercept from baseline fit

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
    basecoefAC(1),basecoefAC(2)];
ub = [20,9,9,9,9,9,... % upper bounds
    0,9,9,180,... % set first value to 0 to mute beating
    0,9,9,180,... % set first value to 0 to mute beating
    0,9,9,180,... % set first value to 0 to mute beating
    basecoefAC(1),basecoefAC(2)];

% Option to fit or just plot initial guesses
fityn = 1; % 1 = do fitting, 0 = plot initial guesses

% Option to use previous fit parameters as initial guesses
%r0 = fitval; % uncomment to use previous fit

if fityn == 1
    % Run lsqnonlin - minimizes error function by adjusting r values
    options=optimoptions(@lsqnonlin);
    options.MaxFunctionEvaluations = 1e6; % let the fit run longer
    options.MaxIterations = 1e6; % let the fit run longer
    options.FunctionTolerance = 1e-12; % make the fit be more accurate
    options.Display = 'off'; % silence console output after fit completes
    fitval = lsqnonlin(fun,r0,lb,ub,options); % fit coefficients vector
    %options.Algorithm = 'levenberg-marquardt'; % trying different algorithm
    %fitval2 = lsqnonlin(fun,r0,lb,ub,options);
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
lifetime1 = 1/fitval(4); % ps
lifetime2 = 1/fitval(6); % ps
FWHM = 2*log(2)*sqrt(1/(2*fitval(2))); % ps

% Plot data & fit with residuals below
figure;
tiledlayout(2,2) % for multi-panel graph up to 2x2
nexttile([1 2]) % top panel
hold on
plot(t,signal,'o')
plot(t,fitcurve,'-')
hold off
legend('Data','Fit')
xlabel('Time (ps)')
ylabel('Signal (AU)')
nexttile([1 2]) % bottom panel
plot(t,resid)
legend('Residuals')
xlabel('Time (ps)')
ylabel('Residuals (AU)')

% Comparing trust-region-reflective (default algo) against levenberg-marq
%fitcurve2 = conv(exp(-fitval2(2)*(tconv-fitval2(1)).^2), ...
%    fitval2(3)*exp(-fitval2(4)*(tconv-fitval2(1)))+ ...
%    fitval2(5)*exp(-fitval2(6)*(tconv-fitval2(1)))+ ...
%    fitval2(7)*exp(-fitval2(8)*(tconv-fitval2(1))).* ...
%    cos(fitval2(9)*tconv-fitval2(10))+ ...
%    fitval2(11)*exp(-fitval2(12)*(tconv-fitval2(1))).* ...
%    cos(fitval2(13)*tconv-fitval2(14))+ ...
%    fitval2(15)*exp(-fitval2(16)*(tconv-fitval2(1))).* ...
%    cos(fitval2(17)*tconv-fitval2(18)));

%resid2 = signal-fitcurve2;

%figure;
%tiledlayout(2,2)
%nexttile([1 1]) % top-left panel
%hold on
%plot(t,signal,'o')
%plot(t,fitcurve,'-')
%hold off
%legend('Data','Fit')
%xlabel('Time (ps)')
%ylabel('Signal (AU)')
%nexttile([1 1]) % top-right panel
%hold on
%plot(t,signal,'o')
%plot(t,fitcurve2,'-')
%hold off
%legend('Data','Fit')
%xlabel('Time (ps)')
%ylabel('Signal (AU)')
%nexttile([1 1]) % bottom-left panel
%plot(t,resid)
%legend('Residuals')
%xlabel('Time (ps)')
%ylabel('Residuals (AU)')
%nexttile([1 1]) % bottom-right panel
%plot(t,resid2)
%legend('Residuals')
%xlabel('Time (ps)')
%ylabel('Residuals (AU)')
