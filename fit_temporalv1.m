%% Fit temporal sweep data as exponentially modified Gaussian
%%% (convolution of Gaussian and exponential decay)

% Load data
data = importdata('779m-706mW-1770.8nm-68mW-ND1-PMT1-10kHz-250kHz.txt');

% Parse data
x = data(:,1); % Save delay position as x
DC = data(:,2);
AC = data(:,3);

%%% Set cutoffs for baseline fitting
% Want to exclude points between xlow and xhigh
xlow = x(10); % xlow set as 10th data point
xhigh = x(end-6); % xhigh set as 7th-to-last data point

%%% Setup excluded points for baseline fitting
excludeDupes = zeros(size(x)); % using length(x) gives a matrix
for i = 1:length(x)
    if x(i) > xlow && x(i) < xhigh
        excludeDupes(i) = i;
    end
end

%%% Remove duplicates from exclude list
excludeunique = unique(excludeDupes(:).'); % adding the .' transposes it into a row vector
exclude1 = nonzeros(excludeunique); % isolate nonzero data positions to exclude

%%% Fit baselines with lines
f1 = fit(x,DC,'poly1','Exclude',exclude1);
f2 = fit(x,AC,'poly1','Exclude',exclude1);

%%% Baseline correction
% Save coefficient values from fit to vector
fitcoefDC = coeffvalues(f1);
fitcoefAC = coeffvalues(f2);

% Create vectors of fitted values
fitDC = polyval(fitcoefDC,x);
fitAC = polyval(fitcoefAC,x);

% Subtract fit vectors from raw data
corrDC = DC-fitDC;
corrAC = AC-fitAC;

% Specify fit options first
guesses = [0.04 6.4 228 0.063];
fo = fitoptions('Method','NonlinearLeastSquares', ...
    'StartPoint',guesses, ... % initial guesses
    'Lower',[0 0 -inf 0]); % set lower bounds
% Define fitting function - exponentially modified Gaussian
% A is amplitude, mx is the center (mean in x) in mm, 
% l is the exponential rate in mm-1, s is the Gaussian stdev in mm
posfit = fittype(['A*0.5*l*exp(0.5*l*(2*mx+l*s^2-2*x))*' ...
    '(1-erf((mx+l*s^2-x)/(s*sqrt(2))))'],'options',fo);

[curve,gof] = fit(x,corrAC,posfit);

%figure;
%plot(x,corrAC)
%hold on
%plot(curve)
%hold off
%xlabel('Position (mm)')
%ylabel('Corrected AC signal (AU)')

% Extract fitted t0 using value of mx
cmmps = 299792458*1E3*1E-12; % c ~ 0.3 mm/ps
clear t;
t(:) = (x(:)-curve.mx)*4/cmmps; % time vector in ps
t = t.'; % transpose to column vector

% Verifying
figure;
plot(t,corrAC)

% Write out corrected data as a txt file
writematrix([t corrAC],'779m_proc.txt')

% Looks like the time is pretty different because the exponentially
% modified Gaussian doesn't peak at t=0 (peaks slightly after)

% Extract lifetime using value of l
lifetime = 4/(cmmps*curve.l); % lifetime in ps


%% Testing log scale

% Seems like MATLAB won't plot on semilog axes unless your data
% vary over an order of magnitude?
%semilogy(x,AC)

% Taking log manually
lDC(:) = log(DC(:)); % 'log' does ln
lAC(:) = log(AC(:));

% Normalize data to 1
NDC(:) = DC(:)/max(DC);
NAC(:) = AC(:)/max(AC);
NCDC(:) = corrDC(:)/max(corrDC);
NCAC(:) = corrAC(:)/max(corrAC);
NlDC(:) = lDC(:)/max(lDC);
NlAC(:) = lAC(:)/max(lAC);

%figure;
%plot(x,NlAC)

NCAC(:) = abs(NCAC(:));
semilogy(t,NCAC)


%% Refitting to verify lifetime calculation
f2o = fitoptions('Method','NonlinearLeastSquares', ...
    'StartPoint',[0.33 0.7 0.4], ... % initial guesses - now only 3
    'Lower',[0 0 0]); % set lower bounds
% Now fitting in time with At, lt, and st
% Using mt as 'problem' variable to hold fixed in fitting at 0
tfit = fittype(['At*0.5*lt*exp(0.5*lt*(2*mt+lt*st^2-2*x))*' ...
    '(1-erf((mt+lt*st^2-x)/(st*sqrt(2))))'],'problem','mt','options',f2o);

[curve2,gof2] = fitval(t,corrAC,tfit,'problem',0); % holding mt=0

%figure;
%plot(t,corrAC)
%hold on
%plot(curve2)
%hold off
%xlabel('Time (ps)')
%ylabel('Corrected AC signal (AU)')

% Slight disagreement if mt floats
% Exact agreement if set mt=0 for tfit
lifetime2 = 1/curve2.lt; % lifetime in ps

%% Method from BonFIRE manuscript
% Fit 1.3 ps to 3.3 ps just with exponential - should be ~ 1.9 ps
% Couldn't image full temporal sweep because they did point scanning

% Make time vector where t0 is peak
[corrACpeak, t0index] = max(corrAC);
clear time
time(:) = (fitval(:)-fitval(t0index))*4/cmmps; % time vector in ps
time = time.';

% Points to include are indices 23,24,25,26
incl = [23 24 25 26];
excl = [linspace(1,min(incl)-1,min(incl)-1) ...
    linspace(max(incl)+1,length(fitval),length(fitval)-max(incl))].';

% Fitting with single exponential function
BFfit = fitval(time,corrAC,'exp1','Exclude',excl);

%figure;
%plot(BFfit,time,corrAC,excl)
%xlabel('Time (ps)')
%ylabel('Corrected AC signal (AU)')
%ylim([-inf 0.25]) % -inf lets MATLAB auto-set the lower limit

lifetime3 = -1/BFfit.b;



%% Simulating fit
x1 = linspace(-1000,1000,1001);
a = 0.04;
b = 10;
c = 1;
d = 0.005;

gaussfn = @(x2) a*exp(-(x2-b).^2/(2*c^2));
expfn = @(x2) d*exp(-d*x2);

% Convolve functions
gexp = @(x1) integral(@(x2) gaussfn(x2).*expfn(x1-x2),-inf,inf);
expg = @(x1) integral(@(x2) expfn(x2).*gaussfn(x1-x2),-inf,inf);
% gexp works better

value = zeros([1,length(x1)]);
value2 = zeros([1,length(x1)]);

for i=1:length(x1)
    value(i) = gexp(x1(i));
    value2(i) = expg(x1(i));
end

figure;
plot(x1,value)
hold on
plot(x1,value2)
hold off