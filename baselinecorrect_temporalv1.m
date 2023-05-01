%% Extract peak height from temporal sweep data
% Load data
data = importdata('779m-706mW-1770.8nm-68mW-ND1-PMT1-10kHz-250kHz.txt');

% Parse data
% Pull columns
x = data(:,1); % Save delay position as x
DC = data(:,2); % Save channel 1 as DC
AC = data(:,3); % Save channel 2 as AC

%%% Set cutoffs for baseline fitting
% Want to exclude points between xlow and xhigh
% xlow set as 10th data point
xlow = x(10);
% xhigh set as 7th-to-last data point
xhigh = x(end-6);

%%% Setup excluded points for baseline fitting

excludeDupes = zeros(size(x)); % using length(x) gives a matrix
for i = 1:length(x)
    if x(i) > xlow && x(i) < xhigh
        excludeDupes(i) = i;
    end
end

%%% Remove duplicates from exclude list
% Adding the .' transposes it into a row vector
excludeunique = unique(excludeDupes(:).');
% Isolate nonzero data positions to exclude
exclude1 = nonzeros(excludeunique);

%%% Fit baselines with lines
% 'Exclude' excludes the data points specified by exclude1
% by position in the vector, not the value of the number
f1 = fit(x,DC,'poly1','Exclude',exclude1);
f2 = fit(x,AC,'poly1','Exclude',exclude1);

% Plot signal and baseline fit with excluded points
figure;
plot(f2,x,AC,exclude1)
xlabel('Position (mm)')
ylabel('AC signal')

%%% Baseline correction
% Save coefficient values from fit to vector
fitcoefDC = coeffvalues(f1);
fitcoefAC = coeffvalues(f2);

% Check formula of fit
%formula(f1)
%formula(f2)

% Create vectors of fitted values
fitDC = polyval(fitcoefDC,x);
fitAC = polyval(fitcoefAC,x);

% Subtract fit vectors from raw data
corrDC = DC-fitDC;
corrAC = AC-fitAC;

%figure;
%plot(x,corrAC)

%%% Write data to txt file
%T = table(x,DC,AC,corrDC,corrAC);
%T1 = T;
%T = table(T1,x,DC,AC,corrDC,corrAC);
%Tfinal = splitvars(T);
%writetable(Tfinal,'processed.txt')


%% Appendix
% Brute-force baseline correction with for loop
% Better to just use polyval
%corrDC = zeros(size(x));
%corrAC = zeros(size(x));
%for i = 1:length(x)
%    pos = x(i);
%    corrDC(i) = DC(i) - f1(pos);
%    corrAC(i) = AC(i) - f2(pos);
%end

% Kwan's way to background-correct
%i_low = 10;
%i_high = length(x)-6;

%x_bs = x([1:i_low i_high:end]);
%DC_bs = DC([]);