%% Loop through folders for .txt files, baseline correct, and combine
%%% into a single txt file
%%% Open Terminal before running
%%% Copy to directory with data in subfolders
Table1 = 0;
Table2 = [];
D = pwd;
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.txt')); % look for .txt files
    C = {T(~[T.isdir]).name}; % files in subfolder
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        %%% Process data in file F
        % Load data
        data = importdata(F);
        
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
        
        %%% Save figures invisibly to be able to check baseline
        %%% correction after the fact and easily
        fig1 = figure('visible','off');
        plot(f1,x,DC,exclude1)
        xlabel('Position (mm)')
        ylabel('DC signal')
        DCcheckname = string(F) + '_ch1_bs_corr.png';
        saveas(fig1,DCcheckname)
        fig2 = figure('visible','off');
        plot(f2,x,AC,exclude1)
        xlabel('Position (mm)')
        ylabel('AC signal')
        ACcheckname = string(F) + '_ch2_bs_corr.png';
        saveas(fig2,ACcheckname)
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

        % Create bookkeeping string array
        if length(x) >= 3
            bookkeeping = strings(length(x),1);
            bookkeeping(1) = 'folder = ' + string(D);
            bookkeeping(2) = 'subfolder = ' + string(N{ii});
            bookkeeping(3) = 'file = ' + string(F);
        end
        
        % Write to individual txt files in the individual folders
        % Includes raw data and fits to check back
        FinalFile = string(F) + '_proc.txt';
        Table1 = table(bookkeeping,x,DC,AC,fitDC,fitAC,corrDC,corrAC);
        writetable(Table1,FinalFile);

        %%% Write all data to one txt file     
        TF = isempty(Table2);
        if TF == 1
            %c = 'hello world' % debugging
            Table2 = strings(length(x),1);
            Table2 = table(Table2,bookkeeping,x,corrDC,corrAC);
        else
            %c = 'goodbye world' % debugging
            Table2 = table(Table2,bookkeeping,x,corrDC,corrAC);
        end
        Table2 = splitvars(Table2); % split concatenated Table2 into single-columns
    end
end
% Write final table
Tfinal = Table2;
Tfinal(:,1) = [];
%%% Note that writetable will overwrite an existing file
writetable(Table2,'batch_processed.txt');



%% Appendix
%%% Combine individual files (unfinished)
% Procfiles = dir(fullfile(D,'*_proc.txt')); % get all _proc.txt files in D
% % Make Tcomb as an empty string array with same length as x
% Tcomb = strings(length(x),1);
% for i = 1:length(Procfiles)
%     procdata = importdata(Procfiles(i,1).name,',');
%     Tcomb = [Tcomb procdata.data]; % right-appends data
% end
% 
% % Get rid of first empty column from Tcomb
% Tfinal = Tcomb;
% Tfinal(:,1) = [];
% writematrix(Tfinal,'batch_processed.txt')