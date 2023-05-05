%%% Loop through folders for .txt files, baseline correct, and 
%%% combine into a single txt file
%%% Version History
%%% v2 - filename parsing for power-normalizing, picking maxima

%%% Open Terminal before running
%%% Copy to directory with data in subfolders

%% Baseline correction, power-normalization, and combining into a single Table
Table1 = []; % for writing to individual _proc.txt files
Table2 = []; % for combining all data into one table for further analysis
D = pwd; % get current directory
S = dir(fullfile(D,'*')); % search current directory
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D
for ii = 1:numel(N) % ii = subfolder number
    T = dir(fullfile(D,N{ii},'*.txt')); % look for .txt files in subfolder
    C = {T(~[T.isdir]).name}; % all .txt files in subfolder ii
    for jj = 1:numel(C) % jj = file number in subfolder ii
        F = fullfile(D,N{ii},C{jj});
        %%% Process data in file F
        data = importdata(F); % load data
        x = data(:,1); % Save delay position as x
        DC = data(:,2); % Save channel 1 as DC
        AC = data(:,3); % Save channel 2 as AC

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

        %%% Save figures invisibly to be able to easily check baseline
        %%% correction after the fact
        fig1 = figure('visible','off');
        plot(f1,x,DC,exclude1)
        xlabel('Position (mm)')
        ylabel('DC signal (AU)')
        DCcheckname = string(F) + '_ch1_bs_corr.png';
        saveas(fig1,DCcheckname)

        fig2 = figure('visible','off');
        plot(f2,x,AC,exclude1)
        xlabel('Position (mm)')
        ylabel('AC signal (AU)')
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

        % Parse filename
        folderinfo = split(F,"/"); % could also use fileparts for this
        fileinfo = split(folderinfo(end),"-");
        subfolderinfo = split(N{ii}," ");

        % Extract experimental parameters from filename of type:
        % [pumpWL]-[pumppower]-[signalWL]-[IRpower]-[ND]-[PMT]-[chopper]-[etc]
        pumpWLstrfull = char(fileinfo(1)); % extract pump WL from file name - "760nm" or "760m"
        if strlength(pumpWLstrfull) == 5
            pumpWLstr = pumpWLstrfull(1:end-2); % remove "nm" suffix - "760"
            pumpWL = sscanf(pumpWLstr,'%f'); % convert to double-precision number
        else
            pumpWLstr = pumpWLstrfull(1:end-1); % remove "m" suffix - "760"
            pumpWL = sscanf(pumpWLstr,'%f'); % convert to double-precision number
        end
        pumppowerstrfull = char(fileinfo(2)); % extract pump power from file name
        pumppowerstr = pumppowerstrfull(1:end-2); % remove "mW"
        pumppower = sscanf(pumppowerstr,'%f');
        signalWLstrfull = char(fileinfo(3)); % extract signal WL from file name
        signalWLstr = signalWLstrfull(1:end-2); % remove "nm"
        signalWL = sscanf(signalWLstr,'%f');
        IRpowerstrfull = char(fileinfo(4)); % extract IR power from file name
        IRpowerstr = IRpowerstrfull(1:end-2); % remove "mW"
        IRpower = sscanf(IRpowerstr,'%f');
        IRWL = 1/((2/signalWL) - (1/1031.2)); % calculate DFG IR WL (nm)
        IRWN = 1E7/IRWL; % calculate DFG IR WN (cm-1)

        % Normalize to 800 mW pump power and 70 mW IR power
        pumpnorm = 800/pumppower; %normalization factor for pump
        IRnorm = 70/IRpower; %normalization factor for IR
        corrnormDC = corrDC*pumpnorm*IRnorm;
        corrnormAC = corrAC*pumpnorm*IRnorm;

        % Create bookkeeping string array (only works if more data points
        % than labels desired)
        if length(x) >= 30
            bookkeeping = strings(length(x),1);
            bookkeeping(1) = 'subfolder = ' + string(N{ii});
            bookkeeping(2) = 'file = ' + string(folderinfo(end));
            bookkeeping(3) = 'pumpWL = ' + string(pumpWL) + ' nm';
            bookkeeping(4) = 'pumppower = ' + string(pumppower) + ' mW';
            bookkeeping(5) = 'signalWL = ' + string(signalWL) + ' nm';
            bookkeeping(6) = 'IRpower = ' + string(IRpower) + ' mW';
            bookkeeping(7) = 'IRWN = ' + string(IRWN) + ' cm-1';
            for i = 5:length(fileinfo)
                bookkeeping(i+3) = string(fileinfo(i));
            end
            bkl = length(fileinfo);
            for i = 1:length(subfolderinfo)
                bookkeeping(i+3+bkl) = string(subfolderinfo(i));
            end

        end

        if length(x) >= 7
            values = zeros(length(x),1);
            values(3) = pumpWL;
            values(4) = pumppower;
            values(5) = signalWL;
            values(6) = IRpower;
            values(7) = IRWN;
        end

        % Save individual data to a table
        Table1 = table(bookkeeping,values,x,DC,AC,fitDC,fitAC,corrDC,corrAC,corrnormDC,corrnormAC);
        % Write to individual txt files in the individual folders
        % Includes raw data and fits to check back
        procfile = string(F) + '_proc.dat'; % writing to .dat so it doesn't get picked up if I run the analysis again
        writetable(Table1,procfile);

        %%% Split into labels and values for easier MATLAB processing
        % Setup file names
        filename = char(folderinfo(end));
        tablename = filename(end-40:end-4); % need char to be able to trim by position
        bookname = "Info_" + string(ii) + "_" + string(jj) + "_" + string(tablename); % need quotes to make a string
        valuename = "Value_" + string(ii) + "_" + string(jj) + "_" + string(tablename);
        xname = "x_" + string(ii) + "_" + string(jj) + "_" + string(tablename); % need strings for table variable names
        DCname = "DC_" + string(ii) + "_" + string(jj) + "_" + string(tablename);
        ACname = "AC_" + string(ii) + "_" + string(jj) + "_" + string(tablename);
        oldnames = ["bookkeeping","values","x","corrnormDC","corrnormAC"]; 
        newnames = [bookname, valuename, xname, DCname, ACname];

        %%% Combine all data into a table to be exported
        TF = isempty(Table2);
        if TF == 1 % for the first dataset
            Table2 = strings(length(x),1); % create table of correct length
            Table2 = table(Table2,bookkeeping,values,x,corrnormDC,corrnormAC);
            Table2(:,1) = []; % remove empty first column from making table
            % Rename variables to be unique by filename
            Table2 = renamevars(Table2,oldnames,newnames);
        else % for successive datasets
            %Table2 = table(Table2,bookkeeping,x,corrnormDC,corrnormAC);
            datatable = table(bookkeeping,values,x,corrnormDC,corrnormAC); % put current data in a new table
            Table2 = [Table2 datatable]; % combine previous data with current data
            Table2 = renamevars(Table2,oldnames,newnames);
        end
    end
end

%%% Write out all data to one txt file
%writetable(Table2,'batch_processed.txt');

%% Load batch-processed data

%%% Figure out for if I need to clear the workspace and don't want to wait
%%% for the analysis to run again
%Table2 = readtable('batch_processed.txt','VariableNamingRule','preserve');
% Close but not quite right


%% Case-by-case analysis
% Create infotable to figure out what files went where
infowidth = width(Table2)/5;
for i = 1:infowidth
    j = 5*i-4;
    infotable(:,i) = Table2(:,j);
end

% Create array with only numbers from Table2
dataarray = zeros(height(Table2),width(Table2));
dataarray(:,:) = table2array(Table2(:,:));

x = dataarray(:,3); % in case I load in Table2 from batch_processed data

%%% File annotation - PMT10 unless otherwise specified
% 1; 1. Rh800 100 uM alignment benchmark
% 6; 10. 6 mM NB 5 uM Rh800
% 11; 11. 5 mM NB 5 uM Rh800
% 16; 12. 9 mM NB control
% 21; 13. 8 mM NB control
% 26; 14. 7 mM NB control
% 31; 15. 6 mM NB control
% 36; 16. 5 mM NB control
% 41; 16. 5 mM NB control repeat
% 46; 17. 3 mM NB control
% 51; 18. 1 mM NB control
% 56; 19. 0.1 mM NB control
% 61; 19. 0.1 mM NB control PMT100
% 66; 2. 5 uM Rh800 BonFIRE control
% 71; 20. 3 mM NB 5 uM Rh800
% 76; 21. 1 mM NB 5 uM Rh800
% 81; 21. 1 mM NB 5 uM Rh800 dt 1s tc 300ms
% 86; 21. 1 mM NB 5 uM Rh800 dt 1s tc 50 ms
% 91; 22. 0.1 mM NB 5 uM Rh800
% 96; 22. 0.1 mM NB 5 uM Rh800 PMT100
% 101; 3. 5 uM Rh800 vsFRET control (Rh800 filter)
% 106; 3. 5 uM Rh800 vsFRET control (vsFRET filter)
% 111; 4. 10 mM NB BonFIRE control PMT1
% 116; 4. 10 mM NB BonFIRE control
% 121; 5. 10 mM NB vsFRET control shutteroff
% 126; 5. 10 mM NB vsFRET control 
% 131; 5. 10 mM NB vsFRET control shutteroff IR on
% 136; 5. 10 mM NB vsFRET control shutteroff PMT100 - comparing to previous
% 141; 5. 10 mM NB vsFRET control shutteroff PMT1 - comparing to previous
% 146; 6. 10 mM NB 5 uM Rh800 PMT1
% 151; 6. 10 mM NB 5 uM Rh800
% 156; 7. 9 mM NB 5 uM Rh800 PMT1
% 161; 7. 9 mM NB 5 uM Rh800 
% 166; 7. 9 mM NB 5 uM Rh800 PMT2
% 171; 7. 9 mM NB 5 uM Rh800 PMT4
% 176; 7. 9 mM NB 5 uM Rh800 PMT6
% 181; 7. 9 mM NB 5 uM Rh800 PMT8
% 186; 8. 8 mM NB 5 uM Rh800 PMT1
% 191; 8. 8 mM NB 5 uM Rh800 PMT10 (initial)
% 196; 8. 8 mM NB 5 uM Rh800 PMT100
% 201; 8. 8 mM NB 5 uM Rh800 PMT10 (final)
% 206; 8. 8 mM NB 5 uM Rh800 PMT2
% 211; 8. 8 mM NB 5 uM Rh800 PMT20
% 216; 8. 8 mM NB 5 uM Rh800 PMT4
% 221; 8. 8 mM NB 5 uM Rh800 PMT40
% 226; 8. 8 mM NB 5 uM Rh800 PMT6
% 231; 8. 8 mM NB 5 uM Rh800 PMT60
% 236; 8. 8 mM NB 5 uM Rh800 PMT8
% 241; 8. 8 mM NB 5 uM Rh800 PMT80
% 246; 9. 7 mM NB 5 uM Rh800

% Looking at PMT gain trend with 8 mM NB 5 uM Rh800 data
pmtgain = [1 2 4 6 8 10 10 20 40 60 80 100];
pmtgainlabel(:) = string("Gain " + pmtgain(:));
pmtarrayAC(:,1) = dataarray(:,190); % PMT1
pmtarrayAC(:,2) = dataarray(:,210); % PMT2
pmtarrayAC(:,3) = dataarray(:,220); % PMT4
pmtarrayAC(:,4) = dataarray(:,230); % PMT6
pmtarrayAC(:,5) = dataarray(:,240); % PMT8
pmtarrayAC(:,6) = dataarray(:,195); % PMT10 initial
pmtarrayAC(:,7) = dataarray(:,205); % PMT10 final
pmtarrayAC(:,8) = dataarray(:,215); % PMT20
pmtarrayAC(:,9) = dataarray(:,225); % PMT40
pmtarrayAC(:,10) = dataarray(:,235); % PMT60
pmtarrayAC(:,11) = dataarray(:,245); % PMT80
pmtarrayAC(:,12) = dataarray(:,200); % PMT100

% Plot PMT gain
%figure;
%plot(x,pmtarrayAC(:,1),x,pmtarrayAC(:,2),x,pmtarrayAC(:,3),x,pmtarrayAC(:,4),x,pmtarrayAC(:,5),x,pmtarrayAC(:,6),x,pmtarrayAC(:,7),x,pmtarrayAC(:,8),x,pmtarrayAC(:,9),x,pmtarrayAC(:,10),x,pmtarrayAC(:,11),x,pmtarrayAC(:,12))
%xlabel('Position (mm)')
%ylabel('Corrected AC signal (AU)')
%legend(pmtgainlabel)

% Integrated PMT gain data
pmtmax = max(pmtarrayAC);
pmtint = sum(pmtarrayAC);

%figure;
%plot(pmtgain,pmtint,'-o') % the "-" means you get the trace and markers together
%xlabel('PMT Gain')
%ylabel('Integrated AC signal (AU)')

% 3D-plot of PMT gain as a function of position and PMT gain
[gain,pos] = meshgrid(pmtgain,x);

%figure;
%surf(gain,pos,pmtarrayAC)
%xlabel('PMT Gain')
%ylabel('Position (mm)')
%zlabel('Corrected AC signal (AU)')

% Isolate [NB] controls
nbconc = [0.1 1 3 5 6 7 8 9 10]; % mM
nbconclabel(:) = string(nbconc(:) + " mM");
nbcontAC(:,1) = dataarray(:,60); % 0.1 mM
nbcontAC(:,2) = dataarray(:,55); % 1 mM
nbcontAC(:,3) = dataarray(:,50); % 3 mM
nbcontAC(:,4) = dataarray(:,40); % 5 mM
nbcontAC(:,5) = dataarray(:,35); % 6 mM
nbcontAC(:,6) = dataarray(:,30); % 7 mM
nbcontAC(:,7) = dataarray(:,25); % 8 mM
nbcontAC(:,8) = dataarray(:,20); % 9 mM
nbcontAC(:,9) = dataarray(:,130); % 10 mM

% Plot NB controls
%figure;
%plot(x,nbcontAC(:,1),x,nbcontAC(:,2),x,nbcontAC(:,3),x,nbcontAC(:,4),x,nbcontAC(:,5),x,nbcontAC(:,6),x,nbcontAC(:,7),x,nbcontAC(:,8),x,nbcontAC(:,9))
%xlabel('Position (mm)')
%ylabel('Corrected AC signal (AU)')
%legend(nbconclabel)

% Max and integrated signal as a function of [NB]
nbmax = max(nbcontAC);
nbint = sum(nbcontAC);
%figure;
%plot(nbconc,nbmax)
%xlabel('[NB] (mM)')
%ylabel('Max AC signal')

%figure;
%plot(nbconc,nbint,'-o') % the "-" means you get the trace and markers together
%xlabel('[NB] (mM)')
%ylabel('Integrated AC signal (AU)')

% Looks like integrated and max signal agree well for single-decay temporal sweeps


%% Appendix
%%% Write out all data to one txt file (replaced)
%Tfinal = Table2;
%Tfinal(:,1) = [];
%%% Note that writetable will overwrite an existing file
%writetable(Tfinal,'batch_processed.txt');

%%% Combine individual files (replaced)
% Procfiles = dir(fullfile(D,'*_proc.txt')); % get all _proc.txt files in D
% % Make Tcomb as an empty string array with same length as x
% Tcomb = strings(length(x),1);
% for i = 1:length(Procfiles)
%     procdata = importdata(Procfiles(i,1).name,',');
%     Tcomb = [Tcomb procdata.data]; % right-appends data
% end

% Import all _proc.txt files written from previous (replaced)
%clear
%D1 = pwd;
%S1 = dir(fullfile(D1,'*'));
%N1 = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D
%for ii = 1:numel(N1) % ii = subfolder number
%    T1 = dir(fullfile(D1,N1{ii},'*_proc.txt')); % look for _proc.txt files in subfolder
%    C1 = {T1(~[T1.isdir]).name}; % all .txt files in subfolder ii
%    for jj = 1:numel(C1) % jj = file number in subfolder ii
%        F1 = fullfile(D1,N1{ii},C1{jj});
%    end
%end

%%% From when I had to split concatenated tables (replaced)
%Table2 = splitvars(Table2); % split concatenated Table2 into single-columns