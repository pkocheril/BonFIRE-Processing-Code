%%% Loop through folders for .txt files, baseline correct, and 
%%% combine into a single txt file
%%% Version History
%%% v2 - filename parsing for power-normalizing, picking maxima
%%% v3 - finding time 0, image work
%%% v4 - convolution temporal fitting, normalizing with ND -- in progress
%%% -- replaced by batchfit_temporal

%%% Copy to directory with data in subfolders
%%% Update guesses as needed - can try fit_temporal.m on one first

%% Batch processing
%%% Baseline correction, power-normalization, temporal sweep fitting,
%%% and writing out processed data in tables
guesses = [0.04 6.4 228 0.063];
Table1 = []; % for writing to individual _proc.txt files
Table2 = []; % for combining all data into one table for further analysis
D = pwd; % get current directory
S = dir(fullfile(D,'*')); % search current directory
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D 
% - done by taking all directories and removing . and ..
for ii = 1:numel(N) % ii = subfolder number
    T = dir(fullfile(D,N{ii},'*.txt')); % look for .txt files in subfolders
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
%%% Comment here to stop figure saving
        saveas(fig1,DCcheckname)

        fig2 = figure('visible','off');
        plot(f2,x,AC,exclude1)
        xlabel('Position (mm)')
        ylabel('AC signal (AU)')
        ACcheckname = string(F) + '_ch2_bs_corr.png';
%%% Comment here to stop figure saving
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

        % Parse filename - could also use fileparts for this
        folderinfo = split(F,"/"); % split path into folders
        fileinfo = split(folderinfo(end),"-"); % split filename by hyphens
        subfolderinfo = split(N{ii}," "); % split subfolder by spaces

        % Extract experimental parameters from filename of type:
        % [pumpWL]-[pumppower]-[signalWL]-[IRpower]-[ND]-[PMTgain]-[chopper]-[PMTBW]-[etc]
        pumpWLstrfull = char(fileinfo(1)); % extract pump WL from file name - "760m"
        pumpWLstr = pumpWLstrfull(1:end-1); % remove "m" suffix - "760"
        pumpWL = sscanf(pumpWLstr,'%f'); % convert to double-precision number
        pumppowerstrfull = char(fileinfo(2)); % extract pump power from file name
        pumppowerstr = pumppowerstrfull(1:end-2); % remove "mW"
        pumppower = sscanf(pumppowerstr,'%f');
        signalWLstrfull = char(fileinfo(3)); % extract signal WL from file name
        signalWLstr = signalWLstrfull(1:end-2); % remove "nm"
        signalWL = sscanf(signalWLstr,'%f');
        IRpowerstrfull = char(fileinfo(4)); % extract IR power from file name
        IRpowerstr = IRpowerstrfull(1:end-2); % remove "mW"
        IRpower = sscanf(IRpowerstr,'%f');
        DFGWL = 1/((2/signalWL) - (1/1031.2)); % calculate DFG IR WL (nm)
        DFGWN = 1E7/DFGWL; % calculate DFG IR WN (cm-1)

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
            bookkeeping(7) = 'DFGWN = ' + string(DFGWN) + ' cm-1';
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
            values(7) = DFGWN;
        end

        % Save individual data to a table
        Table1 = table(bookkeeping,values,x,DC,AC,fitDC,fitAC,corrDC,corrAC,corrnormDC,corrnormAC);
        % Write to individual txt files in the individual folders
        % Includes raw data and fits to check back
        procfile = string(F) + '_proc.dat'; % writing to .dat so it doesn't get picked up if I run the analysis again
%%% Comment here to stop individual data saving
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

        %%% Combine processed data into a table to be exported
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
writetable(Table2,'batch_processed_v4.txt');

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

% Create array with only numeric data from Table2
dataarray = zeros(height(Table2),width(Table2));
dataarray(:,:) = table2array(Table2(:,:));

%x = dataarray(:,3); % in case I load in Table2 from batch_processed data

%%% File annotation
%	1	;	
%	6	;	
%	11	;	
%	16	;	
%	21	;	
%	26	;	
%	31	;	
%	36	;	
%	41	;	
%	46	;	
%	51	;	
%	56	;	
%	61	;	
%	66	;	
%	71	;	
%	76	;	
%	81	;	
%	86	;	
%	91	;	
%	96	;	
%	101	;	
%	106	;	
%	111	;	
%	116	;	
%	121	;	
%	126	;	
%	131	;	
%	136	;	
%	141	;	
%	146	;	
%	151	;	

% Example parsing
%pmtgain = [1 2 4 6 8 10 10 20 40 60 80 100];
%pmtgainlabel(:) = string("Gain " + pmtgain(:));
%pmtarrayAC(:,1) = dataarray(:,190); % PMT1
%pmtarrayAC(:,2) = dataarray(:,210); % PMT2

% Example finding time 0
%cmmps = 299792458*1E3*1E-12; % c in mm/ps
%[pmtarray1peak, t0index] = max(pmtarray(:,1));
%t(:) = (x(:)-x(t0index))*4/cmmps; % time vector in ps

% Example plot
%figure;
%plot(x,pmtarrayAC(:,1),x,pmtarrayAC(:,2),'-o') % the "-" means you get the trace and markers together
%xlabel('Position (mm)')
%ylabel('Corrected AC signal (AU)')
%legend(pmtgainlabel)

% Example peak vs integration data
%pmtmax = max(pmtarrayAC);
%pmtint = sum(pmtarrayAC);

%figure;
%plot(pmtgain,pmtint,'-o') % the "-" means you get the trace and markers together
%xlabel('PMT Gain')
%ylabel('Integrated AC signal (AU)')

% Example 3D-plot
%[gain,pos] = meshgrid(pmtgain,x);
%figure;
%surf(gain,pos,pmtarrayAC)
%xlabel('PMT Gain')
%ylabel('Position (mm)')
%zlabel('Corrected AC signal (AU)')

% Example error bar plot
%figure;
%errorbar(nbconc,meanints,stdevints)
%xlabel('[NB] (mM)')
%ylabel('Mean AC integrated signal (AU)')

%% Processing images
% Example loading measurements from Fiji
%unsortimage = readtable('./Images/BonFIRE/BFimages.csv','VariableNamingRule','preserve');
%imagetable = sortrows(unsortimage,2); % sort by second column (should be filename)

% Isolate RawIntDen (sum of pixel values) in a vector
%imagedata(:,1) = imagetable(:,8); % still a table at this point
%rawintden(:,:) = table2array(imagedata(:,:));

% Separate rawintden into individual columns in an array
%BFimages = zeros(3,9);
%BFimages(:,1) = rawintden(1:3); % 0.5 mM
%BFimages(:,9) = rawintden(4:6); % 10 mM
%BFimages(:,2) = rawintden(7:9); % 1 mM
%BFimages(:,3) = rawintden(10:12); % 3 mM
%BFimages(:,4) = rawintden(13:15); % 5 mM
%BFimages(:,5) = rawintden(16:18); % 6 mM
%BFimages(:,6) = rawintden(19:21); % 7 mM
%BFimages(:,7) = rawintden(22:24); % 8 mM
%BFimages(:,8) = rawintden(25:27); % 9 mM


%% Appendix
%%% Load-in tif images - needed the image processing toolbox (will probably
%%% work now)
%cd Images
%D2 = pwd; % set directory
%S2 = dir(fullfile(D2,'*tempOFF*.tif')); % search D2 for tempOFF.tif files
%N2 = {S2.name}; % list of .tif files in D2
%tempoff = [];
%for ii = 1:numel(N2)
%    image = Tiff(char(N2(ii)),'r');
%    tiffdata = read(image);
%    tempoff = [tempoff tiffdata];
%end

%S3 = dir(fullfile(D2,'*tempON*.tif')); % search D2 for tempON.tif files
%N3 = {S3.name}; % list of .tif files in D2
%tempon = [];
%for ii = 1:numel(N3)
%    image = Tiff(char(N3(ii)),'r');
%    tiffdata = read(image);
%    tempoff = [tempoff tiffdata];
%end

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