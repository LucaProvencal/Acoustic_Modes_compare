function [PeakFreqComp, stat, Fabove, Fbelow] = acoustic_modes_compare(FRF, freqStart, freqStop, peakfindTol, tolInd, xlsExport)
%% Acoustic Modes Comparison Script
%  August 3, 2016
%
% ------------------------------------------------------------------------
% Authors - Luca Provencal (provencal.l@northeastern.edu)
%           Chris Beale (christopher_beale@student.uml.edu)
% ------------------------------------------------------------------------
% 
% This function's purpose is to check an acoustic cavity's day-to-day test 
% results for consistency. It does so by analyzing and processing each 
% day's frequency response function and comparing it with the other days' 
% functions.
%
% ------------------------------------------------------------------------
% I/O
% ------------------------------------------------------------------------
%   Inputs: 
%   
%       FRF - The test data that is to be analyzed. For 
%           acoustic_modes_compare to work, each day's FRF data MUST be 
%           stored in a structure in the following format:
%
%                *** "Struct_Name.Trace_Name.Microphone_Name" ***
% 
%           The fieldnames do not matter, as they are pulled from the 
%           structure using MATLAB's "fieldnames()" function. As long as 
%           the data is in the above format, the function will properly 
%           analyze it. This input is required.
%           
%       freqStart - Lower limit of the frequencies acoustic_modes_compare 
%           will process. Should typically be the lowest frequency that is 
%           excited during acoustic testing. This input is required.
%
%       freqStop - Upper limit of the frequencies acoustic_modes_compare 
%           will process. Similarly, freqStop should be the highest 
%           frequency excited during acoustic testing. This input is 
%           required.
%
%       peakfindTol - The tolerance used to calculate the peaks of each 
%           FRF. Acoustic_modes_compare uses "peakfinder.m" to find these
%           peaks. Essentially, the higher the value for peakfindTol, the
%           smaller the number of peaks found will be. For more info, see
%           "peakfinder.m" (peakfinder calls this value "sel"). This input
%           is not required. It defaults to .5.
%
%       tolInd - Indicates which method of tolerancing is desired: a single
%           tolerance or a tolerance that increases at higher frequencies.
%           A 0 indicates a static tolerance while a 1 indicates a
%           dynamic tolerance. Exactly how the tolerance is increased can 
%           be edited starting at line 179. By default, however, the 
%           tolerance is increased every 5000 Hz. This input is not 
%           required, and has a default value of 0 (static tolerance).
%          
%
%       xlsExport - Indicates whether or not to export the output data to
%           an excel sheet or sheets. A 0 indicates that an excel sheet(s) 
%           will not be generated while a 1 indicates that an excel 
%           sheet(s) will be generated. This input is also not 
%           required, and has a default value of 0 (no export).
%
%   Outputs:
%       
%       PeakFreqComp - Peak-frequency comparison table. The first few 
%           columns will contain frequency values organized so that at a
%           glance, the user can tell how well each day's test matches with
%           the other days. The number of columns in this comparison
%           section is equivalent to the number of tests being compared.
%           This section of PeakFreqComp can be interpreted as follows:
%                                                    
%           T1 Freq (Hz) | T2 Freq (Hz)     A zero acts as a placeholder, 
%           _____________|_____________     indicating the test whose 
%                10      |      0        <- column it is in does not have
%           _____________|_____________     frequencies matching with those 
%                0       |      15          on the same line. 
%           _____________|_____________     
%                30.1    |      30       <- A line without a zero cell 
%           _____________|_____________     indicates all peaks on the line
%                                           are a match.
%
%           PeakFreqComp contains other data columns placed after the 
%           comparison section. The order is as follows:
%           
%           - Boolean: Contains a 1 if two or more frequencies
%           match, and a 0 if there is no match on the line. 
%           
%           - Average: Contains the average frequency value of each
%           row of the comparison section.
%
%           - Standard Deviation: Contains the standard deviation of 
%           each row of the comparison section.
%
%           - Delta F above: Contains the difference of the average of the 
%           current row and the above row.
%
%           - Delta F below: Contains the difference of the average of the 
%           below row and the current row.
%
%           - Trace match column: Shows which traces match on each line. 
%           Format is a string separated by commas.
% 
%       stat - Table of statistics. For each microphone, it displays number 
%           of matched sets, number of non-matched peaks, average delta 
%           value, and the trace with the most matches.
% 
% ------------------------------------------------------------------------
% *** Peak Comparison requires the "peakfinder.m" function to be located  
%     in the MATLAB path. ***
% ------------------------------------------------------------------------
%
% ------------------------------ END INTRO -------------------------------

%% Configure and Validate Inputs
if nargin < 3
   error('Not enough input arguments. Please specify input data as well as start and stop frequency.')
end

if nargin < 4 || isempty(peakfindTol) % executes if peakfindTol is not specified
    peakfindTol = .5; % sets default value of peakfindTol to .5
end

if nargin < 5 || isempty(tolInd) % executes if tolInd is not specified
    tolInd = 0; % sets default method of tolerancing to static
else
    if tolInd ~= 0 && tolInd ~= 1 % returns an error if input value is not 1 or 0
        error('acoustic_modes_compare:InvtolInd','Invalid input. Tolerance indicator must be either 1 for dynamic tolerance or 0 for static tolerance.');
    end
end

if nargin < 6 || isempty(xlsExport) % executes if xlsExport is not specified
    xlsExport = 0; % tells function to not export to Excel by default
else
    if xlsExport ~= 0 && xlsExport ~= 1 % returns an error if input value is not 1 or 0
        error('acoustic_modes_compare:InvxlsExport','Invalid input. xlsExport must be either 1 to export or 0 to not export.');
    end
end

%% Find fieldnames and dF,d and implement Fstart and Fstop 

tNames = fieldnames(FRF); % grabs trace names from input structure                           
mNames = fieldnames(FRF.(tNames{1})); % grabs microphone names from input structure

Nt = length(tNames); % Stores length of tNames in variable Nt for simplicity.

F = FRF.Trace1.(mNames{1})(:, 1); 
dF = mean(diff(F)); % gets dF value to be used in calculating tolerances
FStart = find(F < freqStart + dF & F > freqStart - dF); % finds the index of the input FRF's where the start and stop frequencies occur. 
FStop = find(F < freqStop + dF & F > freqStop - dF);

%% Find peaks and arrange into structure
% Finds peaks for every mic for every trace. Creates a new structure called "peakStruc" to store frequency and amplitude values for peaks only.
% Strucure is in same format as the input "FRF" structure.

peakStruc = [];

for i = 1:Nt

    for jj = 1:numel(mNames)                    
        Trace = FRF.(tNames{i}).(mNames{jj});
        Trace = Trace(FStart:FStop, 1:2);
        [PeakLoc] = peakfinder(Trace(1:end, 2), peakfindTol, -inf, 1, false, false);
        peakStruc.(tNames{i}).(mNames{jj}) = Trace(PeakLoc, 1:2);

    end
    
end

%% Rearrange peakStruc into struct with processable matrices
% Reformats the peakStruc from "Struct_Name.Trace_Name.Microphone_Name" to
% Struct_Name.Microphone_Name, with each microphone field containing a
% matrix with a column for each trace. Each column contains that trace's
% peak frequencies. 

PeakFreqComp = [];

for i = 1:length(mNames)
    
    for j = 1:length(tNames);
        L(i) = length(peakStruc.(tNames{j}).(mNames{i}));
    end
    
    maxLen = max(L);
    zMat = zeros(maxLen,length(tNames));
    
    for jj = 1:length(tNames);
        
        zMat(1:length(peakStruc.(tNames{jj}).(mNames{i})),jj) = ...
            getfield(peakStruc,tNames{jj},mNames{i},{':',1}); 
    end
    
    PeakFreqComp.(mNames{i}) = zMat;
end

%% Create traceMatchall string
% This small bit of code creates a string to be used in the traceMatch
% column. It simply creates a string to be used when all the traces match.
% It is in the format '1,2,3...'

traceMatchall = '';

for i = 1:Nt
    traceMatchall = strcat(traceMatchall, num2str(i), ',');
end

traceMatchall = traceMatchall(1:end -1);

%% Configure Statistics Matrix
% Creates the output "stat" array. Labels the columns and establishes
% column indicators to be used in the main loop.

stat = num2cell(zeros(length(mNames) + 1, 5)); 

statMatch = 2;
statnonMatch = 3;
statAvDelt = 4;
statTrace = 5;
stat(2:end, 1) = mNames;
stat{1, 1} = [];
stat{1, 2} = '# Matched sets';
stat{1, 3} = '# Non-Matching Peaks';
stat{1, 4} = 'Average DeltaVal';
stat{1, 5} = 'Trace w/ Most Matches';

%% Main Loop
% Executes for each microphone. Here is where the PeakFreqComp matrix is
% made.

for i = 1:length(mNames)

t = 2 * dF;             % Tolerance (Hz)    

linecount = 1;          % Where the current line of PeakFreqComp is stored. Also is equivalent to number of while loop iterations. A for loop is unusable here because the length of the loop is not determined until the loop has begun processing PeakFreqComp. 
tCount = 0;             % indicates how many times tolerance has switched. necessary for tolerance increases to be executed within the if tolInd statement.  
endRow = ones(1, Nt);   % established in order to keep dimensions consistent when inserting zeros into PeakFreqComp

Bool = zeros(2, 1);     % Creates empty variables with correct orientations. Values will be added
Delta = zeros(2, 1);    % to these during the main loop.
Average = zeros(2, 1);
StDev = zeros(2, 1);
Fabove = zeros(2, 1);
Fbelow = zeros(2, 1);
traceMatch = {};
statTraceMatch = zeros(1, length(tNames));

    while 1
        
        if tolInd       % Executed if dynamic tolerance is indicated. Will increase the tolerance based on frequency values.
            if any(PeakFreqComp.(mNames{i})(linecount,1) > 5000) && tCount < 1
                t = 6 * dF;
                tCount = 1;
            elseif any(PeakFreqComp.(mNames{i})(linecount,1) > 10000) && tCount < 2
                t = 10 * dF;
                tCount = 2;
            elseif any(PeakFreqComp.(mNames{i})(linecount,1) > 15000) && tCount < 3
                t = 14 * dF;
                tCount = 3;
            end
        end
        
        [tempMin, ~] = min(nonzeros(PeakFreqComp.(mNames{i})(linecount,:)));   % finds min of each row
        minVec = ones(1, Nt) * tempMin;           % creates vector to be subtracted from row.
        tempdiffVec = abs(PeakFreqComp.(mNames{i})(linecount,:) - minVec);    % vector of difference of each cell in the row and the min of the row
        match = tempdiffVec <= t; % calculates logical match (within tolerance) vector. 1 is a match and stays on the line, 0 is no match and must be moved.
        tempTraceMatch = '';      % creates empty string to be used in the traceMatch column.
        
        done = (sum(match) == numel(nonzeros(PeakFreqComp.(mNames{i})(linecount,:))) && ...  % If all nonzero elements on line match within tolerance and the next line is all zeros,
            isequal(PeakFreqComp.(mNames{i})(linecount + 1,:), zeros(1, Nt))) || ...         % OR linecount is the same value as the length of the matrix, execution is terminated for current loop iteration.
            linecount == length(PeakFreqComp.(mNames{i})); 
                                                                                
        if ~ismember(ones(1, Nt), match, 'rows')    % ismember will return 1 if match contains all ones. SO, if everything does not match, then proceed to insert zeros. 
            for j = 1:Nt 
                if match(j) % if the cell matches, keep it on the line. if it doesnt, add a zero.
                    endRow(j) = 0; % required to keep matrix dimensions consistent.
             
                    tempTraceMatch = strcat(tempTraceMatch, num2str(j), ','); % builds string to be used in traceMatch column. Each time a cell matches, the trace that it is in will be added to this string.
               
                else % if the cell does not match, add a zero to its column. Store the last value in endRow and temporarily remove it from matrix in order to keep dimensions consistent since a zero is being added.
                    endRow(j) = PeakFreqComp.(mNames{i})(end, j);
                    
                    PeakFreqComp.(mNames{i})(:,j) = vertcat(PeakFreqComp.(mNames{i})(1:linecount - 1, j)...
                        , 0, PeakFreqComp.(mNames{i})(linecount:end - 1, j));
                end
               
            end
            
            PeakFreqComp.(mNames{i}) = vertcat(PeakFreqComp.(mNames{i}), endRow); % returns end values to columns where a zero was added.
            traceMatch{linecount} = tempTraceMatch(1:end - 1); % creates traceMatch column.
            
        else    
            traceMatch{linecount} = traceMatchall; % if everything matches, use traceMatchall for this line's traceMatch column.

        end    
        
      % get multitrace stats. see which trace had the most matches by
      % create a matchcount for each trace. statTraceMatch will be analyzed
      % outside of the while loop.
      
            if sum(match) > 1
                
            Indexx = find(match == 1);
            
                for jj = 1:length(Indexx)
                    statTraceMatch(Indexx(jj)) = statTraceMatch(Indexx(jj)) + 1;
                end
            end

        % Create Bool column
        if sum(PeakFreqComp.(mNames{i})(linecount,:) == 0) ~= Nt - 1  % If there is not one match, Boolean column for that line = 1.
           Bool(linecount) = 1;                                       % If there is one match, Boolean column = 0.
        else
           Bool(linecount) = 0;
        end
        
        % Create Average column
        Average(linecount) = mean(nonzeros(PeakFreqComp.(mNames{i})(linecount,:))); % Gets average of each line

        % Create stdev column
        StDev(linecount) = std(nonzeros(PeakFreqComp.(mNames{i})(linecount,:)));
        
        % Create Fabove column
        if linecount > 1
            Fabove(linecount) = Average(linecount) - Average(linecount - 1);
        elseif linecount == 1
            Fabove(linecount) = 0;
        end
        
        % Create Fbelow column
        if linecount > 1
            Fbelow(linecount - 1) = Average(linecount) - Average(linecount - 1);
            if done
                Fbelow(linecount) = 0;
            end
        end

        linecount = linecount + 1; % adds one to linecount for next iteration.
        
        if done  % If the done condition is true, terminate execution for this microphone.
            break
        end
        
    end
    
    % In order to account for added zeros, endRow was used to keep matrix
    % dimensions consistent. The result of this is a number of appended 
    % zeros at the end of the matrix. This piece of code finds the first
    % all zero row and removes it and every row after it. 
    
    endLine = zeros(1, Nt); 
    [~, Locb] = ismember(endLine, PeakFreqComp.(mNames{i}), 'rows');
    PeakFreqComp.(mNames{i}) = PeakFreqComp.(mNames{i})(1:Locb - 1, :);
    
    % The lines below build the PeakFreqComp matrix for this microphone 
    % from the additional columns created during the while loop. Converts 
    % it into a cell array in order to add the traceMatch column to the 
    % array.
    
    PeakFreqComp.(mNames{i}) = [PeakFreqComp.(mNames{i}) Bool Average StDev Fabove Fbelow]; 
    PeakFreqComp.(mNames{i}) = num2cell(PeakFreqComp.(mNames{i}));
    
    PeakFreqComp.(mNames{i}) = [PeakFreqComp.(mNames{i}) traceMatch'];
    
    %% Get Stats
    % Gets statistics and inserts them into "stat" array.
    
    [~, statTraceMatchInd] = max(statTraceMatch); % finds the index of the trace with the most matches for the current microphone.
    
    stat(i + 1, statMatch) = num2cell(sum(Bool));
    stat(i + 1, statnonMatch) = num2cell(sum(Bool == 0));
    stat(i + 1, statAvDelt) = num2cell(mean(nonzeros(Fabove)));
    stat(i + 1, statTrace) = tNames(statTraceMatchInd);
end

%% Create Green circle and red circle strucs
% This section creates structures of data to be plotted as red circles and
% green circles for each microphone. A green circle is an matched peak, 
% while a red circle is an unmatched peak.

GCfreq = [];
RCfreq = [];

for i = 1:Nt % For each trace
    for j = 1:length(mNames) % For each mic
        
        tempGCamp = []; % Create blank temporary variables
        tempRCamp = [];
        
        GCfreqIndex = find(cell2mat(PeakFreqComp.(mNames{j})(:, Nt + 1))); % Returns the index where the Boolean column == 1. These are the indices of rows with matches.
        RCfreqIndex = find(cell2mat(PeakFreqComp.(mNames{j})(:, Nt + 1)) == 0); % Returns the index where the Boolean column == 0. These are the indices of rows with non-matches.
        
        tempGCfreq = cell2mat(PeakFreqComp.(mNames{j})(GCfreqIndex, i)); % Uses the freqIndex variables to grab the value on each indexed row. Some of these will be zeros, since even though a peak may match there still may be one trace on the line that is a zero.
        tempRCfreq = cell2mat(PeakFreqComp.(mNames{j})(RCfreqIndex, i));
     
        GCfreq.(tNames{i}).(mNames{j}) = nonzeros(tempGCfreq); % Removes these potential zeros and stores the resulting vector in its respective GC or RC structure.
        RCfreq.(tNames{i}).(mNames{j}) = nonzeros(tempRCfreq);
        
        % Amplitude values for the circles will be necessary for plotting.
        % These loops will find the amplitude value for each GC and RC
        % frequency.
        
        for k = 1:length(GCfreq.(tNames{i}).(mNames{j})) % For each Green circle
            [~, GCampIndex] = ismember(GCfreq.(tNames{i}).(mNames{j})(k), FRF.(tNames{i}).(mNames{j})); % Gets the index where the exact GC frequency is located in the original FRF structure.
            tempGCamp(k) = FRF.(tNames{i}).(mNames{j})(GCampIndex, 2); % Uses this index to return the ampitude and store it in tempGCamp.
        end
        
        for ii = 1:length(RCfreq.(tNames{i}).(mNames{j})) % Same process is followed as the above GC for loop.
            [~, RCampIndex] = ismember(RCfreq.(tNames{i}).(mNames{j})(ii), FRF.(tNames{i}).(mNames{j}));
            tempRCamp(ii) = FRF.(tNames{i}).(mNames{j})(RCampIndex, 2);
        end
        
        GCfreqAmp.(tNames{i}).(mNames{j}) = [GCfreq.(tNames{i}).(mNames{j}) tempGCamp']; % Stores both the frequency and amplitude of GC's and RC's in their respective structures.
        RCfreqAmp.(tNames{i}).(mNames{j}) = [RCfreq.(tNames{i}).(mNames{j}) tempRCamp']; % These will be used to plot the circles.
    end
end


%% Plot Circles + Traces

for i = 1:length(mNames)
    
    figure
    for j = 1:Nt
        
        h(i,j) =  subplot(Nt, 1, j); % Stores the subplots in order to link axes.
        hold on
        title([strrep(mNames{i}, '_', ' '), ', ', tNames{j}])
        xlabel('Frequency (Hz)')
        ylabel('FRF Amplitude (dB)')
        
        plot(FRF.(tNames{j}).(mNames{i})(FStart:FStop, 1), FRF.(tNames{j}).(mNames{i})(FStart:FStop, 2), '.-', 'linewidth', 2) % Plots both traces
        
        % Plot green circles
        if ~isempty(GCfreqAmp.(tNames{j}).(mNames{i})) % logic is necessary in case there is a trace that never has a mismatch. If this happens, there will be no red circles, which causes an error when plotting. 
            plot(GCfreqAmp.(tNames{j}).(mNames{i})(:,1), GCfreqAmp.(tNames{j}).(mNames{i})(:,2), 'go', 'linewidth', 2)
        end
        
        % Plot red circles
        if ~isempty(RCfreqAmp.(tNames{j}).(mNames{i}))
            plot(RCfreqAmp.(tNames{j}).(mNames{i})(:,1), RCfreqAmp.(tNames{j}).(mNames{i})(:,2), 'ro', 'linewidth', 2)
        end
        
        legend('FRF Data', 'Matching Peak', 'Unmatched Peak') % Creates legend. May return a warning if there is a trace with no mismatches (or matches) since there will be only 2 data sets on the plot. Code will still execute regardless.
    end
    
    linkaxes(h(i,:),'x'); % Links axes to all obey the following xlim command. 
    xlim([freqStart 2000]); % Sets the default window on each plot. Can be adjusted to preference.
end

%% Export struc to excel UNFINISHED
if xlsExport
% this must be configured.    
% xlswrite('PeakFreqComp.Ext_M1.xlsx', PeakFreqComp.Ext_M1(:,1:end - 1))

fprintf('Done.\n\n')

end