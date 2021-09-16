function [oMeasure, oErrors, oFrameErrors, oBlackErrors] =...
    PerformanceTRA(gtPath, resPath, K, threshold)

% File to save the performance evaluation to, so that it does not have to
% be recomputed the next time.
resFile = fullfile(resPath, 'AOGMR_log.txt');

% File with AOGM results for an empty set of tracks.
resFileBlack = fullfile(gtPath, 'TRA', 'AOGMR_BLACK_log.txt');

% Compute the AOGM results for an empty set of tracks if necessary.
if ~exist(resFileBlack, 'file')
    AOGMR_BLACK(fullfile(gtPath, 'TRA'));
else
    delete(strtrim(resFileBlack));
    AOGMR_BLACK(fullfile(gtPath, 'TRA'));
end

% Compute performance data if necessary.
if exist(resFile, 'file')
    delete(strtrim(resFile));
end
if ~isfolder(resPath)
    ExportCellsTif(aSeqPaths, aTestVer)
end
AOGMR(resPath, fullfile(gtPath, 'TRA'), threshold)

% Read the files with AOGM-results.
[resMeasure, oErrors, oFrameErrors] = ReadAOGMFile(K, resFile);
[blackMeasure, oBlackErrors] = ReadAOGMFile(K, resFileBlack);

% Compuate the TRA measure.
oMeasure = max(1 - resMeasure / blackMeasure, 0);
end

function [oMeasure, oErrors, oFrameErrors] = ReadAOGMFile(K, aPath)
% Reads log files with AOGM-results.
%
% Inputs:
% aImData - ImageData object for the image sequence.
% aPath - Path of the log file with AOGM-results.
%
% Outputs:
% oMeasure - AOGM measure.
% oErrors - Error matrix with a single row. The matrix has the same
%           structure as the output oErrors returned by PerformanceTRA.
% oFrameErrors - Matrix with the errors in each frame. The matrix has the
%                same structure as the output oFrameErrors returned by
%                PerformanceTRA.

% Read the saved file with performance data as a single string.
fid = fopen(aPath, 'r');
results = fscanf(fid, '%c', inf);
fclose(fid);

% Extract the AOGM measure.
oMeasure = str2double(regexp(results,...
    '(?<=AOGM (measure|value): )[\d\.eE+]*', 'match'));

% Split the file content into lines.
lines = regexp(results, '\r\n', 'split');

% The names of the different error types, which appear on the rows above
% the errors of the given types.
headings = {...
    'False Negative Vertices'
    'False Positive Vertices'
    'Splitting Operations'
    'Edges To Be Added'
    'Redundant Edges To Be Deleted'
    'Edges with Wrong Semantics'}';

% Indices of lines with headings.
headLines = nan(size(headings));
for i = 1:length(headLines)
    headMatch = regexp(lines, ['.*' headings{i} '.*'], 'once');
    index = find(~cellfun(@isempty, headMatch));
    if ~isempty(index)
        headLines(i) = index;
    end
end

% Indices of lines describing errors.
times = regexp(lines, '(?<=T=)\d+', 'match', 'once');
times = cellfun(@str2double, times) + 1;

% Array of line indices where the lines which do not describe errors are
% NaNs.
timeLines = 1:length(times);
timeLines(isnan(times)) = nan;

% Count the number of lines describing errors after each heading.
oErrors = zeros(1,6);
oFrameErrors = zeros(K, 6);
for i = 1:length(headLines)  % Go through the error types.
    % Index of the line directly after the last error of the given type.
    nextHeadLine = min([headLines(headLines > headLines(i)) max(timeLines)+1]);
    % The frames that errors occurred in.
    errorTimes = times(timeLines > headLines(i) & timeLines < nextHeadLine);
    
    % Count the number of errors in each frame.
    oErrors(i) = length(errorTimes);
    for j = 1:length(errorTimes)
        oFrameErrors(errorTimes(j),i) = oFrameErrors(errorTimes(j),i) + 1;
    end
end
end