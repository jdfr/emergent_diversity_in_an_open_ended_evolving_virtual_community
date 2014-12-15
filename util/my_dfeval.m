function varargout = my_dfeval(dfcn, varargin)
%customized version of dfeval

error(nargoutchk(1,inf,nargout, 'struct'))
error(nargchk(2,inf,nargin, 'struct'))% for argOutLp

%% Create, submit and return
jobObj = my_dfevalasync(dfcn,nargout,varargin{:});

waitForState(jobObj,'finished');
data = getAllOutputArguments(jobObj);

errorMessages = get(jobObj.Tasks, {'ErrorMessage'});
% Retain only the non-empty error messages, if any.
errorMessages = errorMessages(~cellfun(@isempty, errorMessages));

if ~isempty(errorMessages)
    % Indent the task error so that is visually separated from the error
    % we are throwing here.
    indMsg = regexprep(errorMessages{1}, '^(.)', '    $1', 'lineanchors');
    error('distcomp:dfeval:TaskErrored', ...
        ['Job %d encountered the following error:\n%s\n', ...
        'Data must be manually retrieved from the job, \n', ...
        'and it also needs to be manually destroyed.'], ...
        get(jobObj,'ID'), indMsg);
end

if size(data,2) < nargout
    % We only reach this situation in corner cases such as when the task
    % function throws errors with empty error messages.
    error('distcomp:dfeval:TooManyOutputArgs', ...
        ['Requested number of output arguments is greater than the number\n' ...
        'of output arguments returned by the tasks.  Data must be manually\n' ...
        'retrieved from job %d.\n', ...
        'The job also needs to be manually destroyed.'], get(jobObj,'ID'));
end

% The data comes back in a NxM array and it needs to be reformatted
% into M Nx1 outputs.
varargout = cell(1, nargout);
for i = 1:nargout
    varargout{i} = data(:, i);
end

% No error messages, and we were able to retrieve all the data that the
% user asked for, so we destroy the job.
jobObj.destroy;
