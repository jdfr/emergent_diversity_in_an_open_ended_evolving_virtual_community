function [ended states legal] = waitForJobs(ftst, jobs, decisionFunction, terclus) %#ok<INUSD>
%check whether jobs are finished or not. Try it until custom condition is
%met

  fprintf('WAITING FOR JOBS... ');

  if isempty(jobs)
    ended  = false(size(jobs));
    states =  cell(size(jobs));
    legal  =  true(size(jobs));
    return
  end
  legal = cellfun(@(x) isa(x.job, 'distcomp.job') || ...
                       isa(x.job, 'distcomp.simplejob') || ...
                          (isstruct(x.job) && ...
                           isfield(x.job, 'fakeJob') && ...
                           x.job.fakeJob), ...
                  jobs);
  if ~all(legal)
    notlegal = jobs(~legal);
    cls      = cellfun(@(x)[class(x) ', '], notlegal, 'uniformoutput', false);
    cls{end} = cls{end}(1:(end-2));
    cls      = horzcat(cls{:});
    mprintf(ftst, 'ERROR!!!! Some jobs are not ''distcomp.job'' but: %s\n', cls);
    st       = dbstack;
    cls      = cell(size(st));
    for k=1:numel(st)
      cls{k} = sprintf('   FILE: %s\n   NAME: %s\n   LINE: %d\n', st(k).file, st(k).name, st(k).line);
    end
    cls      = horzcat(cls{:});
    mprintf(ftst, '   THE STACK FOR THIS WARNING IS:\n%s', cls);
    if numel(jobs(legal))==0
      ended  = true(size(jobs));
      states = repmat({'notjob'}, size(jobs));
      return;
    end
  end
  %timeStart = clock;
  legalc =  arrayfun(@(x)x, legal, 'uniformoutput', false);
  while true 
    states  = cellfun(@safeQueryState, jobs, legalc, 'uniformoutput', false);
    ended   = cellfun(@stateEnd, states);
    if decisionFunction(ended) %do this until custom condition is met
      break;
    else
      %pause(terclus.pauseTime); %DAN: better than pooling the
      %cluster for the job state is to leave the cluster to awake you
      for i=1:size(jobs, 2)
        waitForState(jobs{i}.job, 'finished', 300);
      end
    end
  end
  
fprintf('DONE!\n');
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = safeQueryState(x, legal)
if legal
  s = x.job.State;
else
  s = 'notjob';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ended = stateEnd(state)
  switch state
    case {'finished', 'failed', 'destroyed', 'unavailable', 'notjob'}
      ended = true;
    otherwise
      ended = false;
  end
