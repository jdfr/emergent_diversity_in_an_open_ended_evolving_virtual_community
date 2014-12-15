function doForSeveralSimulationsParallel(mode, fun, addargs, basedir, path, recursive, reorderfun, torquedata)
%for a series of simulations, launch parallel workers to perform tasks in
%each one of them

if ~exist('recursive', 'var')
  recursive = false;
end
if ~exist('path', 'var')
  path = '';
end
if ~exist('reorderfun', 'var')
  reorderfun = [];
end
if ~exist('torquedata', 'var')
  torquedata = [];
end

doForSeveralSimulations(@doParallel, basedir, torquedata, {mode, path, fun, array2cell(addargs)}, recursive, reorderfun);


function torquedata = doParallel(basedir, torquedata, mode, path, fun, addargs)

oldtorquedata = torquedata;

if strcmp(mode, 'torque')
  if not(isfield(torquedata, 'jobArgs'))
    torquedata.jobArgs = {};
  end
  if not(isfield(torquedata, 'jobMgr'))
    torquedata.jobMgr = findResource('scheduler', 'type', 'torque');
  end
  if isfield(torquedata, 'tag')
    torquedata.jobArgs(end+[1 2]) = {'Tag', torquedata.tag};
  end
  torquedata.jobArgs(end+[1 2]) = {'PathDependencies',  genpath(path)};
end

switch mode
  case 'embedded'
    fun(basedir, addargs{:});
  case 'torque'
    job = my_dfevalasyncJM(...
        torquedata.jobMgr, ...
        fun, ...
        1, ...
        {basedir}, ...
        addargs{:}, ...
        torquedata.jobArgs{:} ...
        ); %#ok<NASGU>
    %error('mode not implemented yet: %s', mode);
  otherwise
    error('mode not understood: %s', mode);
end

torquedata = oldtorquedata;
