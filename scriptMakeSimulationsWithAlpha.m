function joba = scriptMakeSimulationsWithAlpha(addparams, subdir, path, torqueArgs, randpermute, ntimes, alphas, remote)%(basedir, subdir, alphas, ntimes, opts)
%dispatch a batch of simulations with varying alphas and other parameters.
%If the simulation dies too early (90 generations) it is respawned until it
%passed beyond that wall.
%joba = scriptMakeSimulationsWithAlpha({'generaciones', '1000', 'stopIfTooTall', '100000'}, 10, (1:-0.05:0.5)');  
%joba = scriptMakeSimulationsWithAlpha({}, 'testLONG', 'rene/src', {}, true, 10, [0.75 1]);  

%subdirectory of rene/dataset where simulations are
%generated
if not(exist('subdir', 'var'))
  path = 'test110NB';
end

%path of the source code
if not(exist('path', 'var'))
  path = 'rene/src';
end

%arguments qsub
if not(exist('torqueArgs', 'var'))
  torqueArgs = '';
end

%permute at random the files
if not(exist('randpermute', 'var'))
  randpermute = false;
end

%config params for the simulation (to be passed over to main)
if not(exist('addparams', 'var'))
  addparams = {};
end

%number of simulations per alpha
if not(exist('ntimes', 'var'))
  ntimes = 10;
end

%for each value in this vector, we will do "ntimes" simulations with this
%alpha value
if not(exist('alphas', 'var'))
  alphas = (0.5:0.05:1)';
end

%do it remotely
if not(exist('remote', 'var'))
  remote = true;
end

alphas = alphas(:);
alphas = repmat(alphas, 1, ntimes)';
alphas = alphas(:);
nums   = (1:numel(alphas))';
%well mixed, please
if islogical(randpermute)
  if randpermute
    perm   = randperm(numel(alphas));
    alphas = alphas(perm);
    nums   = nums(perm);
  end
else
  if numel(randpermute)==numel(alphas)
    alphas = alphas(randpermute);
    nums   = nums(randpermute);
  elseif not(isempty(randpermute))
    error('randpermute not understood!!!!!');
  end
end

alphas = array2cell(alphas);
nums   = array2cell(nums);
addparams = repmat({addparams}, size(alphas));
subdir = repmat({subdir}, size(alphas));
remotes = repmat({remote}, size(alphas));
jobArgs = {'Tag', 'alphaSims', 'PathDependencies',  genpath(path)};

jobMgr = findResource('scheduler', 'type', 'torque');

if not(isempty(torqueArgs))
  oldargs = get(jobMgr, 'SubmitArguments');
  set(jobMgr, 'SubmitArguments', torqueArgs);
end

if remote
  joba = my_dfevalasyncJM(jobMgr, @scriptMakeSim, 1, subdir, alphas, nums, addparams, remotes, jobArgs{:});
else
  for k=1:numel(alphas)
    scriptMakeSim(subdir{k}, alphas{k}, nums{k}, addparams{k}, remotes{k});
  end
end

if not(isempty(torqueArgs))
  set(jobMgr, 'SubmitArguments', oldargs);
end

function res = scriptMakeSim(subdir, alpha, num, addparams, remote)

res = [];

newdir = [getNameByTime '_N' num2str(num, '%03g')];

base = 'rene';

origdir = 'rene/dataset/G0.5F1/dos';%/1=0.001_1';

opts       = {'ensayos',             '1', ...
              'prob',           '[0 0.05 0.05 0.05 0.05 0.05 0.05 0]', ...
              'generaciones',   '20000', ...
              ...'plottingFactor', '110', ...
              'alphaG',  mat2str(alpha), ...
              'alphaF', '1', ...
              'tam_poblacion',  '1', ...
              'stopIfTooTall',  '20000', ...
              'tooLongEvTime',  '15000', ...
              'initialPos',     '5000', ...
              'PRESIONES',      '0.001', ...
              'terclus.ranges', '[1 1]', ...
              'initialGeneration', '0', ...
              'initialPopulation', '{''G''}', ...
              '[correc.zeroHeightOK]',      'true', ...
              '[preserveUniqueTree]', 'false', ...
              '[T.sameLevelRule]', 'false', ...
              '[save.disp.gnm]', 'true', ... %way too costly
              '[save.disp.mng]', 'true', ... %not so costly, but useless if not recording any of the other two
              '[save.disp.pht]', 'true', ... %way too costly
              '[plotAtEnd]', '[true true true]', ...
              'plotDendroTree',    'false'};

subdir = [base filesep 'dataset' filesep subdir filesep 'G' num2str(alpha)];

% if remote
%   system('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/matlab/sys/opengl/lib/glnxa64');
% end
goon = true;
numIntentos = 1;
while goon
  params = main('embedded', {subdir, 'replay'}, {origdir}, newdir, 'seedRandom', {sum(clock)*1e5}, '[numTrials]', mat2str(numIntentos), opts{:}, addparams{:});
  nomdir = params.completedirs{1};
  poph = loadSimulation(nomdir, false);
  maxgen = max(poph.generation);
  clear poph;
  goon = maxgen<90;
  numIntentos = numIntentos + 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ahora = getNameByTime

ahora = datestr(now, 'yyyy_mmm_dd_HH_MM_SS');

