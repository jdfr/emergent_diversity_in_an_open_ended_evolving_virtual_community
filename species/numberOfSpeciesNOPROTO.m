function numberOfSpeciesNOPROTO(basedir,GenI,GenF,parallel, prefix, path, options)
if not(exist('parallel', 'var'))
  parallel = [];
end
if not(exist('prefix', 'var'))
  prefix = '';
end
if not(exist('path', 'var'))
  path = '';
end
if not(exist('options', 'var'))
  options = struct;
end
opts = {'oldWayToFindPrototype', true; 'thresholdMode', 'absolute'; 'threshold', 1};
for k=1:size(opts,1)
  if not(isfield(options, opts{k,1}))
    options.(opts{k,1}) = opts{k,2};
  end
end

if isempty(GenF)
  generations = GenI(:)';
else
  generations = (GenI:GenF);
end

if isempty(parallel)
  for i = generations%GenI : GenF
    doAGeneration(basedir, i, prefix, options);
  end
elseif islogical(parallel)
  error('''Parallel'' argument must be the number of jobs!!!!!');
else
  generations = onlyNonProcessed(basedir, prefix, generations);
  if isempty(generations)
    fprintf('numberOfSpeciesNOPROTO prefix %s already done for %s!!!\n', prefix, basedir);
    return
  else
    fprintf('NUMGENERATIONS: %d\n', numel(generations));
  end
  if iscell(parallel)
    if isempty(parallel{1})
      parallel = numel(generations);
    else
      parallel = ceil(numel(generations)/parallel{1});
    end
  end
  jobArgs = {'Tag', 'numberOfSpecies', 'PathDependencies',  {path}};
  gs = generations';%(GenI:GenF)';
  [ranges rangesizes] = calculateRanges(numel(gs), parallel);
  gens = mat2cell(gs, rangesizes, 1);
  bass = repmat({basedir}, size(gens));
  prefixes = repmat({prefix}, size(gens));
  opts = repmat({options}, size(gens));
  fprintf('GOING PARALLEL (%d JOBS)...\n', numel(bass));
  if numel(generations)<10
    fprintf('only these gens: %s\n', any2str(gs));
  end
  dummy = my_dfeval(@doSeveralGenerations, bass, gens, prefixes, opts, jobArgs{:}); %#ok<NASGU>
  fprintf('RECEIVING PARALLEL...\n');
end

function  gens = onlyNonProcessed(basedir, prefix, gens)
calc = false(size(gens));
pr = [basedir filesep prefix];
for k=1:numel(gens)
    n = num2str(gens(k), '%03g');
    calc(k) = not( exist([pr 'sim' n '.mat'], 'file') && exist([pr 'numGroups' n '.mat'], 'file') && exist([pr 'individualsSpecies' n '.mat'], 'file') );
end
gens = gens(calc);

function dummy = doSeveralGenerations(basedir, gens, prefix, options)
for i=gens(1):gens(end)
  doAGeneration(basedir, i, prefix, options);
end
dummy = gens;


function doAGeneration(basedir, i, prefix, options)
    tosave1 = ['individualsSpecies' num2str(i, '%03g')];
    tosave4 = ['sim' num2str(i, '%03g')];
    tosave5 = ['numGroups' num2str(i, '%03g')];
    if exist([basedir filesep prefix tosave1 '.mat'], 'file') && exist([basedir filesep prefix tosave4 '.mat'], 'file') && exist([basedir filesep prefix tosave5 '.mat'], 'file')
      return;
    end
    fprintf('G: %d\n', i);
    load(sprintf('%s%s%sZ%03d.mat',basedir,filesep,prefix,i));
    Z = eval(sprintf('Z%03d', i));
    load(sprintf('%s%s%smatrixDist%03d.mat',basedir,filesep,prefix,i));
    matrixDist = eval(sprintf('matrixDist%03d', i));
    load(sprintf('%s%s%sindicesClasificados%03d.mat',basedir,filesep,prefix,i));
    indicesClasificados = eval(sprintf('indicesClasificados%03d', i));

    [individualsSpecies sim numGroups] = numberGroupNOPROTO(basedir,i,Z,matrixDist,indicesClasificados, options); %#ok<NASGU>

    

    eval([tosave1 ' = individualsSpecies;']);
    eval([tosave4 ' = sim;']);
    eval([tosave5 ' = numGroups;']);

    save([basedir filesep prefix tosave1 '.mat'],tosave1);
    save([basedir filesep prefix tosave4 '.mat'],tosave4);
    save([basedir filesep prefix tosave5 '.mat'],tosave5);

    clear('individualsSpecies','sim','numGroups', tosave1,tosave4,tosave5);
