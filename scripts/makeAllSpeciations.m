function makeAllSpeciations(basedir, geni, genf, path, allparallel, toUse, subprefix, optionsHC, optionsNG)
if not(exist('path', 'var'))
  path = pwd;
end
if not(exist('allparallel', 'var'))
  allparallel = false;
end
if not(exist('subprefix', 'var'))
  subprefix = '';
end
if not(exist('optionsHC', 'var'))
  optionsHC = [];
end
if not(exist('optionsNG', 'var'))
  optionsNG = [];
end

if isempty(optionsHC)
  optionsHC = {};
else
  optionsHC = {optionsHC};
end
if isempty(optionsNG)
  optionsNG = {};
else
  optionsNG = {optionsNG};
end

if allparallel
  jobArgs={'Tag', 'makeAllSpeciations', 'PathDependencies',  {path}};
end
data = {'nummuts'  'NM_';...
        'edit'     'E1_';...
        'miniedit' 'E2_';...
        'morpho'   'MP_'};
if exist('toUse', 'var')
  data = data(toUse,:);
end
modes = data(:,1);
prefixes = data(:,2);
if not(isempty(subprefix))
  prefixes = cellfun(@(x)[x(1:2) subprefix x(3)], prefixes, 'uniformoutput', false);
end
parallel = true;
ngs = genf-geni+1;
nj_species = {[]};%{20};%{20};%{[]};
fprintf('Processing speciations for %s\n', basedir);
fprintf('generating additional data\n');
if allparallel
  fprintf('POSSIBLY GENERATING ADDITIONAL DATA...\n');
  lastgen = generateAdditionalForHierarchical(basedir);
  fprintf('DONE...\n');
  jobs = cell(numel(modes),1);
else
  fprintf('POSSIBLY GENERATING ADDITIONAL DATA (parallel)...\n');
  lastgen = my_dfeval(@generateAdditionalForHierarchical, {basedir}, jobArgs{:});
  lastgen = lastgen{1};
  fprintf('DONE...\n');
end
if isempty(genf)
  genf = lastgen;
end
for k=1:numel(modes)
  if exist([basedir filesep prefixes{k} 'numGroups.mat'], 'file')
    fprintf('MODE %s ALREADY DONE!!!!\n', modes{k});
    continue;
  end
  fprintf('speciation: %s, clustering\n', modes{k});
  HierarchicalClusteringALL(basedir,[],modes{k},geni,genf,parallel,prefixes{k}, path, optionsHC{:});
  fprintf('speciation: %s, numSpecies\n', modes{k});
  numberOfSpeciesNOPROTO(basedir,geni,genf,nj_species, prefixes{k}, path, optionsNG{:});
  if not(allparallel)
    fprintf('COMPRESSION OF GEMAFILES...\n');
    compressGemaFiles(basedir, geni, genf, prefixes{k}, 12);
  else
    fprintf('ASYNC COMPRESSION OF GEMAFILES...\n');
    jobs{k} = my_dfevalasync(@compressGemaFiles2, 1, {basedir}, {geni}, {genf}, prefixes(k), {12}, jobArgs{:});
  end
end
for k=1:numel(modes)
  if isfloat(jobs{k})
    fprintf('STRANGE, NO GEMAJOB FOR PREFIX %s!!!\n', prefixes{k});
    continue;
  end
  fprintf('WAITING FOR ASYNC GEMAJOB %k (PREFIX %s)\n', k, prefixes{k});
  waitForState(jobs{k}, 'finished');
  nt = numel(jobs{k}.Tasks);
  erase = true;
  for z=1:nt
    if not(isempty(jobs{k}.Tasks(z).ErrorMessage))
      f = jobs{k}.Tasks(z).Function;
      if not(ischar(f))
        f = func2str(f);
      end
      fprintf('Job %d, task %d, with Tag %s and function %s with Error!!!!\n', k, z, jobs{k}.Tag, f);
      showError(jobs{k}.Tasks(z).Error);
      erase = false;
      break;
    end
  end
  if erase    
    destroy(jobs{k});
  end
end

function dummy = compressGemaFiles2(varargin)
dummy = 0;
compressGemaFiles(varargin{:});