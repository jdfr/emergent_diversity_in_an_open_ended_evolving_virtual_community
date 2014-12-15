function job = recalculateClustersScript(indexOpts, indexType, option, range)

if not(exist('option', 'var'))
  option=0;
end

if not(exist('indexOpts', 'var')) || isempty(indexOpts)
  indexOpts = 5:8;
end

if not(exist('indexType', 'var')) || isempty(indexType)
  indexType = 1:3;
end

if ispc
  base = 'lsystemdani\dataset';
  subs = {'G1F1\uno\1=0.001_3', 'G0.75F1\uno\1=0.001_1', 'G0.5F1\dos\1=0.001_1'};
else
  base = 'rene/dataset';
  subs = {'G1F1/uno/1=0.001_3', 'G0.75F1/uno/1=0.001_1', 'G0.5F1/dos/1=0.001_1'};
end

dirs = cellfunc(@(x)[base filesep x], subs);
labels = cellfunc(@(x)x(1:find(x==filesep, 1)), subs);

numberOfJobs = [];
path = genpath('rene/src');

jobArgs={'Tag', 'makeAllClusterings', 'PathDependencies',  path};

combinations = cell(numel(dirs)*numel(indexOpts)*numel(indexType), 5);
combinations(1:end,1) = {numberOfJobs};
combinations(1:end,2) = {path};
comblabels   = cell(size(combinations,1), 1);
k=1;
for a=1:numel(dirs)
  for b=1:numel(indexType)
    for c=1:numel(indexOpts)
      combinations{k,3} = dirs(a);
      comblabels{k}     = labels{a};
      combinations{k,4} = indexOpts(c);
      combinations{k,5} = indexType(b);
      k=k+1;
    end
  end
end

if option>1
  opts = allOptions;
  opts = opts(indexOpts,2);
  posfix = cellfunc(@(x)[x.thresholdMode mat2str(x.threshold)], opts);
  switch option
    case 2
      for k=1:numel(posfix)
        fprintf('GOING FOR %s...\n', posfix{k});
        figuras = makeFigurasStructNEW(dirs, @(x)['R_' x posfix{k} '_']); %#ok<NASGU>
        save([base filesep 'figuras_' posfix{k} '.mat'], 'figuras');
      end
    case 4
      for k=1:numel(posfix)
        fprintf('GOING FOR %s...\n', posfix{k});
        oldprefix = 'absolute1';
        newprefix = posfix{k};
        prefixfun = @(x)['R_' x newprefix '_'];
        makeFigurasStructOnlyNewGroups(base, oldprefix, newprefix, prefixfun);
      end
    case 3
      if not(ischar(range))
        range = posfix{range};
      end
      figuras = load([base filesep 'figuras_' range '.mat']);
      figuras = figuras.figuras;
            nms = 'diversity';
            tit = sprintf('Threshold criterion: %s', range);
            h=figurassNuevo(figuras, 1:3, 'ed', true, [0 0 1 0 0]); %#ok<NASGU>
            set(get(gca, 'YLabel'), 'string', [nms ' (edit dist.)']);
            title(tit);
            %imgs{end+1,1} = imcapture(h, 'all', 600); close(h);
            h=figurassNuevo(figuras, 1:3, 'edred', true, [0 0 1 0 0]); %#ok<NASGU>
            set(get(gca, 'YLabel'), 'string', [nms ' (red. edit dist.)']);
            title(tit);
            %imgs{end+1,1} = imcapture(h, 'all', 600); close(h);
            h=figurassNuevo(figuras, 1:3, 'morpho', true, [0 0 1 0 0]); %#ok<NASGU>
            set(get(gca, 'YLabel'), 'string', [nms ' (Jaccard distance)']);
            title(tit);
            %imgs{end+1,1} = imcapture(h, 'all', 600); close(h);
            %img = vertcat(imgs{:});
            %imwrite(img, names{k}{1});
    otherwise
      error('jarllll!!!!');
  end
  job = [];
  return
end

if option==1
  for k=1:numel(range)
    recalculateClustersScriptMini(combinations{range(k),:});
    job = [];
  end
elseif option==0
  job = my_dfevalasync(@recalculateClustersScriptMini, 1, combinations(:,1), combinations(:,2), combinations(:,3), combinations(:,4), combinations(:,5), jobArgs{:});
elseif option==-1
  for k=1:numel(range)
    recalculateNumGroupsOnly(combinations{range(k),3}, combinations{range(k),4}, combinations{range(k),5});
  end
  job = [];
else
  error('jarlll');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opts = allOptions
opts = {...
  01 struct('oldWayToFindPrototype', false, 'thresholdMode', 'absolute',   'threshold', 1); ...
  02 struct('oldWayToFindPrototype', false, 'thresholdMode', 'abslast',    'threshold', 1); ...
  03 struct('oldWayToFindPrototype', false, 'thresholdMode', 'absolute',   'threshold', 0.1); ...
  04 struct('oldWayToFindPrototype', false, 'thresholdMode', 'abslast',    'threshold', 0.1); ...
  05 struct('oldWayToFindPrototype', false, 'thresholdMode', 'relval',     'threshold', 0.05); ...
  06 struct('oldWayToFindPrototype', false, 'thresholdMode', 'relvallast', 'threshold', 0.05); ...
  07 struct('oldWayToFindPrototype', false, 'thresholdMode', 'relval',     'threshold', 0.01); ...
  08 struct('oldWayToFindPrototype', false, 'thresholdMode', 'relvallast', 'threshold', 0.01); ...
  09 struct('oldWayToFindPrototype', false, 'thresholdMode', 'absolute',   'threshold', 0.01); ...
  10 struct('oldWayToFindPrototype', false, 'thresholdMode', 'abslast',    'threshold', 0.01); ...
  11 struct('oldWayToFindPrototype', false, 'thresholdMode', 'absolute',   'threshold', 0.05); ...
  12 struct('oldWayToFindPrototype', false, 'thresholdMode', 'abslast',    'threshold', 0.05); ...
  13 struct('oldWayToFindPrototype', false, 'thresholdMode', 'reldiff',     'threshold', 0.05); ...
  14 struct('oldWayToFindPrototype', false, 'thresholdMode', 'reldifflast', 'threshold', 0.05); ...
  15 struct('oldWayToFindPrototype', false, 'thresholdMode', 'reldiff',     'threshold', 0.01); ...
  16 struct('oldWayToFindPrototype', false, 'thresholdMode', 'reldifflast', 'threshold', 0.01); ...
  };


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dummy = recalculateClustersScriptMini(numberOfJobs, path, dirs, indexOpts, indexType)

dummy = 1;

opts = allOptions;
opts = opts(indexOpts,2);
for z=1:numel(opts)
  newpref = [opts{z}.thresholdMode mat2str(opts{z}.threshold)];
  data = {...%'nummuts'  'NM_';...
          'morpho'   'MP_' ['R_MP_' newpref '_'];...
          'edit'     'E1_' ['R_E1_' newpref '_'];...
          'miniedit' 'E2_' ['R_E2_' newpref '_'];...
          };
  data = data(indexType,:);
  prefixesOLD = data(:,2);
  prefixesNEW = data(:,3);
  for k=1:numel(dirs)
    recalculateClusters(dirs{k}, prefixesOLD, prefixesNEW, numberOfJobs, path, opts{z});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dummy = recalculateNumGroupsOnly(dirs, indexOptsNew, indexType)

dummy = 1;

opts = allOptions;
opts = opts(:,2);
for z=1:numel(indexOptsNew)
  oldpref = 'absolute1';%[opts{indexOptsOld(z)}.thresholdMode mat2str(opts{indexOptsOld(z)}.threshold)];
  newpref = [opts{indexOptsNew(z)}.thresholdMode mat2str(opts{indexOptsNew(z)}.threshold)];
  data = {...%'nummuts'  'NM_';...
          'morpho'   ['R_MP_' oldpref '_'] ['R_MP_' newpref '_'];...
          'edit'     ['R_MP_' oldpref '_'] ['R_E1_' newpref '_'];...
          'miniedit' ['R_MP_' oldpref '_'] ['R_E2_' newpref '_'];...
          };
  data = data(indexType,:);
  prefixesOLD = data(:,2);
  prefixesNEW = data(:,3);
  for k=1:numel(dirs)
    fprintf('GOING FOR (%s), %s...\n', newpref, dirs{k});
    sim = load([dirs{k} filesep prefixesOLD{k}  'sim.mat']);
    numGroups = cellfunc(@calculateNumGroups, sim.sim, repmat(opts(indexOptsNew(z)), size(sim.sim))); %#ok<NASGU>
    save([dirs{k} filesep prefixesNEW{k}  'numGroups.mat'], 'numGroups');
  end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeFigurasStructOnlyNewGroups(base, oldprefix, newprefix, prefixfun)
clusterings = {'morphological', prefixfun('MP_') 11; ...
               'edit distance', prefixfun('E1_') 12; ...
               'reduced edit distance', prefixfun('E2_') 13; ...
               };
prefixes = clusterings(:,2);
prefixesi = vertcat(clusterings{:,3});
clusterings = clusterings(:,1);

figuras = load([base filesep 'figuras_' oldprefix '.mat']);
figuras = figuras.figuras;


indexes = 11:13;

for k=1:numel(indexes)
  for z=1:size(figuras.tabla,2)
    nm = [figuras.tabla{1,z} filesep prefixes{k} 'numGroups.mat'];
    numGroups = load(nm);
    numGroups = numGroups.numGroups;
    if iscell(numGroups)
      numGroups = [numGroups{:}]';
    end
    fprintf('(%s,%s)=%s\n', mat2str(indexes(k)), mat2str(z), nm);
    figuras.tabla{indexes(k),z} = numGroups;
    clear numGroups;
  end
end


save([base filesep 'figuras_' newprefix '.mat'], 'figuras');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figuras = makeFigurasStructNEW(basedirs, prefixfun, colors, indexes)

tipos = {'very harsh', 1; ...
         'harsh', 0.8; ...
         'harsh', 0.75; ...
         'semiharsh', 0.7;...
         'mild', 0.5; ...
         };
tiposG = [tipos{:,2}];
tipos = tipos(:,1);
campos = {01, 'nombres'; ...
          02, 'color'; ...
          03, 'order'; ...
          04, 'legend'; ...
          05, 'alphaparam'; ...
          06, 'betaparam'; ...
          07, 'generations'; ...
          08, 'popsizes'; ...
          09, 'numRamas'; ...
          10, 'speciesGenerations'; ...
          11, 'numClusters_morpho'; ...
          12, 'numClusters_edit'; ...
          13, 'numClusters_editreduced'; ...
          };
clusterings = {'morphological', prefixfun('MP_') 11; ...
               'edit distance', prefixfun('E1_') 12; ...
               'reduced edit distance', prefixfun('E2_') 13; ...
               };
prefixes = clusterings(:,2);
prefixesi = vertcat(clusterings{:,3});
clusterings = clusterings(:,1);
        
if not(exist('colors', 'var'))
  colors = mat2cell(lines(numel(basedirs)), ones(numel(basedirs),1), 3);
end
if not(exist('indexes', 'var'))
  indexes = array2cell(1:numel(basedirs));
end
if not(iscell(indexes))
  indexes = array2cell(indexes);
end
tabla = cell(numel(campos), numel(basedirs));

tabla(1,:) = basedirs(:)';
tabla(2,:) = colors(:)';
tabla(3,:) = indexes(:)';


for k=1:numel(basedirs)
  b = basedirs{k};
  params = load([b filesep '..' filesep 'estado.mat']);
  params = params.ACIParams;
  g = find(tiposG==params.alphaG);
  tabla{4,k} = tipos{g};
  tabla{5,k} = params.alphaG;
  tabla{6,k} = params.alphaF;
  fprintf('Going for %s\n', b);
  fn = [b filesep 'poph.mat'];
  if exist(fn, 'file')
    poph = load(fn);
    poph = poph.poph;
  else
    poph = loadSimulation(b);
  end
  mng = min(poph.generation);
  mxg = max(poph.generation);
  gs = (mng:mxg)';
  tabla{7,k} = gs;
  pops = zeros(size(gs));
  bios = zeros(size(gs));
  fprintf('   pops+bios...\n');
  for z=1:numel(gs)
    thisg = find(poph.generation==gs(z));
    thisg = thisg(poph.nPixelsInSoil(thisg)==0);
    pops(z) = numel(thisg);
    bios(z) = sum(poph.nbranches(thisg));
  end
  tabla{8,k} = pops;
  tabla{9,k} = bios;
    gss = (1:mxg)';
    tabla{10,k} = gss;
    for z=1:numel(prefixes)
      fprintf('   numSpecies %s...\n', clusterings{z});
      gemafile = [b filesep prefixes{z} 'numGroups.mat'];
      if exist(gemafile, 'file')
        numGroups = load(gemafile);
        numGroups = numGroups.numGroups;
        tabla{prefixesi(z),k} = [numGroups{:}]';
      else
        error('We need this file!!!!: <%s>', gemafile);
      end
    end
end

figuras = struct('campos', {campos}, ...
                 'tipos', {tabla(4,:)}, ...
                 'tabla', {tabla});
