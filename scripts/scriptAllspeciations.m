function varargout = scriptAllspeciations(option, whattodo)

%sobre las plantas:
%-Las simulaciones con G0.75 son las que producen mas variedad dentro de
%las que son de tipo harsh. No obstante, hay una probabilidad elevada de
%que las plantas peque�as sucumban a las grandes y cese el crecimiento 
%(solo van quedando unos pocos arboles grandes dispersos, que mueren poco
%a poco)
%-si se quitan las duplicaciones tandem y duplicacion en nivel, los arboles
%se hacen mas plumosos, pero a costa de una probabilidad mucho mas elevada
%de desaparicion a largo plazo (supongo que por tener que mantener un
%mayor ramaje y mayor probabilidad de curvarse por debajo de tierra)

if nargin<2
  whattodo = 2;
end

if numel(whattodo) >1
  for k=1:numel(whattodo)
    scriptAllspeciations(option, whattodo(k));
  end
end

%alldirs = getAllDirs('rene/dataset', true);

if ischar(option)
  dirs = {option};
elseif iscell(option)
  dirs = option;
else

  switch option
    case -2.1
  dirs = { ...
    'rene/src/morpho/2009_Jul_29_12_16_39' ; ...
    'rene/src/morpho/2009_Jul_29_12_21_07' ; ...
    'rene/src/morpho/2009_Jul_29_14_22_12' ; ...
    };
    case -2
  dirs = { ...
    'rene/src/2009_Jul_29_12_16_39/3/1=0.001_1' ; ...
    'rene/src/2009_Jul_29_12_21_07/2/1=0.001_1' ; ...
    'rene/src/2009_Jul_29_14_22_12/5/1=0.001_1' ; ...
    'rene/src/miniedit/2009_Jul_29_12_16_39/1=0.001_1' ; ...
    'rene/src/miniedit/2009_Jul_29_12_21_07/1=0.001_1' ; ...
    'rene/src/miniedit/2009_Jul_29_14_22_12/1=0.001_1' ; ...
    'rene/src/mutdist/2009_Jul_29_12_16_39/3/1=0.001_1' ; ...
    'rene/src/mutdist/2009_Jul_29_12_21_07/2/1=0.001_1' ; ...
    'rene/src/mutdist/2009_Jul_29_14_22_12/5/1=0.001_1' ; ...
    };
    case -1
      dirs = {'rene/complexityfig/G0.7F0.5/2009_Jul_05_17_08_32/1=0.0001_4'};
    case 0
      dirs = { ...
        ...%'rene/dataset/G0.5F0.1/dos/1=0.001_1' ; ...
        ...%'rene/dataset/G0.5F0.5/dos/1=0.001_1' ; ...
        'rene/dataset/G0.5F0.5/uno/1=0.001_1' ; ...
        'rene/dataset/G0.5F0.5/uno/1=0.001_2' ; ...
        'rene/dataset/G0.5F0.5/uno/1=0.001_3' ; ...
        ...%'rene/dataset/G0.5F1/dos/1=0.001_1' ; ...
        'rene/dataset/G0.5F1/uno/1=0.001_1' ; ...
        'rene/dataset/G0.5F1/uno/1=0.001_2' ; ...
        'rene/dataset/G0.5F1/uno/1=0.001_3' ; ...
        'rene/dataset/G0.7F0.5/uno/1=0.001_1' ; ...
        'rene/dataset/G0.7F0.5/uno/1=0.001_2' ; ...
        'rene/dataset/G0.7F0.5/uno/1=0.001_3' ; ...
        'rene/dataset/G0.7F1/uno/1=0.001_1' ; ...
        'rene/dataset/G0.7F1/uno/1=0.001_2' ; ...
        'rene/dataset/G0.8F0.1/dos/1=0.001_1' ; ...
        'rene/dataset/G0.8F0.5/uno/1=0.001_1' ; ...
        'rene/dataset/G0.8F1/uno/1=0.001_1' ; ...
        'rene/dataset/G1F0.5/uno/1=0.001_1' ; ...
        'rene/dataset/G1F0.5/uno/1=0.001_2' ; ...
        'rene/dataset/G1F0.5/uno/1=0.001_3' ; ...
        'rene/dataset/G1F1/uno/1=0.001_1' ; ...
        'rene/dataset/G1F1/uno/1=0.001_2' ; ...
        'rene/dataset/G1F1/uno/1=0.001_3' ; ...
        };
    case 1
      dirs = {...
        ...%'rene/dataset/G0.5F0.1/dos/1=0.001_1' ; ...
        ...%'rene/dataset/G0.5F0.5/dos/1=0.001_1' ; ...
        'rene/dataset/G0.5F0.5/uno/1=0.001_1' ; ...
        'rene/dataset/G0.5F0.5/uno/1=0.001_2' ; ...
        'rene/dataset/G0.5F0.5/uno/1=0.001_3' ; ...
        ...%'rene/dataset/G0.5F1/dos/1=0.001_1' ; ...
        'rene/dataset/G0.5F1/uno/1=0.001_1' ; ...
        'rene/dataset/G0.5F1/uno/1=0.001_2' ; ...
        'rene/dataset/G0.5F1/uno/1=0.001_3' ; ...
        };
    case 2
      dirs = { ...
        'rene/dataset/G0.7F0.5/uno/1=0.001_1' ; ...
        'rene/dataset/G0.7F0.5/uno/1=0.001_2' ; ...
        'rene/dataset/G0.7F0.5/uno/1=0.001_3' ; ...
        'rene/dataset/G0.7F1/uno/1=0.001_1' ; ...
        'rene/dataset/G0.7F1/uno/1=0.001_2' ; ...
        'rene/dataset/G0.8F0.1/dos/1=0.001_1' ; ...
        'rene/dataset/G0.8F0.5/uno/1=0.001_1' ; ...
        'rene/dataset/G0.8F1/uno/1=0.001_1' ; ...
        };
    case 3
      dirs = {
        'rene/dataset/G1F0.5/uno/1=0.001_1' ; ...
        'rene/dataset/G1F0.5/uno/1=0.001_2' ; ...
        'rene/dataset/G1F0.5/uno/1=0.001_3' ; ...
        'rene/dataset/G1F1/uno/1=0.001_1' ; ...
        'rene/dataset/G1F1/uno/1=0.001_2' ; ...
        'rene/dataset/G1F1/uno/1=0.001_3' ; ...
        };
    case 4
      dirs = {      'rene/dataset/G0.8F0.5/uno/1=0.001_1' };
    case 5
      dirs = {  'rene/dataset/G0.8F1/uno/1=0.001_1' };
    case 6
  dirs = { ...
    'rene/dataset/G0.5F0.1/dos/1=0.001_1' ; ...
    'rene/dataset/G0.5F0.5/dos/1=0.001_1' ; ...
    'rene/dataset/G0.5F0.5/uno/1=0.001_1' ; ...
    'rene/dataset/G0.5F0.5/uno/1=0.001_2' ; ...
    'rene/dataset/G0.5F0.5/uno/1=0.001_3' ; ...
    'rene/dataset/G0.5F1/T3dos/1=0.001_1' ; ...
    'rene/dataset/G0.5F1/dos/1=0.001_1' ; ...
    'rene/dataset/G0.5F1/tres/1=0.001_1' ; ...
    'rene/dataset/G0.5F1/uno/1=0.001_1' ; ...
    'rene/dataset/G0.5F1/uno/1=0.001_2' ; ...
    'rene/dataset/G0.5F1/uno/1=0.001_3' ; ...
    'rene/dataset/G0.75F0.1/uno/1=0.001_1' ; ...
    'rene/dataset/G0.75F0.5/uno/1=0.001_1' ; ...
    'rene/dataset/G0.75F1/uno/1=0.001_1' ; ...
    'rene/dataset/G0.7F0.5/uno/1=0.001_1' ; ...
    'rene/dataset/G0.7F0.5/uno/1=0.001_2' ; ...
    'rene/dataset/G0.7F0.5/uno/1=0.001_3' ; ...
    'rene/dataset/G0.7F1/uno/1=0.001_1' ; ...
    'rene/dataset/G0.7F1/uno/1=0.001_2' ; ...
    'rene/dataset/G0.8F0.1/dos/1=0.001_1' ; ...
    'rene/dataset/G0.8F0.1/tres/1=0.001_1' ; ...
    'rene/dataset/G0.8F0.1/uno/1=0.001_1' ; ...
    'rene/dataset/G0.8F0.5/uno/1=0.001_1' ; ...
    'rene/dataset/G0.8F1/T3uno/1=0.001_1' ; ...
    'rene/dataset/G1F0.1/uno/1=0.001_1' ; ...
    'rene/dataset/G1F0.5/uno/1=0.001_1' ; ...
    'rene/dataset/G1F0.5/uno/1=0.001_2' ; ...
    'rene/dataset/G1F0.5/uno/1=0.001_3' ; ...
    'rene/dataset/G1F1/T3uno/1=0.001_1' ; ...
    'rene/dataset/G1F1/uno/1=0.001_1' ; ...
    'rene/dataset/G1F1/uno/1=0.001_2' ; ...
    'rene/dataset/G1F1/uno/1=0.001_3' ; ...
    'rene/dataset/G0.75F1/T3uno/1=0.001_1';
    };
    case 6.1
  dirs = { ...
    'rene/dataset/G0.5F0.1/dos/1=0.001_1' ; ...
    'rene/dataset/G0.5F0.5/dos/1=0.001_1' ; ...
    'rene/dataset/G0.5F0.5/uno/1=0.001_1' ; ...
    'rene/dataset/G0.5F0.5/uno/1=0.001_2' ; ...
    'rene/dataset/G0.5F0.5/uno/1=0.001_3' ; ...
    'rene/dataset/G0.5F1/T3dos/1=0.001_1' ; ...
    'rene/dataset/G0.5F1/T3dosbis/1=0.001_1' ; ...
    'rene/dataset/G0.5F1/dos/1=0.001_1' ; ...
    'rene/dataset/G0.5F1/tres/1=0.001_1' ; ...
    'rene/dataset/G0.5F1/uno/1=0.001_1' ; ...
    'rene/dataset/G0.5F1/uno/1=0.001_2' ; ...
    'rene/dataset/G0.5F1/uno/1=0.001_3' ; ...
    'rene/dataset/G0.75F0.1/dos/1=0.001_1' ; ...
    'rene/dataset/G0.75F0.1/uno/1=0.001_1' ; ...
    'rene/dataset/G0.75F0.5/uno/1=0.001_1' ; ...
    'rene/dataset/G0.75F1/T3uno/1=0.001_1' ; ...
    'rene/dataset/G0.75F1/uno/1=0.001_1' ; ...
    'rene/dataset/G0.7F0.5/uno/1=0.001_1' ; ...
    'rene/dataset/G0.7F0.5/uno/1=0.001_2' ; ...
    'rene/dataset/G0.7F0.5/uno/1=0.001_3' ; ...
    'rene/dataset/G0.7F1/uno/1=0.001_1' ; ...
    'rene/dataset/G0.7F1/uno/1=0.001_2' ; ...
    'rene/dataset/G0.8F0.1/dos/1=0.001_1' ; ...
    'rene/dataset/G0.8F0.1/tres/1=0.001_1' ; ...
    'rene/dataset/G0.8F0.1/uno/1=0.001_1' ; ...
    'rene/dataset/G0.8F0.5/uno/1=0.001_1' ; ...
    'rene/dataset/G0.8F1/T3uno/1=0.001_1' ; ...
    'rene/dataset/G0.8F1/uno/1=0.001_1' ; ...
    'rene/dataset/G1F0.1/uno/1=0.001_1' ; ...
    'rene/dataset/G1F0.5/uno/1=0.001_1' ; ...
    'rene/dataset/G1F0.5/uno/1=0.001_2' ; ...
    'rene/dataset/G1F0.5/uno/1=0.001_3' ; ...
    'rene/dataset/G1F1/T3uno/1=0.001_1' ; ...
    'rene/dataset/G1F1/T3unobis/1=0.001_1' ; ...
    'rene/dataset/G1F1/uno/1=0.001_1' ; ...
    'rene/dataset/G1F1/uno/1=0.001_2' ; ...
    'rene/dataset/G1F1/uno/1=0.001_3' ; ...
    };
    case 6.5
      dirs = {...
    'rene/dataset/G1F1/T3unobis/1=0.001_1';
    'rene/dataset/G0.5F1/T3dosbis/1=0.001_1';
    'rene/dataset/G0.8F1/uno/1=0.001_1' ; ...
    'rene/dataset/G0.75F0.1/dos/1=0.001_1' ; ...
    };
    case 7
      dirs = {...
        'genomas\porhacer\G1F1\2009_Jul_29_14_22_12\5\1=0.001_1'; ...
        'genomas\porhacer\G0.8F1\2009_Jul_29_12_21_07\2\1=0.001_1'; ...
        'genomas\porhacer\G0.5F1\2009_Jul_29_12_16_39\3\1=0.001_1'; ...
        };
    case 8
      dirs = {...
        'rene/src/2009_Jul_29_12_16_39/3/1=0.001_1'; ...
        'rene/src/2009_Jul_29_12_21_07/2/1=0.001_1'; ...
        'rene/src/2009_Jul_29_14_22_12/5/1=0.001_1'; ...
        };
    case 9
      dirs = {...
        'rene/dataset/G0.8F1/T3uno/1=0.001_1'; ...
        'rene/dataset/G1F1/T3uno/1=0.001_1'; ...
        };
    case 10
      dirs = {...
        'rene/dataset/G1F1/uno/1=0.001_1'; ...
        'rene/dataset/G0.75F1/uno/1=0.001_1'; ...
        'rene/dataset/G0.8F1/uno/1=0.001_1'; ...
        };
    case 11
      dirs = {...
        'rene/dataset/G1F1/uno/1=0.001_1'; ...
        'rene/dataset/G0.75F1/uno/1=0.001_1'; ...
        'rene/dataset/G0.5F1/tres/1=0.001_1'; ...
        };
    case 12
      dirs = {...
        'rene/dataset/G1F1/uno/1=0.001_1'; ...
        'rene/dataset/G0.75F1/uno/1=0.001_1'; ...
        'rene/dataset/G0.5F1/dos/1=0.001_1'; ...
        };
    case 13
      dirs = {...
        'rene/dataset/G1F1/uno/1=0.001_3'; ...
        'rene/dataset/G0.75F1/uno/1=0.001_1'; ...
        'rene/dataset/G0.5F1/dos/1=0.001_1'; ...
        };
  end
end


switch whattodo
  case 0
    %make clusterings for all these simulations
    for k=1:numel(dirs)
      makeAllSpeciations(dirs{k}, 1, [], 'rene/src', true, [2 3 4]);
    end
  case 0.1
    %make clusterings for all these simulations
    for k=1:numel(dirs)
      makeAllSpeciations(dirs{k}, 1, [], 'rene/src', true, [2]);
    end
  case 1
    %make clusterings for all these simulations
    for k=1:numel(dirs)
      makeAllSpeciations(dirs{k}, 1, [], 'rene/src', true);
    end
  case 2
    %test if the clusterings are finished for all these simulations
    prf={'NM_';'E1_';'E2_';'MP_'};
    for k=1:numel(dirs)
      faltan = false(size(prf));
      for z=1:numel(prf)
        nm = [dirs{k} filesep prf{z} 'numGroups.mat'];
        nm2 = [dirs{k} filesep prf{z} 'individualsSpeciesJust.mat'];
        faltan(z) = not(exist(nm, 'file')) || not(exist(nm2, 'file'));
      end
      if any(faltan)
        fprintf('NOT FINISHED FOR %s, <%s> REMAINS TO BE DONE!!\n', dirs{k}, any2str(prf(faltan)));
      else
        fprintf('FINISHED FOR %s!!!\n', dirs{k});
      end
    end
  case 3
    %make the "figuras" struct
    if nargout<1
      error('We need you to provide at least one output argument!!!');
    end
    varargout{1} = makeFigurasStructDOIT(dirs);
  case 4
    %make the minimized genome study
    for k=1:numel(dirs)
      nm = [dirs{k} filesep 'reducedGenomeStudy.mat'];
      if not(exist(nm, 'file'))
        fprintf('GOING FOR %s\n', nm);
        res = analyzeGenomeReduced(dirs{k}, [], false); %#ok<NASGU>
        save(nm, 'res');
        clear res;
      end
    end
  case 5
    %remake the broken matrixDist files, this time making sure the fucking
    %var is below 2GB (even if the compressed file is less than 50MB, just
    %crazy, mathworks guys suck!!!!
    data = {'nummuts'  'NM_';...
            'edit'     'E1_';...
            'miniedit' 'E2_';...
            'morpho'   'MP_'};
    modes = data(:,1);
    prefixes = data(:,2);
    parallel = true;
    options = struct('justMatrixDist', true, 'alsoClasificados', false);
    for k=1:numel(dirs)
      d = dir(dirs{k});
      dn = {d.name}';
      genf = [];
      for z=1:numel(prefixes)
        fl = [prefixes{z} 'matrixDist.mat'];
        fprintf('looking for bad %s...\n', [dirs{k} filesep fl]);
        idx = find(strcmp(fl, dn));
        if not(isempty(idx)) && d(idx).bytes==128
          if isempty(genf)
            fprintf('finding out the number of generations...\n');
            genf = lastGenMAT(dirs{k});
          end
          fprintf('speciation: %s, clustering\n', modes{z});
          HierarchicalClusteringALL(dirs{k},[],modes{z},1,genf,parallel,prefixes{z}, 'rene/src', options);
          compressGemaFiles(dirs{k}, 1, genf, prefixes{z}, 4);
          makeAllSpeciations(dirs{k}, 1, [], 'rene/src', true);
        end
      end
    end
  case 5.2
    %remake the broken matrixDist files, this time making sure the fucking
    %var is below 2GB (even if the compressed file is less than 50MB, just
    %crazy, mathworks guys suck!!!!
    data = {...%'nummuts'  'NM_';...
            'edit'     'E1_';...
            'miniedit' 'E2_';...
            'morpho'   'MP_'};
    modes = data(:,1);
    prefixes = data(:,2);
    parallel = true;
    options = struct('justMatrixDist', true, 'alsoClasificados', true);
    for k=1:numel(dirs)
      d = dir(dirs{k});
      dn = {d.name}';
      genf = [];
      for z=1:numel(prefixes)
        fl = [prefixes{z} 'matrixDist.mat'];
        fprintf('looking for %s...\n', [dirs{k} filesep fl]);
        if not(exist([dirs{k} filesep fl], 'file'))
          fprintf('IT DOES NO EXIST!!!\n');
          if isempty(genf)
            fprintf('finding out the number of generations...\n');
            genf = lastGenMAT(dirs{k});
          end
          fprintf('matrixDist: %s, clustering\n', modes{z});
          HierarchicalClusteringALL(dirs{k},[],modes{z},1,genf,parallel,prefixes{z}, 'rene/src', options);
          compressGemaFiles(dirs{k}, 1, genf, prefixes{z}, 4);
          makeAllSpeciations(dirs{k}, 1, [], 'rene/src', true);
        end
      end
    end
  case 6
    %make phylotypic speciations
    data = {...
            'morpho'   'MP_' 'PMP_'; ...
            'miniedit' 'E2_' 'PE2_';...
            'edit'     'E1_' 'PE1_';...
            };
    umbrales = [10 1];
    modes = data(:,1);
    prefixesALL = data(:,2);
    prefixes = data(:,3);
    for k=1:numel(dirs)
      nm = [dirs{k} filesep 'poph.mat'];
      if exist(nm, 'file')
        poph = load(nm);
        poph = poph.poph;
      else
        poph = loadSimulation(dirs{k});
      end
      for z=1:numel(modes)
        for q=1:numel(umbrales)
          fprintf('Processing mode %s with threshold %d for %s\n', modes{z}, umbrales(q), dirs{k});
          phylotypicSpeciation(dirs{k}, poph, modes{z}, umbrales(q), prefixesALL{z}, prefixes{z});
        end
      end
      clear poph;
    end
  case 7
    %make the "figuras" struct for phylotypic speciations
    if nargout<2
      error('We need you to provide two output arguments!!!');
    end
    umbral = 1;
    varargout{1} = makeFigurasStructDOITOTHER(dirs, 'phylotipic', umbral);
    umbral = 10;
    varargout{2} = makeFigurasStructDOITOTHER(dirs, 'phylotipic', umbral);
  case 7.5
    %make the "figuras" struct for phylotypic speciations
    if nargout<1
      error('We need you to provide one output arguments!!!');
    end
    varargout{1} = makeFigurasStructDOITOTHER(dirs, 'raw');
  case 8
    %calculate the complexity
    for k=1:numel(dirs)
      generateSimulationCompressions(dirs{k}, [], [], 9);
    end
  case {9, 9.1}
    %options 12, 12.1 in compressGemaFiles
    opt = whattodo+3 ;
    %call compressGemaFiles for old simulations
    for k=1:numel(dirs)
      command =  ['awk ''END {print $1}'' ' dirs{k} filesep 'poblacion.txt'];
      [res res]= system(command);
      ngens = str2double(res);
      fprintf('TO %d in %s\n', ngens, dirs{k});
      compressGemaFiles(dirs{k}, 1, ngens, [], opt);
    end
  case 10
    %remake the clustering, but with Rene's criterion (the prototype's
    %dispersion is measured as the mean distance to other members of the
    %cluster, not the absolute sum of distances).
    error('to be elaborated...');
    numberOfJobs = [];
    path = 'rene/src';
    recalculateClustersScript(numberOfJobs, path, dirs, indexOpts, indexType);
end

function lastGen = lastGenMAT(basedir)
poph = load([basedir filesep 'poph.mat']);
lastGen = max(poph.poph.generation);
clear poph;

function figuras = makeFigurasStructDOIT(basedirs, colors, indexes)

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
          14, 'numClusters_mutations'; ...
          };
clusterings = {'morphological', 'MP_' 11; ...
               'edit distance', 'E1_' 12; ...
               'reduced edit distance', 'E2_' 13; ...
               'mutation distannce', 'NM_' 14; ...
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
      end
    end
end

figuras = struct('campos', {campos}, ...
                 'tipos', {tabla(4,:)}, ...
                 'tabla', {tabla});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function figuras = makeFigurasStructDOITOTHER(basedirs, mode, umbral, colors, indexes)

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
clusterings = {'morphological', 'PMP_' 11; ...
               'edit distance', 'PE1_' 12; ...
               'reduced edit distance', 'PE2_' 13; ...
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
    switch mode
      case 'phylotipic'
        for z=1:numel(prefixes)
          fprintf('   numSpecies %s...\n', clusterings{z});
          gemafile = [b filesep mat2str(umbral) prefixes{z} 'phyloNumSpecies.mat'];
          if exist(gemafile, 'file')
            numSpecies = load(gemafile);
            numSpecies = numSpecies.numSpecies;
            numSpecies = cellfun(@(x)sum(x>0), numSpecies(:));
            tabla{prefixesi(z),k} = numSpecies(2:end);
          end
        end
      case 'raw'
          fprintf('   raw diversity...\n');
          tabla{prefixesi(1),k} = zeros(mxg,1);
          tabla{prefixesi(2),k} = zeros(mxg,1);
          tabla{prefixesi(3),k} = zeros(mxg,1);
          for t=1:mxg
            np = sprintf('P%03g', t);
            P = load([b filesep np '.mat']);
            P = P.(np);
              ramasok = find(P.counts(:,3)==0);
              raster = P.raster(ramasok);
              genome = P.genome(ramasok);
              ming = P.mingnm(ramasok);
             rmp = true(size(ramasok));
             for za=1:numel(rmp)-1;
               zz=raster{za};
               for zb=za+1:numel(rmp)
                 if isequal(zz, raster{zb})
                   rmp(za) = false;
                   break;
                 end
               end
             end
             tabla{prefixesi(1),k}(t) = sum((rmp));
             tabla{prefixesi(2),k}(t) = numel(unique(genome));
             tabla{prefixesi(3),k}(t) = numel(unique(ming));
          end
    end
end

figuras = struct('campos', {campos}, ...
                 'tipos', {tabla(4,:)}, ...
                 'tabla', {tabla});
               
