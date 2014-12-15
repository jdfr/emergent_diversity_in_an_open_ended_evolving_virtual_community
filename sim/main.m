function ACIParams = main(parmode, mode, varargin)

MAINSUB = 'resultados';

if ~exist('parmode', 'var')
  parmode = 'embedded';
end

if iscell(mode)
    MAINSUB = mode{1};
    mode    = mode{2};
end

MAINSUBD = MAINSUB;
if ~isempty(MAINSUBD)
  MAINSUBD = [MAINSUBD filesep];
end

switch mode
  case 'replay' % additionalArgs = {subdirectory,newsubdirectory}
    subdirectory    = varargin{1};
    newsubdirectory = varargin{2};
    if iscell(subdirectory)
      subdirectory = subdirectory{1}; 
      estadofile        = [subdirectory,filesep,'estado.mat'];
    else
      if isempty(MAINSUBD)
        estadofile      = [subdirectory,filesep,'estado.mat'];
      else
        estadofile      = [MAINSUBD,filesep , subdirectory,filesep,'estado.mat'];
      end
    end
    load(estadofile,'-mat');      % carga los parámetros del fichero estado
    ACIParams = putDefaultParameters(ACIParams); %#ok<NODEF>
    ACIParams = modifyParams(ACIParams, varargin(3:end)); %#ok<NODEF>
    rand(ACIParams.randmethod, ACIParams.seedRandom); %#ok<NODEF>
    if isempty(MAINSUBD)
      ACIParams.nomdir = newsubdirectory;
    else
      ACIParams.nomdir = [MAINSUBD filesep newsubdirectory];
    end
    if exist(ACIParams.nomdir, 'dir')
        rmdir(ACIParams.nomdir, 's');
    end
    mkdir(ACIParams.nomdir);
    save([ACIParams.nomdir,filesep,'estado.mat'], 'ACIParams','-mat');
  case 'continue' % additionalArgs = {subsubdirectory,newsubdirectory,presion,numgenerations}
    subsubdirectory = varargin{1};
    subdirectory    = [subsubdirectory filesep '..'];
    newsubdirectory = varargin{2};
    presion         = varargin{3};
    numgenerations  = varargin{4};
    estadofile      = [MAINSUBD,subdirectory,filesep,'estado.mat'];
    load(estadofile,'-mat');      % carga los parámetros del fichero estado
    ACIParams.nomdir = [MAINSUBD newsubdirectory];
    if exist(ACIParams.nomdir, 'dir')
        rmdir(ACIParams.nomdir, 's');
    end
    mkdir(ACIParams.nomdir);
    poph                          = loadSimulation([MAINSUBD,subsubdirectory], false, false, false);
    lastGen                       = max(poph.generation);
    thisGen                       = poph.generation==lastGen;
    initialPopulation             = poph.genome(thisGen);
    initialPos                    = poph.xpos(thisGen);
    clear poph thisGen;
    ACIParams.initialGeneration   = lastGen;
    ACIParams.randomPlacementForInitialGeneration = false;
    ACIParams.generaciones        = numgenerations;
    ACIParams.PRESIONES           = presion;
    ACIParams.num_iteraciones     = numel(presion);
    ACIParams.seedRandom          = makeRandomSeed;
    ACIParams.randmethod          = 'twister';
    save([ACIParams.nomdir,filesep,'estado.mat'], 'ACIParams','-mat');
    ACIParams.initialPopulation   = initialPopulation;
    ACIParams.initialPos          = initialPos;
    ACIParams.tam_poblacion       = numel(initialPopulation);
    ACIParams.terclus.ranges      = calculateRanges(ACIParams.tam_poblacion, ACIParams.terclus.numRanges);
    clear initialPopulation initialPos;
    ACIParams = putDefaultParameters(ACIParams); %#ok<NODEF>
    ACIParams = modifyParams(ACIParams, varargin(5:end));
    rand(ACIParams.randmethod, ACIParams.seedRandom); %#ok<NODEF>
  case {'newsim', 'newsimSeeded'} % additionalArgs = pairs of strings 'option', 'value'
    switch mode
      case 'newsim'
        otherArgs = varargin;
      case 'newsimSeeded'
        seed = varargin{1};
        otherArgs = varargin(2:end);
    end
    ACIParams                     = struct;
    ACIParams.seedRandom          = makeRandomSeed;
    ACIParams.randmethod          = 'twister';
%---------------DEFINE PARAMETERS-----------------------------
    % PARAMETROS GENERALES
    ACIParams.PRESIONES           = 0.0001;       %[1 0.2 0.1 0.0001 0.5];% [0.75 0.5 1/3 0.2 0.1 0.0001];
    ACIParams.num_iteraciones     = numel(ACIParams.PRESIONES);
    ACIParams.ensayos             = 10; %nº de veces que se hace el experimento para cada presion
    ACIParams.generaciones        = 500;%1000;%500;         % número de intentos para encontrar un árbol de una relación de aspecto dada
    ACIParams.tam_poblacion       = 1;1000;%1000;%300;         % número de individuos en la población
    %this flag tells whether to use a variable population size
    ACIParams.variablePopSize     = true;
    %if variablePopSize==true, this variable determines the way to
    %calculate the number of sons:
    %      -mixed: the fitness is powered to alphaF, and the result is
    %              divided in integer and fractional parts. The integer is
    %              the number of sons, and the fraction is the probability
    %              to have an additional son
    %      -round, ceil, floor: the fitness is rounded, ceiled or floored,
    %                           and the result is the amount of sons
    %      -round#INTERVAL, ceil#INTERVAL, floor#INTERVAL: the fitness is
    %                           rounded, ceiled or floored, then perturbed
    %                           by a random number uniformly distribute in
    %                           INTERVAL, and the result is the amount of
    %                           sons
    ACIParams.variableNumSonsMode = 'round#[0 0.6]'; %'mixed'; 'round'; 'ceil'; 'floor';
    %if true, and there is just one tree, it is pardoned from death if it
    %spawns no descendants. Useful to no abort simulations too early
    ACIParams.preserveUniqueTree  = false; 
    %if variablePopSize==true and variableNumSonsMode=='mixed', this is the
    %minimal probability to leave a descendant in the next population, even
    %if the tree has zero fitness 
    ACIParams.minimalDescendProbability = 0.5; 
    %if variablePopSize==true, this is the exponent to power the fitness
    %before calculating the number of sons
    ACIParams.alphaF              = 0.1;
    ACIParams.ramasobj            = 4;         % número deseado de ramas en el árbol
    ACIParams.presion             = []; %0.1;        % presion selectiva (0: no ha presión, 1: escalado normal)
    ACIParams.nomdir              = directorio(MAINSUB);

    % PARAMETROS DE ACI
    %modes for insertion operator:
    %     -'one_level': the operator can insert '+', '-' and 'G'
    %                   everywhere, but '[]' only at root level, to prevent
    %                   nesting
    %     -'normal': the operator can insert '+', '-', 'G', and '[]'
    %                everywhere
    %     -insert_list: the operator can insert any string in he cell list
    %                   "insert_list" everywhere 
    %false, it cannot insert nested brackets
    ACIParams.mode_op_insercion  = 'normal'; %'normal'; 'one_level'; {{'+', '-'}};
    %if true, only insert symbols inside brackets
    ACIParams.insertOnlyNested   = false; %true; %false;
    %if true, the deletion operator can delete non-empty brackets. If
    %false, it can only delete empty brackets
    ACIParams.op_eliminacion_extended  = true;%false;
    %the way to select a symbol for mutation: by symbol (first select a
    %symbol, then select a position in the genome that has that symbol) or
    %by position (first select a position in the genome, then delete it)
    ACIParams.selectSymbolForMutation  = 'byPos'; %'byPos'; %'bySymbol';
    %the mutations will be applied in normal or random order, or only just
    %one of them
    ACIParams.mutationOrder            = 'justone';%'normal'; %'random'; %'justone';
    ACIParams.allowPrefixChangeDir     = false;
   
    ACIParams.maxLongs.maxInsert  = []; %this is the max long of sub-string for mutations that insert sub-strings
    ACIParams.maxLongs.maxString  = []; %this is the max long of sub-string for mutations that insert sub-strings
    ACIParams.maxLongs.maxGs      = []; %this is the max number of G's in a string
    ACIParams.maxLongs.penalizeFitness = true; %flag to appy penalities based on these values to the fitness function
    ACIParams.maxLongs.penalizeOperators = true; %flag to appy penalities based on these values to the results of mutation operators
    %if the tree has more than this number of negative pixels (pixels under
    %soil), it will not be allowed to leave descend nor be placed on the landscape
    ACIParams.negYThreshold       = 0; %0; %inf; %10;
    ACIParams.long_cadena         = 5;
    ACIParams.n_anidamientos      = 4;
    ACIParams.prob                =  [...
                  0 ...      %      1   elitismo             =  
                0.05...      %      2   op_alteracion       A simbolo  posicion  
                0.05...      %      3   op_dupaleatoria      R segmento posicion
                0.05...      %      4   op_dupnivel          L segmento posicion
                0.05...      %      5   op_dupsecuencia      C segmento posicion
                0.05...      %      6   op_eliminacion       D posicion
                0.05...      %      7   op_insercion         I simbolo  posicion
                0.05...      %      8   op_transferencia     T segmento posicion
                ];   % probabilidad relativa de cada operador
    
    % diferentes tipos de transparencias:
    %       'canonical': todo tapa (tree)
    %       'probabilistic': transparencia probabilistica
    %       'transparentBranches': ramas transparentes (treeRamasTransparentes) 
    %   PARA AÑADIR ALEATORIEDAD EN LA GENERACION DE RAMAS:
    %             ACIParams.T = {'transparentBranches', 'random'}
    ACIParams.T.shadowing         = 'transparentBranches'; %'canonical'; %'probabilistic'; %'transparentBranches';
    ACIParams.T.randomize         = false;
    ACIParams.T.angle             = -5;
    ACIParams.heightT             = [];%1000;%[];%500; %max tree height
    %if heightT is not empty, this flag says what to do with too tall
    %trees: to place them in the landscape (so they will shadow smaller
    %trees, or not placing them (thus, smaller trees have more chances to
    %get light without being shadowed by too tall trees).
    ACIParams.tooTallArePlaced    = false; %true; %false;
    %if any tree becomes taller than this value, the simulation ends
    ACIParams.stopIfTooTall       = [];
    ACIParams.alphaG              = 0.35; % exponente del denominador del fitness
    ACIParams.factorG             = 1;  % factor    del denominador del fitness
    ACIParams.divideFitness       = 'byBranches'; %'byLeafs'; %'byBranches'; %dividir el fitness por el nº de ramas o por el nº de hojas
    ACIParams.N                   = inf; %3000; %ACIParams.tam_poblacion * 15; % anchura del landscape: tamaño de la poblacion * 15.
    ACIParams.initialPos          = [];
    %ACIParams.initialPopulation   = [];
    
%     if ACIParams.N < 3000
%         ACIParams.N = 3000;
%     end
    %way to place the tree inside the landscape
    ACIParams.treePlacement       = 'uniformDistByRoot'; %'uniformDistByRoot' %'random'; %'covering', 'gaussianDist', 'uniformDist'
    %This parameter has several meanings:
    %      -if treePlacement=='gaussianDist', the standard deviation is the
    %       tree width multiplied by this parameter
    %      -if treePlacement=='uniformDist', the maximum width
    ACIParams.treePlacementParam  = 0.2;
    ACIParams.treePlacementOffset = 5;
    
    %the name of the fitness function to be used
    ACIParams.fitnessFunction     = 'fitness2';
    %if fitnessFunction=='fitnessSimilarity', this parameter is the tree
    %that must be matched. This parameter can be:
    %       -a string: it is interpreted as the genome
    %       -a cell array holding as many cells as arguments are returned
    %        by ls2: it is interpreted as the tree itself to be matched
    ACIParams.treeToMatch         = '';
    
    ACIParams.plotting            = true; % guardar 2 pngs por generación con el bosque y los mejores: ocupa mucha memoria!!!
    ACIParams.plottingFactor      = 1;
    ACIParams.plottingSample      = ACIParams.plotting && false; %guardamos el sample de cada generación solo si este flag está activo
    ACIParams.plotSampleMode      = 'firstIndexes'; %'best'; %'firstIndexes'; %tipo de sample: mejores o los primeros de cada generacion
    ACIParams.plotDendroTree      = false; % guardar 2 pngs por generación con el bosque y los mejores: ocupa mucha memoria!!!
    ACIParams.tooLongEvTime       = inf;
    
    % PARAMETROS DE PARALELISMO
    ACIParams.terclus            = struct;
    ACIParams.terclus.doit       = true;
    ACIParams.terclus.tag        = 'iLSystemNEW';
    ACIParams.terclus.jmArgs     = {'LookupURL',         'n0001'};
    ACIParams.terclus.jobArgs    = {'Tag',               ACIParams.terclus.tag, ...
                                      ...'FileDependencies', my_filenamedirrec('*.m'), ...
                                      'PathDependencies', {'srctree'}, ...
                                     };
    ACIParams.terclus.pauseTime         = 1;
    ACIParams.terclus.numRanges         = 1;
    ACIParams.terclus.ranges            = calculateRanges(ACIParams.tam_poblacion, ACIParams.terclus.numRanges);
    ACIParams.terclus.concJobs          = 1;%2; %number of concurrent jobs in the cluster
    ACIParams.terclus.concurrentWorkers = 3;15;%2; %number of concurrent workers in a job
    
    ACIParams = modifyParams(ACIParams, otherArgs);
    
    if strcmp(mode, 'newsimSeeded')
      if ~iscell(seed)
        seed = {seed};
      end
      switch numel(seed)
        case 1
          seed = repmat(seed, ACIParams.tam_poblacion, 1);
        case ACIParams.tam_poblacion
          seed = reshape(seed, [], 1);
        otherwise
          error('If there are more than one seed, the population size and the number of seed must match!!!');
      end
      ACIParams.initialPopulation = seed;
      clear seed otherArgs;
    end

    rand(ACIParams.randmethod, ACIParams.seedRandom);
    fprintf('seedRandom = %20.20g\n', ACIParams.seedRandom);
    
    save([ACIParams.nomdir,filesep,'estado.mat'], 'ACIParams','-mat');
%---------------END DEFINE PARAMETERS-------------------------
    otherwise
        error('mode not understood <%s>', any2str(mode));
end

printReadableEstado([ACIParams.nomdir,filesep,'readableEstado.txt'], ACIParams, 'ACIParams');

%place it here to not save it to estado.mat
switch parmode
    case 'local'
      ACIParams.terclus.jobMgr = findResource('scheduler', 'type', 'local');
    case 'remote'
      ACIParams.terclus.jobMgr = findResource('scheduler', 'type', 'torque', ACIParams.terclus.jmArgs{:});
      %ACIParams.terclus.jobMgr = findResource('scheduler', 'type', 'jobmanager', ACIParams.terclus.jmArgs{:});
    case 'embedded'
      ACIParams.terclus.jobMgr = [];
    otherwise
      error('First argument must be either ''local'' or ''remote'', but it is ''%s''!!!\n', num2str(parmode));
end

if not(isfield(ACIParams, 'plotAtEnd'))
  ACIParams.plotAtEnd = [false false false];
end

ACIParams.completedirs = cell(ACIParams.num_iteraciones*ACIParams.ensayos,1);
zk=1;
for iter=1:ACIParams.num_iteraciones
    ACIParams.num_iteracion= iter;
    ACIParams.presion = ACIParams.PRESIONES(iter);
    fprintf('Iteration nº %g, pressure %.10f\n', iter, ACIParams.presion);
    
    for i = 1:ACIParams.ensayos
        ACIParams.ensayo       = i;
        fprintf('ENSAYO %g\n', i);
        subnomdir = ACI(ACIParams);
        ACIParams.completedirs{zk} = subnomdir;
        if ACIParams.plotAtEnd(1)
          makeAllGraphics(subnomdir, true);
        end
        if false%ACIParams.plotAtEnd(2)
          makeLongGenomeGraphics(subnomdir, inf, true);
        end
        zk=zk+1;
    end
end
if ACIParams.plotAtEnd(3)
  showFitnessPerGenRec(true, true, false, ACIParams.nomdir);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seed = makeRandomSeed
clk  = clock;
seed = mod(floor(sum(clk([1 6 5 4 3 2]).*[1e10 1e8 1e6 1e4 1e2 1])), 2^32);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nomdir, ahora] = directorio(nombre)

if ~isempty(nombre)
  nombre = [nombre filesep];
end

ahora = datestr(now, 'yyyy_mmm_dd_HH_MM_SS');
nomdir = [nombre,ahora];
mkdir(nomdir);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ACIParams = modifyParams(ACIParams, args)
k=1;
while k<numel(args)
  if ~ischar(args{k})
    error('I expect ''option'', ''value'' pairs!!!');
  end
  try
    if (args{k}(1)=='[') && (args{k}(end)==']')
      args{k} = args{k}(2:end-1);
    else
      eval(['ACIParams.' args{k} ';']);
    end
    if ischar(args{k+1})
      eval(['ACIParams.' args{k} ' = ' args{k+1} ';']);
    elseif iscell(args{k+1}) && (numel(args{k+1})==1)
      points = find(args{k}=='.');
      if isempty(points)
        ACIParams.(args{k})                                   = args{k+1}{1};
      else
        starts = [1 (points+1)];
        ends   = [(points-1) numel(args{k})];
        strs   = arrayfun(@(s,e) args{k}(s:e), starts, ends, 'uniformoutput', false);
        switch numel(starts)
          case 2
            ACIParams.(strs{1}).(strs{2})                     = args{k+1}{1};
          case 3
            ACIParams.(strs{1}).(strs{2}).(strs{3})           = args{k+1}{1};
          case 4
            ACIParams.(strs{1}).(strs{2}).(strs{3}).(strs{4}) = args{k+1}{1};
          otherwise
            error('error');
        end
      end
    else
      error('error');
    end
    if any(strcmp(args{k}, {'tam_poblacion', 'terclus.numRanges'}))
      ACIParams.terclus.ranges = calculateRanges(ACIParams.tam_poblacion, ACIParams.terclus.numRanges);
    end
  catch ME
    error('Sorry, but the pair <%s>, <%s> canot be evaluated!!! Parent Error:\n-----------\n%s\n-------------\n', args{k}, any2str(args{k+1}), showError(ME));
  end
  k=k+2;
end

