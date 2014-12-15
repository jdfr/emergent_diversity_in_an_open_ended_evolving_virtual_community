function HierarchicalClusteringALL(basedir,poph,mode,GenIn,GenFin,parallel,prefix, path,options)
%HierarchicalClusteringALL('bestia\variable_pop_retune\afinando1\T2_A20\continua_guai_2009_Jun_24_12_31_04\1=0.0001_1', poph, 'morpho', 440, 440, false, '', '');
if not(exist('parallel', 'var'))
  parallel = false;
end
if not(exist('prefix', 'var'))
  prefix = '';
end
if not(exist('path', 'var'))
  path = pwd;
end
if not(exist('options', 'var'))
  options = struct;
end

if islogical(options)
  options = struct('justMatrixDist', options);
end

opts = {'justMatrixDist', false; 'alsoClasificados', false};
for k=1:size(opts,1)
  if isfield(options, opts{k,1})
    options.(opts{k,1}) = opts{k,2};
  end
end

if isempty(GenFin)
  generations = GenIn;
else
  generations = (GenIn:GenFin)';
end

% poph = loadSimulation(basedir);

if parallel
  jobArgs={'Tag', 'clusteringEDITDISTANCE', 'PathDependencies',  {path}};
  generations = onlyNonProcessed(basedir, prefix, generations);
  if isempty(generations)
    fprintf('Clustering mode %s already done for %s!!!\n', mode, basedir);
    return
  else
    fprintf('NUMGENERATIONS: %d\n', numel(generations));
  end
  gens = array2cell(generations);
  pophs = cell(size(gens));
  basedirs = repmat({basedir}, size(gens));
  prefixes = repmat({prefix}, size(gens));
  modes = repmat({mode}, size(gens));
  optionss = repmat({options}, size(gens));
  fprintf ('GOING PARALLEL (%d JOBS)\n', numel(modes));
  if numel(generations)<10
    fprintf('only these gens: %s\n', any2str(generations));
  end
  dummy = my_dfeval(@doAGeneration, modes, pophs, basedirs, gens, prefixes, optionss, jobArgs{:});
  fprintf('RECEIVING PARALLEL\n');
else
  if isempty(poph)
     poph = loadSimulation(basedir);
  end
  switch mode
    case 'nummuts'
      zs=load([basedir filesep 'mutdists.mat']);
      poph.tree.numMuts=zs.numMuts;
      poph.tree.mutDists=zs.mutDists;
      clear zs;
  end
  for a = 1:numel(generations)
    doAGeneration(mode, poph, basedir, generations(a), prefix, options);
  end
end

close all;

function  gens = onlyNonProcessed(basedir, prefix, gens)
calc = false(size(gens));
pr = [basedir filesep prefix];
for k=1:numel(gens)
    n = num2str(gens(k), '%03g');
    calc(k) = not( exist([pr 'Z' n '.mat'], 'file') && exist([pr 'matrixDist' n '.mat'], 'file') && exist([pr 'indicesClasificados' n '.mat'], 'file') );
end
gens = gens(calc);

function dummy = doAGeneration(mode, poph, basedir, a, prefix, options)
dummy = a;
    tosave1 = ['Z' num2str(a, '%03g')];
    tosave2 = ['matrixDist' num2str(a, '%03g')];
    tosave3 = ['indicesClasificados' num2str(a, '%03g')];
    ex1 = exist([basedir filesep prefix tosave1 '.mat'], 'file');
    ex2 = exist([basedir filesep prefix tosave2 '.mat'], 'file');
    ex3 = exist([basedir filesep prefix tosave3 '.mat'], 'file');
    if options.justMatrixDist
      ex = ex2;
      if options.alsoClasificados
        ex = ex && ex3;
      end
    else
      ex = ex1 && ex2 && ex3;
    end
    if ex
      fprintf('GENERATION %d ALREADY DONE!!\n', a);
      return;
    end
%     numIndividuals = zeros(1,a+1);
%     for i = 1 : a+1
%         numIndividuals(i) = sum(poph.generation == i-1); % numIndividuals(1) = numero de individuos en la generacion 0
%     end
% 
if isempty(poph)
  if exist([basedir filesep 'poph.mat'], 'file')
    poph = load([basedir filesep 'poph.mat']);
    poph = poph.poph;
  else
    poph = loadSimulation(basedir);
  end
  switch mode
    case 'nummuts'
      zs=load([basedir filesep 'mutdists.mat']);
      poph.tree.numMuts=zs.numMuts;
      poph.tree.mutDists=zs.mutDists;
      clear zs;
  end
end
a
thisGen = find(poph.generation==a);
numIndividual = numel(thisGen);

nPixelsInSoil = poph.nPixelsInSoil;
    ramasNegativas = nPixelsInSoil(thisGen);%sum(numIndividuals(1:a))+1:sum(numIndividuals(1:a+1)));

    nonEmptyNeg = true; %no clasifico vacios y negativos

    switch mode
      case 'miniedit'
        genome = poph.minGenome(thisGen);
      case 'edit'
        genome = poph.genome(thisGen);
      case 'morpho'
        load(sprintf('%s/P%03d.mat',basedir,a));
        P = eval(sprintf('P%03d', a));
        raster = P.raster;
        dimensions = P.dimensions;
        offset = P.offset;
    end
    % las generaciones empiezan en 0, tal y como aparecen en poblacion.txt
    indicesClasificados = [];

    if nonEmptyNeg

        negativos = find(ramasNegativas > 0);

        indicesNoClasificados = negativos;
        indicesClasificados = setdiff(1:numIndividual,indicesNoClasificados);
        numIndividual = length(indicesClasificados);
        switch mode
          case {'edit', 'miniedit'}
            genome(indicesNoClasificados,:)=[];
          case 'morpho'
            raster(indicesNoClasificados,:)=[];
            dimensions(indicesNoClasificados,:)=[];
            offset(indicesNoClasificados,:)=[];
        end
    end

    dist = zeros(1,(numIndividual*(numIndividual-1))/2);
    k = 1;
    if numIndividual > 1
      switch mode
        case {'edit', 'miniedit'}
          for i = 1 : numIndividual-1
              if mod(i,100)==0
                fprintf('%s %d/%d\n', mode, i,numIndividual);
              end
              for j = i+1 : numIndividual
                  dist(k) = strdstMX(genome{i}, genome{j});%strdist(genome{i}, genome{j});
                  k = k + 1;
              end
          end
        case 'morpho'
          for i = 1 : numIndividual-1
              if mod(i,100)==0
                fprintf('%s %d/%d\n', mode, i,numIndividual);
              end
              raster2 = raster{i};
              offset2 = offset(i,:);
              for j = i+1 : numIndividual
                  raster1 = raster{j};
                  offset1 = offset(j,:);
                  dist(k) = 1 - similitud2(raster1,offset1,raster2,offset2);
                  k = k + 1;
              end
          end
        case 'nummuts'
          for i = 1 : numIndividual-1
              if mod(i,100)==0
                fprintf('%s %d/%d\n', mode, i,numIndividual);
              end
              for j = i+1 : numIndividual
                  dist(k) = mutDist(poph, a, indicesClasificados(i), indicesClasificados(j));
                  k = k + 1;
              end
          end
      end
        matrixDist = squareform(dist);
        if options.justMatrixDist
          Z=[];
        else
          Z = linkage(dist);
        end
        % fig = figure;%('visible', 'off');
        %[H,T,perm]=dendrogram(Z,0);
        %set(H,'ButtonDownFcn', {@showInf,raster,dimensions,offset,H,perm,indicesClasificados}, 'Interruptible', 'off');
        % set(gca,'XTick',[]);
    else
        Z = [];
        matrixDist = [];
    end
    % saveas(fig,[basedir,'hierarchicalClustering' num2str(Generacion, '%03g')],'png');
    % close(fig);

    eval([tosave1 ' = Z;']);
    eval([tosave2 ' = matrixDist;']);
    eval([tosave3 ' = indicesClasificados;']);
    
    if options.justMatrixDist
      save([basedir filesep prefix tosave2 '.mat'],tosave2);
      if options.alsoClasificados
        save([basedir filesep prefix tosave3 '.mat'],tosave3);
      end
    else
      save([basedir filesep prefix tosave1 '.mat'],tosave1);
      save([basedir filesep prefix tosave2 '.mat'],tosave2);
      save([basedir filesep prefix tosave3 '.mat'],tosave3);
    end

    clear('Z', 'matrixDist','indicesClasificados', tosave1,tosave2,tosave3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sim = similitud2(raster1,offset1,raster2,offset2)%, show)
    pixel1 = [raster1{:}];
    pixel2 = [raster2{:}];
    v = isempty(pixel1) + isempty(pixel2);
if v==2
    sim = 1;
elseif v==1
    sim = 0;
else
    
    if size(pixel1,1) == 1 
        pixel1 = [];

        pixel1(:,1) = raster1{1};
        pixel1(:,2) = raster1{2};
    end
    if size(pixel2,1) == 1 
        pixel2 = [];

        pixel2(:,1) = raster2{1};
        pixel2(:,2) = raster2{2};
    end
    
    pixel1(:,1) = pixel1(:,1) - offset1(1);
    pixel1(:,2) = pixel1(:,2) - offset1(2);

    pixel2(:,1) = pixel2(:,1) - offset2(1);
    pixel2(:,2) = pixel2(:,2) - offset2(2);
    
    sim = similitudFast(pixel1, pixel2);
    %sim = size(intersect(pixel1,pixel2,'rows'),1)/size(union(pixel1,pixel2,'rows'),1);
end

% if isempty([raster1{:}]) && isempty([raster2{:}])
%     sim = 1;
% elseif (isempty([raster1{:}]) && ~isempty([raster2{:}])) || (~isempty([raster1{:}]) && isempty([raster2{:}]))
%     sim = 0;
% else
%     pixel1 = [raster1{:}];
%     pixel2 = [raster2{:}];
% 
%     pixel1(:,1) = pixel1(:,1) - offset1(1);
%     pixel1(:,2) = pixel1(:,2) - offset1(2);
% 
%     pixel2(:,1) = pixel2(:,1) - offset2(1);
%     pixel2(:,2) = pixel2(:,2) - offset2(2);
% 
%     %     sraster1 = size(pixel1,1);
%     %     sraster2 = size(pixel2,1);
% 
%     sim = size(intersect(pixel1,pixel2,'rows'),1)/size(union(pixel1,pixel2,'rows'),1);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=mutDist(poph, gen, idx1, idx2)
d=poph.tree.mutDists{gen+1}(idx1, idx2);
