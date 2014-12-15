function data = morphoOnTheFly(mode, poph, gens, data, fun)
% A function to calculate the phenotypes of each population and apply some
% per-generation function

% poph es la salida de loadsimulation
% gens = [mingen maxgen]

if ischar(poph)
  poph = loadSimulation(poph);
end

if isempty(gens)
    gens = (min(poph.generation):max(poph.generation))';
end

indivStruct = struct('genome', {{[]}}, 'fitness', {nan}, 'randN',{nan}, 'raster', {{[]}}, 'dimensions', {[nan nan]}, 'offset', {[nan nan]}, 'maxiy', {{[]}}, 'isLeaf', {{[]}}, 'counts', {[nan nan nan]}, 'mingnm', {{[]}});

mode(1) = lower(mode(1));
switch mode(1)
  case 'r' %remote
    jobMgr = findResource('scheduler', 'type', 'torque', 'LookupURL', 'n0001');
  case 'l'
    jobMgr = findResource('scheduler', 'type', 'local');
end

numWorkers = 15;

startTime = clock;

Pant = [];

hasTree = isfield(poph, 'tree');

data = fun('init', data, poph, gens, mode);

y=1;
for r=1:numel(gens)
  k = gens(r);
  thisGen = poph.generation==k;
  genomes = poph.genome(thisGen);
  fitness = poph.fitness(thisGen);
  xpos    = poph.xpos(thisGen);
  P       = repmatStruct(indivStruct, [numel(genomes), 1]);
  P.genome(1:end)  = genomes(1:end);
  P.mingnm(1:end)   = poph.minGenome(thisGen);
  P.fitness(1:end) = fitness(1:end);
  P.randN(1:end)   = xpos(1:end);
  clear genomes fitness xpos;
  fprintf('Before computing generation %03d...', k);
  if hasTree && (~isempty(Pant))
    thisGT      = find(poph.tree.gD==k);
    nomutate    = cellfun(@(x) isempty(x) || ((numel(x)>=1) && all(x=='=')), poph.tree.change(thisGT));
    nochange    = nomutate;
    recompute   = not(nochange);
    if any(recompute)
      parents     = poph.tree.iA(thisGT(recompute)); %#ok<FNDSB>
      sons        = poph.tree.iD(thisGT(recompute));
      arethesame  = strcmp(P.mingnm(sons), Pant.mingnm(parents));
      nochange(recompute) = arethesame;
    end
    if any(nochange)
      toCompute   = find(~nochange);
      parents     = poph.tree.iA(thisGT(nochange)); %#ok<FNDSB>
      sons        = poph.tree.iD(thisGT(nochange));
      P.raster(sons,:)     = Pant.raster(parents,:);
      P.isLeaf(sons,:)     = Pant.isLeaf(parents,:);
      P.maxiy(sons,:)      = Pant.maxiy(parents,:);
      P.offset(sons,:)     = Pant.offset(parents,:);
      P.counts(sons,:)     = Pant.counts(parents,:);
      P.dimensions(sons,:) = Pant.dimensions(parents,:);
    else
      toCompute = 1:numel(P.genome);
    end
    allparents = poph.tree.iA(thisGT);
  else
    toCompute = 1:numel(P.genome);
  end
  if ~isempty(toCompute)
    switch mode(1)
      case 'e' %'embedded'
        tic;
        initt = cputime;
        for zz=1:numel(toCompute)
          z = toCompute(zz);
          ttoca = toc;
          passta = cputime;
          [P.raster{z}, P.offset(z,:), P.counts(z,:), P.dimensions(z,:), P.isLeaf{z}, P.maxiy{z}] = ls2(P.genome{z}, poph.ACIParams.T);
          ttoc = toc;
          passt = cputime;
          fprintf('Mira: %d, %d/%d. toc: %07.1f/%07.1f, cputime: %07.1f/%07.1f\n', z, zz, numel(toCompute), ttoc-ttoca, ttoc, passt-passta, passt-initt);
        end
      case {'r', 'l'} %'remote'
        [limitsByW indsByW] = calculateRanges(numel(toCompute), numWorkers);
        nulls               = indsByW==0;
        limitsByW(nulls)    = [];
        indsByW(nulls)      = []; %#ok<NASGU>
        indsByW             = arrayfun(@(x)x, indsByW, 'uniformoutput', false)';
        genomeChunks        = arrayfun(@(s,e)P.genome(toCompute(s:e)), limitsByW(:,1), limitsByW(:,2), 'uniformoutput', false);
        [timeSpent, raster,dimensions,offset,maxiy,leaf,counts] = my_dfeval(jobMgr, @evaluar, genomeChunks, indsByW, repmat({poph.ACIParams.T}, size(genomeChunks)), ...
              'Tag', 'PmatEvaluating', 'PathDependencies', {'srctree'}); %#ok<ASGLU,NASGU>
        P.raster(toCompute,:)     = vertcat(raster{:});
        P.isLeaf(toCompute,:)     = vertcat(leaf{:});
        P.maxiy(toCompute,:)      = vertcat(maxiy{:});
        P.offset(toCompute,:)     = vertcat(offset{:});
        P.counts(toCompute,:)     = vertcat(counts{:});
        P.dimensions(toCompute,:) = vertcat(dimensions{:});
        clear timeSpent raster dimensions offset maxiy leaf counts limitsByW indsByW;
    end
  end
  if hasTree && (~isempty(Pant))
    data = fun('gen', data, poph, gens, mode, thisGen, P, k, y, Pant, nomutate, toCompute, allparents);
  else
    data = fun('gen', data, poph, gens, mode, thisGen, P, k, y, []);
  end
  y=y+1;
  Pant = P;
  clear('P', 'toCompute');
  fprintf(' DONE!!! ELAPSED TIME: %.3f\n', etime(clock, startTime));
end

data = fun('end', data, poph, gens, mode);
