function calculatePmat(mode, poph, gens, destdir)

% mode es remote o embedded
% poph es la salida de loadsimulation
% gens = [mingen maxgen]
% destdir es el directorio donde deben de escribirse los P.mat


indivStruct = struct('genome', {{[]}}, 'fitness', {nan}, 'randN',{nan}, 'raster', {{[]}}, 'dimensions', {[nan nan]}, 'offset', {[nan nan]}, 'maxiy', {{[]}}, 'isLeaf', {{[]}}, 'counts', {[nan nan nan]});

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

for k=gens(1):gens(2)
  thisGen = poph.generation==k;
  genomes = poph.genome(thisGen);
  fitness = poph.fitness(thisGen);
  xpos    = poph.xpos(thisGen);
  clear thisGen;
  P       = repmatStruct(indivStruct, [numel(genomes), 1]);
  P.genome(1:end)  = genomes(1:end);
  P.fitness(1:end) = fitness(1:end);
  P.randN(1:end)   = xpos(1:end);
  clear genomes fitness xpos;
  fprintf('Before computing generation %03d...', k);
  if hasTree && (~isempty(Pant))
    thisGT      = find(poph.tree.gD==k);
    nochange    = cellfun(@(x) isempty(x) || ((numel(x)>=1) && all(x=='=')), poph.tree.change(thisGT));
    thisGT      = thisGT(nochange);
    if ~isempty(thisGT)
      toCompute   = find(~nochange);
      parents     = poph.tree.iA(thisGT); %#ok<FNDSB>
      sons        = poph.tree.iD(thisGT);
      P.raster(sons,:)     = Pant.raster(parents,:);
      P.isLeaf(sons,:)     = Pant.isLeaf(parents,:);
      P.maxiy(sons,:)      = Pant.maxiy(parents,:);
      P.offset(sons,:)     = Pant.offset(parents,:);
      P.counts(sons,:)     = Pant.counts(parents,:);
      P.dimensions(sons,:) = Pant.dimensions(parents,:);
    else
      toCompute = 1:numel(P.genome);
    end
    clear thisGT nochange;
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
  tosave = ['P' num2str(k, '%03g')];
  fprintf(' DONE!!! Saving %s...', tosave);
  eval([tosave ' = P;']);
  save([destdir filesep tosave],tosave);
  Pant = P;
  clear('P', 'toCompute', tosave);
  fprintf(' DONE!!! ELAPSED TIME: %.3f\n', etime(clock, startTime));
end

