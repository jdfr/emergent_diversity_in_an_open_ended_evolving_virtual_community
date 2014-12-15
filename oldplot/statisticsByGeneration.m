function statisticsByGeneration(poph, restrictgens)
%compile and show statistics by generation

idxD = poph.tree.idxD;
idxA = poph.tree.idxA;
gD   = poph.tree.gD;

mingen = min(gD);
maxgen = max(gD);

if nargin>1
  toUse = [];
  switch numel(restrictgens)
    case 1
      if restrictgens<maxgen
        toUse  = gD<=restrictgens;
        maxgen = restrictgens;
      end
    case 2
      if restrictgens(1)>restrictgens(2)
        error('Are you crazy? The second parameter is [startGen endGen], so keep startGen lower or equal than endGen, stupid!');
      end
      if (restrictgens(1)>mingen) || (restrictgens(2)<maxgen)
        toUse  = (gD>=restrictgens(1)) & (gD<=restrictgens(2));
        maxgen = min(restrictgens(2), maxgen);
        mingen = max(restrictgens(1), mingen);
      end
    case 0
      %do nothing
    otherwise
      error('Second argument must have either 1 or 2 elements, not %d!!!\n', numel(restrictgens));
  end
  if ~isempty(toUse)
    idxD   = idxD(toUse);
    idxA   = idxA(toUse);
    gD     = gD(toUse);
  end
end

%make index matrix of parents
if (numel(poph.tree.idxD)==numel(idxD)) && isfield(poph.tree, 'parents')
  parents = poph.tree.parents;
else
  %get parent indexes
%   try
%     sparseIndexes = sparse(1,idxD+1,1:numel(idxD),1, max(idxD)+1, numel(idxD));
%     parents       = full(sparseIndexes(idxA+1)); 
%     clear sparseIndexes;
%   catch ME %#ok<NASGU>
%     if strcmpi(ME.identifier, 'MATLAB:nomem')
      [nevermind, parents] = ismember(idxA, idxD); %#ok<ASGLU>
      parents = parents';
      clear nevermind;
%     else
%       rethrow(ME);
%     end
%   end
end

%sparseIndexes = sparse(1,poph.idx,1:numel(poph.idx));

if (numel(poph.tree.idxD)==numel(idxD)) && isfield(poph.tree, 'sons')
  sons = poph.tree.sons;
else
  sons    = cell(size(parents));
  for z=1:int32(numel(parents))
    pz = parents(z);
    if pz>0
      sons{pz}(end+1,1) = z;
    end
  end
end


ngens = maxgen-mingen+1;

figure;
toUse     = (poph.generation>=mingen) & (poph.generation<=maxgen);
data      = [poph.generation(toUse), poph.fitness(toUse)];
hist3(data, {mingen:maxgen, linspace(min(data(:,2)), max(data(:,2)), 20)});
xlabel('generations'); ylabel('fitnesses'); zlabel('amounts'); grid on;
title('histogram of fitness by generation');

%ngens=ngens;

%fitnesses      = cell(ngens,1);
maxfitness     = zeros(ngens,1);
stdfitness     = zeros(ngens,1);
meanfitness    = zeros(ngens,1);
medianfitness  = zeros(ngens,1);
coefvarfitness = zeros(ngens,1);
%numsons        = cell(ngens,1);
haveonesonratio  = zeros(ngens,1);
havesseveralsonsratio  = zeros(ngens,1);
maxsons        = zeros(ngens,1);
%meansons       = zeros(ngens,1);
mediansons     = zeros(ngens,1);

maxIdx = max(poph.idx);

g = mingen;
for k=1:ngens
  thisGen           = gD==g;
  idxDthisGen       = idxD(thisGen);
  idxDthisGen       = idxDthisGen(idxDthisGen<=maxIdx);
  fitnesses         = poph.fitness(ismember(poph.idx, idxDthisGen));% poph.fitness(sparseIndexes(idxDthisGen));
  if numel(fitnesses)>0
    maxfitness(k)     = max(fitnesses);
    stdfitness(k)     = std(fitnesses);
    meanfitness(k)    = nanmean(fitnesses);
    medianfitness(k)  = nanmedian(fitnesses);
    coefvarfitness(k) = stdfitness(k)/meanfitness(k);
  end
  numsons           = cellfun(@numel, sons(thisGen));
  if sum(numsons)>0
    haveonesonratio(k)  = sum(numsons==1)/numel(numsons);
    havesseveralsonsratio(k)  = sum(numsons>1)/numel(numsons);
    maxsons(k)        = max(numsons);
    mediansons(k)     = median(numsons);
  end
  g=g+1;
end

gens = mingen:maxgen;

figure;
plot(gens(1:end-1), maxfitness(1:end-1), gens(1:end-1), stdfitness(1:end-1), gens(1:end-1), meanfitness(1:end-1), gens(1:end-1), medianfitness(1:end-1), gens(1:end-1), coefvarfitness(1:end-1));
legend({'max fitness', 'standard dev. of fitness', 'mean fitness', 'median fitness', 'coef. of var. (std/mean)'});
xlabel('generations'); ylabel('fitness'); grid on;
title('Fitness statistics by generation');

figure;
plot(gens(1:end-1), maxsons(1:end-1), gens(1:end-1), mediansons(1:end-1), gens(1:end-1), haveonesonratio(1:end-1), gens(1:end-1), havesseveralsonsratio(1:end-1));
legend({'max amount of sons', 'median amount of sons', 'ratio of individuals with 1 son', 'ratio of individuals with several sons'});
xlabel('generations'); ylabel('sons'); grid on;
title('Statistics of amount-of-sons-per-individual by generation');

clear fitnesses maxfitness meanfitness medianfitness coefvarfitness numsons maxsons meansons meadinsons gens;

fitnesses        = poph.fitness(toUse);
generations      = poph.generation(toUse);

figure;
hold on;
colors = cool(ngens);
kk=1;
fitnessbins= linspace(min(fitnesses), max(fitnesses), 20);
for k=mingen:maxgen
  thisGen = generations==k;
  n = hist(fitnesses(thisGen), fitnessbins);
  line(fitnessbins,n, 'Color', colors(kk,:));
  kk=kk+1;
end
grid on;
xlabel('fitness');
ylabel('count');
title('Shift of fitness histogram by generation');

figure;
hold on;
kk=1;
fitnessesbinedges = [-inf (fitnessbins(1:end-1)+diff(fitnessbins)/2) inf];
colors = cool(numel(fitnessbins));
%colors = colors(end:-1:1,:);
genbins= mingen:(maxgen-1);
for k=2:numel(fitnessesbinedges)
  thisFit = (fitnesses>fitnessesbinedges(k-1)) & (fitnesses<=fitnessesbinedges(k));
  n = hist(generations(thisFit), genbins);
  line(genbins, n, 'Color', colors(kk,:));
  kk=kk+1;
end
grid on;
xlabel('generation');
ylabel('count');
title('histogram of generations by fitness');
legend(arrayfun(@(x)sprintf('fitness %d', x), fitnessbins, 'uniformoutput', false));

