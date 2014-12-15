function discrepan = redoFitness(poph, generation,drawToScreen)
%this function is used to recalc the fitness for a generation in a
%simulation. It is useful to see what is going on

params = putDefaultParameters(poph.ACIParams);

toUse = poph.generation==generation;
P = struct('genome',      {poph.genome(toUse,:)}, ...
           'fitness',     {poph.fitness(toUse,:)}, ...
           'randN',       {poph.xpos(toUse,:)}, ...
           'raster',      {[]}, ...
           'dimensions',  {[]}, ...
           'offset',      {[]}, ...
           'maxiy',       {[]}, ...
           'isLeaf',      {[]}, ...
           'counts',      {[]});
heights       = poph.height(toUse,:);
widths        = poph.width(toUse,:);
nbranches     = poph.nbranches(toUse,:);
nleafs        = poph.nleafs(toUse,:);
nPixelsInSoil = poph.nPixelsInSoil(toUse,:);
nleafsOnTop   = poph.nleafsOnTop(toUse,:);
clear toUse;
P.raster      = cell(size(P.genome));
P.maxiy       = cell(size(P.genome));
P.isLeaf      = cell(size(P.genome));
P.dimensions  = zeros(numel(P.genome),2);
P.offset      = zeros(numel(P.genome),2);
P.counts      = zeros(numel(P.genome),3);
nchar = 0;
for k=1:numel(P.genome)
  [P.raster{k}, P.offset(k,:), P.counts(k,:), P.dimensions(k,:), P.isLeaf{k}, P.maxiy{k}]=ls2(P.genome{k}, params.T);
  if heights(k)      ~= safeMax(P.maxiy{k});            error('recorded and calculated <heihgt>        disagree for individual %d!!!\nValues %d and %d', k, heights(k),       safeMax(P.maxiy{k}));            end
  if widths(k)       ~= P.dimensions(k,2);              error('recorded and calculated <width>         disagree for individual %d!!!\nValues %d and %d', k, widths(k),        P.offset(k,2));                  end
  if nbranches(k)    ~= P.counts(k,1);                  error('recorded and calculated <nbranches>     disagree for individual %d!!!\nValues %d and %d', k, nbranches(k),     P.counts(k,1));                  end
  if nleafs(k)       ~= P.counts(k,2);                  error('recorded and calculated <nleafs>        disagree for individual %d!!!\nValues %d and %d', k, nleafs(k),        P.counts(k,2));                  end
  if nPixelsInSoil(k) ~= P.counts(k,3);                 error('recorded and calculated <nPixelsInSoil> disagree for individual %d!!!\nValues %d and %d', k, nPixelsInSoil(k), P.counts(k,3));                  end
  if nleafsOnTop(k)  ~= sum(P.isLeaf{k}(P.maxiy{k}>0)); error('recorded and calculated <nleafsOnTop>   disagree for individual %d!!!\nValues %d and %d', k, nleafsOnTop(k),   sum(P.isLeaf{k}(P.maxiy{k}>0))); end
  if nchar>0
    fprintf(repmat('\b', 1, nchar));
  end
  str = sprintf('Generated raster for tree %d', k);
  nchar = numel(str);
  fprintf(str);
end
randN         = P.randN;
fitnessBefore = P.fitness;
[P colors]    = fitness2(P,params,generation,[],drawToScreen,randN); %#ok<NASGU>
fitnessAfter  = P.fitness;

discrepan = find(fitnessBefore~=fitnessAfter);

function m = safeMax(array)
if isempty(array);  m=int32(0);  else  m=max(array); end
