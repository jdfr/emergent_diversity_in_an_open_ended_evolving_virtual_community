function poph = loadSimulation(basedir, loadTree, makeParents, makeSons)
%Loads the results of a genetic algorithm simulation. Inputs:
%   -basedir: the base directory of the simulation (i.e., the directory
%             where poblacion.txt and arbol.txt are placed)
%   -loadTree: If true (default), load arbol.txt along with poblacion.txt
%   -makeParents: If true (loadTree implied), load parents's indexes in
%                 the struct poph
%   -makeSons: If true (loadTree and makeParents implied), load sons's
%              indexes in the struct poph
%Output: population history structure. Fields:
%  -poph.generation: array containing the generation for each individual
%  -poph.rangeid: array containing the range identifier for each individual
%  -poph.individual: array containing the index (inside its corresponding
%                    generation) for each individual
%  -poph.fitness: array containing the fitness for each individual
%  -poph.xpos: array containing the x position of each individual
%  -poph.genome: array containing the genome for each individual
%  -poph.idx: array containing an unique identifier for each individual,
%             based on the combination (poph.generation, poph.rangeid,
%             poph.individual).
%             BEWARE: rangeid SHOULD NEVER BE 0, IF NOT, THE ALGORITHMS
%                     USING poph WILL NOT WORK PROPERLY
%  -poph.triplet2idx: function to convert from triplet
%                     (generation, rangeid, index) to identifiers in
%                     poph.idx
%  -poph.tree: structure containing genealogy:
%    -poph.tree(.gD, .rD, iD, idxD): generation, rangeid, index and
%                                    identifier for each individual
%    -poph.tree(.gA, .rA, iA, idxA): generation, rangeid, index and
%                                    identifier for the ancestor of the
%                                    each individual identified by
%                                    poph.tree.idxD
%    -poph.tree.change: explanation of changes from the ancestor's genome
%                       (identified by corresponding poph.tree.idxA) to
%                       the individual's genome (identified by
%                       corresponding poph.tree.idxD)

if ~exist('loadTree', 'var')
  loadTree = true;
end

if (~exist('makeParents', 'var'))
  makeParents = loadTree;
end

if (~exist('makeSons', 'var'))
  makeSons = makeParents;
end
if iscell(basedir)
  [poph gens]       = gluePopulations(basedir);
else
  poph              = reapPopulationHistory(basedir);
end
if loadTree
  if iscell(basedir)
    poph.tree       = glueTrees(basedir, gens);
  else
    poph.tree       = reapGenealogyData(basedir);
  end
    %descendants may potentially have larger indexes nd1 and nd2 than
    %population and ancestors, so nd1 and nd2 are calculated for them first
    [poph.tree.idxD poph.nd1 poph.nd2]  = triplet2idx(poph.tree.gD,    poph.tree.rD, poph.tree.iD);
    poph.tree.idxA                      = triplet2idx(poph.tree.gA,    poph.tree.rA, poph.tree.iA,    poph.nd1, poph.nd2);
    poph.idx                            = triplet2idx(poph.generation, poph.rangeid, poph.individual, poph.nd1, poph.nd2);
    if makeParents
        %     try
        %       %get parent indexes
        %       sparseIndexes         = sparse(1,poph.tree.idxD+1,1:numel(poph.tree.idxD));
        %       poph.tree.parents     = int32(full(sparseIndexes(poph.tree.idxA+1)));
        %       clear sparseIndexes;
        %     catch ME
        %       clear sparseIndexes;
        %       if strcmpi(ME.identifier, 'MATLAB:nomem')
        %         %that sparse matrix was too much, poor old MATLAB! Let's do it the
        %         %memory-safe way
        [nevermind, poph.tree.parents] = ismember(poph.tree.idxA, poph.tree.idxD); %#ok<ASGLU>
        poph.tree.parents              = poph.tree.parents';
        clear nevermind;
        %       else
        %         rethrow(ME);
        %       end
        %     end
    end
    if makeSons
        parents = poph.tree.parents;
        sons    = cell(size(parents));
        for z=1:int32(numel(parents))
            pz = parents(z);
            if pz>0
                sons{pz}(end+1,1) = z;
            end
        end
        poph.tree.sons = sons;
    end
    %this is to alleviate the problem of file size
    poph = fillIdenticals(poph, 'genome');
    if isfield(poph, 'minGenome')
      poph = fillIdenticals(poph, 'minGenome');
    end
else
    [poph.idx poph.nd1 poph.nd2]        = triplet2idx(poph.generation, poph.rangeid, poph.individual);
end
poph.triplet2idx                      = triplet2idxFUN(poph.nd1, poph.nd2);

if iscell(basedir)
  load([basedir{1} filesep '..' filesep 'estado.mat'], 'ACIParams');
else
  load([basedir    filesep '..' filesep 'estado.mat'], 'ACIParams');
end
poph.ACIParams = ACIParams;

function FUN = triplet2idxFUN(nd1,nd2)
FUN = @(gens,rinds,inds)triplet2idx(gens,rinds,inds,nd1,nd2);

function [idxs nd1 nd2] = triplet2idx(gens,rinds,inds,nd1,nd2)
if nargin<4
    nd1          = nextDigit(inds);
end
idxs           = inds + nd1 * rinds;
if nargin<5
    nd2          = nextDigit(idxs);
end
idxs           = idxs + nd2 * gens;

function nd = nextDigit(values)
nd = power(10, (1+ceil(log10(max(values)))));

function g = reapGenealogyData(basedir)
%load in structure 'g' the contents of arbol.txt

%given a list of values, this function provides the minimum power of ten
%that is greater than the maximum value
% nextDigit = @(values) power(10, (1+ceil(log10(max(values)))));

farb = fopen([basedir filesep 'arbol.txt'], 'r');
datos = textscan(farb, '%f %f %f %s %f %f %f %f %f', 1); %#ok<ASGLU>
fclose(farb);

bufsize = 8192*2-5;

farb = fopen([basedir filesep 'arbol.txt'], 'r');
if isempty(datos{end}) || isnan(datos{end})
    datos = textscan(farb, '%f %f %f %s %f %f %f %*[^\n]', 'BufSize', bufsize);
    fclose(farb);
    g = struct('gD', [], ... %D==Descendant
        'rD', [], ...
        'iD', [], ...
        ...%           'idxDescendant', [], ...
        'idxD', [], ...
        'change', [], ...
        'gA', [], ... %A==Ancestor
        'rA', [], ...
        'iA', [], ...
        'idxA', [] ...
        );
    g.gD      = datos{1};
    g.rD      = datos{2};
    g.iD      = datos{3};
    g.change  = datos{4};
    g.gA      = datos{5};
    g.rA      = datos{6};
    g.iA      = datos{7};
else
    datos = textscan(farb, '%f %f %f %s %f %f %f %f %f %*[^\n]');
    fclose(farb);
    g = struct('gD', [], ... %D==Descendant
        'rD', [], ...
        'iD', [], ...
        ...%           'idxDescendant', [], ...
        'idxD', [], ...
        'change', [], ...
        'gA', [], ... %A==Ancestor
        'rA', [], ...
        'iA', [], ...
        'idxA', [], ...
        'distOrigen', [], ...
        'disruption', [] ...
        );
    g.gD      = datos{1};
    g.rD      = datos{2};
    g.iD      = datos{3};
    g.change  = datos{4};
    g.gA      = datos{5};
    g.rA      = datos{6};
    g.iA      = datos{7};
    g.distOrigen      = datos{8};
    g.disruption      = datos{9};

end

function poph = reapPopulationHistory(basedir)

namefile = [basedir filesep 'poblacion.txt'];

compfile = [basedir filesep 'compression_9.mat'];
if exist(compfile, 'file')
  comp = load(compfile);
  comp = comp.compression.compression;
  compressionField = {'compression9',  {comp}};
else
  compressionField = {};
end

fpob  = fopen(namefile, 'r');
datos = textscan(fpob, '%f %f %f %f %f %f %f %f %f %f %f %s %s', 1); clear nevermind; %#ok<ASGLU>
fclose(fpob);
fpob  = fopen(namefile, 'r');

bufsize = 8192*8-5;

if iscell(datos{end}) && not(isempty(datos{end})) && not(isempty(datos{end}{1}))

      datos = textscan(fpob, '%f %f %f %f %f %f %f %f %f %f %f %s %s', 'BufSize', bufsize);
      if isempty(datos) || isempty(datos{1})
          error('Data not properly read!!!!');
      end
      fclose(fpob);
      poph  = struct('generation',    (datos{1}), ...
          'rangeid',       (datos{2}), ...
          'individual',    (datos{3}), ...
          'fitness',       (datos{4}), ...
          'xpos',          (datos{5}), ...
          'height',        (datos{6}), ...
          'width',         (datos{7}), ...
          'nbranches',     (datos{8}), ...
          'nleafs',        (datos{9}), ...
          'nPixelsInSoil', (datos{10}), ...
          'nleafsOnTop',   (datos{11}), ...
          'minGenome',     {datos{12}}, ...
          'genome',        {datos{13}}, ...
          'idx',           {[]}, ...
          compressionField{:}, ...
          'numFieldsPop',  {[]}, ...
          'tree',          {[]}, ...
          'nd1',           {[]}, ...
          'nd2',           {[]}, ...
          'triplet2idx',   {[]}  ...
          );
  
else
  
  fpob  = fopen(namefile, 'r');
  datos = textscan(fpob, '%f %f %f %f %f %f %f %f %f %f %f %[^\n]', 1, 'BufSize', bufsize); clear nevermind; %#ok<ASGLU>
  fclose(fpob);
  fpob  = fopen(namefile, 'r');

  %this is to handle both old and new file formats
  if isempty(datos{end}) || (~ischar(datos{end}{1}))
      datos = textscan(fpob, '%f %f %f %f %f %[^\n]');
      if isempty(datos) || isempty(datos{1})
          error('Data not properly read!!!!');
      end
      fclose(fpob);
      poph  = struct('generation',   (datos{1}), ...
          'rangeid',      (datos{2}), ...
          'individual',   (datos{3}), ...
          'fitness',      (datos{4}), ...
          'xpos',         (datos{5}), ...
          'genome',       {datos{6}}, ...
          'idx',          {[]}, ...
          'numFieldsPop', {[]}, ...
          'tree',         {[]}, ...
          'nd1',          {[]}, ...
          'nd2',          {[]}, ...
          'triplet2idx',  {[]}  ...
          );
  else
      datos = textscan(fpob, '%f %f %f %f %f %f %f %f %f %f %f %[^\n]', 'BufSize', bufsize);
      if isempty(datos) || isempty(datos{1})
          error('Data not properly read!!!!');
      end
      fclose(fpob);
      poph  = struct('generation',    (datos{1}), ...
          'rangeid',       (datos{2}), ...
          'individual',    (datos{3}), ...
          'fitness',       (datos{4}), ...
          'xpos',          (datos{5}), ...
          'height',        (datos{6}), ...
          'width',         (datos{7}), ...
          'nbranches',     (datos{8}), ...
          'nleafs',        (datos{9}), ...
          'nPixelsInSoil', (datos{10}), ...
          'nleafsOnTop',   (datos{11}), ...
          'genome',        {datos{12}}, ...
          'idx',           {[]}, ...
          'numFieldsPop',  {[]}, ...
          'tree',          {[]}, ...
          'nd1',           {[]}, ...
          'nd2',           {[]}, ...
          'triplet2idx',   {[]}  ...
          );
  end
end


poph.numFieldsPop  = find(strcmp('numFieldsPop', fieldnames(poph)), 1, 'first')-1;

function [poph genss] = gluePopulations(basedirs)

doNotErrorForDifferences = ~((numel(basedirs)>2) && strcmp(basedirs{end}, 'doError'));
if strcmp(basedirs{end}, 'doError')
  basedirs = basedirs(1:end-1);
end

%doNotErrorForDifferences = (numel(basedirs)>2) && strcmp(basedirs{end}, 'doNotError');
% if doNotErrorForDifferences
%   basedirs = basedirs(1:end-1);
% end
genss = cell(numel(basedirs)-1,1);

poph = reapPopulationHistory(basedirs{1});
for k=2:numel(basedirs)
  poph2          = reapPopulationHistory(basedirs{k});
  [poph gens]    = glueTwoPopulations(poph, poph2, doNotErrorForDifferences);
  genss{k-1}     = gens;
end

function [poph gens] = glueTwoPopulations(poph, poph2, doNotErrorForDifferences)

maxgen1 = max(poph.generation);
mingen2 = min(poph2.generation);

if maxgen1~=mingen2
  error('maxgen1 <%d> and mingen2 <%d> do not agree!!!!', maxgen1, mingen2);
end

%check that last and first populations are the same
idx1    = find(poph.generation==maxgen1);
idx2    = find(poph2.generation==mingen2);
fnames  = fieldnames(poph);

nErr    = 0;
lastK   = -inf;
nIErr   = 0;

for k=1:numel(idx1)
  for m=1:poph.numFieldsPop-1
    name = fnames{m};
    if ~isequal(poph.(name)(idx1(k)), poph2.(name)(idx2(k)))
      str = sprintf('Discrepancy detected for individual %d of last generation %d in first simulation, for field %s. value1=%s, value2=%s\n', k, maxgen1, name, any2str(poph.(name)(idx1(k))), any2str(poph2.(name)(idx2(k))));
      if doNotErrorForDifferences
        fprintf(str);
        nErr = nErr+1;
        if lastK~=k
          nIErr = nIErr+1;
        end
        lastK = k;
      else
        error(str); %#ok<SPERR>
      end
    end
  end
end

if nErr>0
  fprintf('   TOTAL NUMBER OF DISCREPANCIES: %d\n', nErr);
  fprintf('   TOTAL NUMBER OF INDIVIDUALS WITH DISCREPANCIES: %d\n', nIErr);
end

%remove last population population in first simulation, and shift the generation
%number in the rest
for m=1:poph.numFieldsPop-1
  name               = fnames{m};
  poph.(name)(idx1)  = [];
end
newmingen2 = mingen2;

clear idx1;
%glue both populations together

% %remove first population in second simulation, and shift the generation
% %number in the rest
% for m=1:poph.numFieldsPop-1
%   name               = fnames{m};
%   poph2.(name)(idx2) = [];
% end
%
% newmingen2       = min(poph2.generation);
% poph2.generation = poph2.generation+(-newmingen2+maxgen1+1);

clear idx2;

%glue both populations
for m=1:poph.numFieldsPop-1
  name         = fnames{m};
  poph.(name)  = [poph.(name); poph2.(name)];
  poph2.(name) = [];
end

gens = struct('maxgen1', maxgen1, 'mingen2', mingen2, 'newmingen2', newmingen2);

function tree = glueTrees(basedirs, genss)
tree    = reapGenealogyData(basedirs{1});
doNotErrorForDifferences = ~((numel(basedirs)>2) && strcmp(basedirs{end}, 'doError')); %#ok<NASGU>
if strcmp(basedirs{end}, 'doError')
  basedirs = basedirs(1:end-1);
end
% doNotErrorForDifferences = (numel(basedirs)>2) && strcmp(basedirs{end}, 'doNotError');
% if doNotErrorForDifferences
%   basedirs = basedirs(1:end-1);
% end

for k=2:numel(basedirs)
  tree2 = reapGenealogyData(basedirs{k});
  tree  = glueTwoTrees(tree, tree2, genss{k-1}, basedirs{k});
end

function tree = glueTwoTrees(tree, tree2, gens, basedir)

fnames  = fieldnames(tree);

toDelete1 = find(tree.gD>gens.maxgen1);
if any(toDelete1)
  for m=1:numel(fnames)
    name                     = fnames{m};
    if ~isempty(tree.(name))
      tree.(name)(toDelete1) = []; %#ok<FNDSB>
    end
  end
end

clear toDelete1;

mingD = min(tree2.gD);

if mingD~=gens.mingen2
  error('arbol.txt and poblacion.txt do not agree for <%s>!!!!', basedir);
end

toDelete2 = find(tree2.gD==mingD);
for m=1:numel(fnames)
  name                     = fnames{m};
  if ~isempty(tree.(name))
    tree2.(name)(toDelete2) = []; %#ok<FNDSB>
  end
end

clear toDelete2;

tree2.gA = tree2.gA+(-mingD+gens.maxgen1);
tree2.gD = tree2.gD+(-mingD+gens.maxgen1);
for m=1:numel(fnames)
  name           = fnames{m};
  if ~isempty(tree.(name))
    tree.(name)  = [tree.(name); tree2.(name)];
    tree2.(name) = [];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function poph = fillIdenticals(poph, fld)

g = poph.(fld);
poph.(fld) = [];
p = poph.tree.parents;

identical = cellfun('prodofsize', g)==1;
identical(identical) = strcmp('I', g(identical));
%taking advantage of the fact that individuals are ordered by time
for k=1:numel(g)
  if identical(k)
    g(k) = g(p(k));
  end
end

poph.(fld) = g;
