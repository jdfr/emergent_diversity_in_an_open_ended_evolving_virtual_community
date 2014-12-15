function recalculateClusters(basedir, prefixes, newprefixes, parallelNS, pathNS, optionsNS)

genI = 1;

lastl = getLastLine([basedir filesep 'poblacion.txt']);
if isempty(lastl)
  error('Mu raro, no hay líneas con caracteres en el fichero!!!!!!');
end
genF = sscanf(lastl, '%f', 1);
if isempty(genF)
  error('Mu raro, mira linea: <%s>!!!!!!', lastl);
end
clear lastl;

filesToDecompress = {'Z', 'matrixDist', 'indicesClasificados'};
filesToGenerate = {'individualsSpecies', 'sim', 'numGroups'};

gens = genI:genF;

for k=1:numel(prefixes)
  paths = cellfun(@(x) [basedir filesep newprefixes{k} x '.mat'], filesToGenerate, 'uniformoutput', false);
  if not(all(cellfun(@(x) exist(x, 'file'), paths)))
    %copy compressed gema files, giving them the new prefix
    for z=1:numel(filesToDecompress)
      fOLD = [basedir filesep    prefixes{k} filesToDecompress{z} '.mat'];
      fNEW = [basedir filesep newprefixes{k} filesToDecompress{z} '.mat'];
      if not(exist(fNEW, 'file'))
        copyfile(fOLD, fNEW);
      end
    end
    %decompress gema files
    decompressGemaFiles(basedir, newprefixes{k}, gens, gens, filesToDecompress);
    %recalculate clusters
    numberOfSpeciesNOPROTO(basedir,gens,[],parallelNS, newprefixes{k}, pathNS, optionsNS);
    %remove decompressed files
    removeFiles(basedir, newprefixes{k}, filesToDecompress, gens);
    %compress new clustering files
    compressGemaFiles(basedir, genI, genF, newprefixes{k}, 2);
  end
end

function decompressGemaFiles(basedir, prefix, indexes, gens, files)

for k=1:numel(files)
  if not(exist([basedir filesep prefix files{k} num2str(gens(end), '%03g') '.mat'], 'file'))
    a=load([basedir filesep prefix files{k} '.mat']);
    a=a.(files{k}); %#ok<NASGU>
    for z=1:numel(indexes)
      tosave=[files{k} num2str(gens(z), '%03g')];
      eval([tosave ' = a{z};']);
      save([basedir filesep prefix tosave '.mat'],tosave);
      clear(tosave);
    end
    clear('a');
  end
end

function removeFiles(basedir, prefix, files, gens)

for k=1:numel(files)
  for z=1:numel(gens)
    delete([basedir filesep prefix files{k} num2str(gens(z), '%03g') '.mat']);
  end
  delete([basedir filesep prefix files{k} '.mat']);
end


