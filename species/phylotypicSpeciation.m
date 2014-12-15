function phylotypicSpeciation(basedir, poph, mode, umbral, prefixMatrixDist, prefix)
if not(exist('prefix', 'var'))
  prefix = '';
end
if not(exist('prefixMatrixDist', 'var'))
  prefixMatrixDist = '';
end
if not(exist('path', 'var'))
  path = pwd;
end

names = {'phyloNumSpecies', 'phyloIndsSpecs',     'phyloProtos', 'phyloSpeciesInd'};
vars  = {'numSpecies',      'individualsSpecies', 'protos',      'speciesInd'};

already = alreadyExist(basedir, prefix, umbral, names);

if already
  fprintf('ALREADY DONE !!!!\n');
  return
end

nm = [basedir filesep prefixMatrixDist 'matrixDist.mat'];
if not(exist(nm, 'file'))
  fprintf('THE FILE %s DOES NOT EXIST!!!\n', nm);
  return;
end
matrixDist = load(nm);
matrixDist = matrixDist.matrixDist;
nm = [basedir filesep prefixMatrixDist 'indicesClasificados.mat'];
if not(exist(nm, 'file'))
  fprintf('THE FILE %s DOES NOT EXIST!!!\n', nm);
  return;
end
indicesClasificados = load(nm);
indicesClasificados = indicesClasificados.indicesClasificados;
maxgen = numel(matrixDist);
ret= speciesPredecesorNOPROTO(basedir,poph,umbral,maxgen, mode, matrixDist, indicesClasificados); %#ok<NASGU>

pr = [basedir filesep mat2str(umbral) prefix];
for k=1:numel(names)
  fprintf('SAVING %s...\n', vars{k});
  eval(sprintf('%s = ret.%s;', vars{k}, vars{k}));
  save([pr names{k} '.mat'], vars{k});
  clear(vars{k});
end

function already = alreadyExist(basedir, prefix, umbral, names)

already = true;
pr = [basedir filesep mat2str(umbral) prefix];
for k=1:numel(names)
  already = exist([pr names{k} '.mat'], 'file');
  if not(already)
    return
  end
end