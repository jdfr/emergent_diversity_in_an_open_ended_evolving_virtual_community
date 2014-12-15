function faltan= gemaCalcFinished(basedir, geni, genf, mode, prefix, justonce)
if ischar(mode)
  switch mode
    case 'clustering'
      files = {'Z', 'matrixDist', 'indicesClasificados'};
    case 'numspecies'
      files = {'individualsSpecies', 'sim', 'numGroups'};
  end
elseif iscell(mode)
  files = mode;
end
  
if not(exist('prefix', 'var'))
  prefix = '';
end
if not(exist('justonce', 'var'))
  justonce = false;
end
faltan={};
if not(exist(basedir, 'dir'))
  error('directory does not exist: %s', basedir);
end
d = dir(basedir);
d = {d.name}';
for k=geni:genf
  for z=1:numel(files)
    fil = [prefix files{z} num2str(k, '%03g') '.mat'];
    if not(any(strcmp(fil, d)))
      faltan = [faltan; {fil}]; %#ok<AGROW>
      if justonce
        return
      end
    end
  end
end
