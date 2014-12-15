function [varargout]= compressGemaFiles(basedir, geni, genf, prefix, option)
% compressGemaFiles('newdists\mutdist\2009_Jul_29_12_21_07\1=0.001_1', 1, 557, '', 12); fprintf('hecho1\n');
% compressGemaFiles('newdists\miniedit\2009_Jul_29_12_21_07\1=0.001_1', 1, 557, '', 12); fprintf('hecho2\n');
% compressGemaFiles('newdists\miniedit\2009_Jul_29_12_16_39\1=0.001_1', 1, 95, '', 12); fprintf('hecho3\n');
% compressGemaFiles('newdists\mutdist\2009_Jul_29_12_16_39\1=0.001_1', 1, 95, '', 12); fprintf('hecho4\n');
% compressGemaFiles('newdists\editdist\2009_Jul_29_12_16_39\1=0.001_1', 1, 95, '', 12); fprintf('hecho5\n');
% compressGemaFiles('newdists\miniedit\2009_Jul_29_14_22_12\1=0.001_1', 1, 500, '', 12); fprintf('hecho6\n');
% compressGemaFiles('newdists\mutdist\2009_Jul_29_14_22_12\1=0.001_1', 1, 500, '', 12); fprintf('hecho7\n');
% compressGemaFiles('newdists\editdist\2009_Jul_29_14_22_12\1=0.001_1', 1, 500, '', 12); fprintf('hecho8\n');


%prf={'NM_';'E1_';'E2_';'MP_'};
%b='rene/dataset/G1F1/uno/1=0.001_1'; topg=500; prf={'NM_';'E1_';'E2_';'MP_'}; for k=1:numel(prf); fprintf('GOING FOR %s\n', prf{k}); compressGemaFiles(b, 1, topg,prf{k}, 12); end; 

switch option
  case 1
    files = {'Z', 'matrixDist', 'indicesClasificados'}; 
  case 2
    files = {'sim', 'individualsSpecies', 'numGroups'};
  case 12
    files = { 'Z', 'matrixDist', 'indicesClasificados', 'sim', 'individualsSpecies', 'numGroups'}; 
  case 12.1
    files = { 'matrixDist', 'indicesClasificados', 'sim', 'individualsSpecies', 'numGroups'}; 
  case 3
    files = {'comLZfigureLongitud'};
  case 4
    files = {'matrixDist'};
  case 5
    files = {'matrixDist', 'indicesClasificados'};
  otherwise
    error('jarlll');
end

if exist([basedir filesep prefix  files{end} '.mat'],'file')
  fprintf('ALREADY DONE!!!!');
  return
end


faltan= gemaCalcFinished(basedir, geni, genf, files, prefix, true);

if not(isempty(faltan))
  fprintf('some files are missing: %s', faltan{1});
  return
end

generations = (geni:genf)';
allfiles = cell(numel(files)*numel(generations),1);
vx = [files(:)'; repmat({{cell(numel(generations), 1)}}, 1, numel(files))];
res = struct(vx{:});
res.generations = generations;
m=1;
for z=1:numel(files)
  fl = files{z};
  isIndSpecs = strcmp(fl, 'individualsSpecies');
  fprintf('LOADING DATA FOR %s\n', fl);
  for k=geni:genf
    nm = [fl num2str(k, '%03g')];
    fil = [basedir filesep prefix nm '.mat'];
    allfiles{m} = fil;
    m = m+1;
    if not(isIndSpecs)
      zs = load(fil);
      res.(fl){k} = zs.(nm);
      clear zs;
    end
  end
  if strcmp(files{z}, 'numGroups')
    fprintf('CALULATING ABRDIGED INDIVIDUALSSPECIES\n');
    individualsSpeciesJust = cell(size(res.numGroups));
    for k=geni:genf
      nm = ['individualsSpecies' num2str(k, '%03g')];
      fil = [basedir filesep prefix nm '.mat'];
      zs = load(fil);
      individualsSpeciesJust{k}= zs.(nm){res.numGroups{k}};
      clear zs;
    end
    save([basedir filesep prefix 'individualsSpeciesJust.mat'], 'individualsSpeciesJust');
    clear individualsSpeciesJust;
  end
  eval(sprintf('%s = res.%s;', files{z}, files{z}));
  if strcmp(files{z}, 'matrixDist')
    zs = whos('matrixDist');
    if zs.bytes>(0.95*(2^31))
      fprintf('We have to shrink the size of the variable. going single...\n');
      for qw=1:numel(matrixDist)
        matrixDist{qw} = single(matrixDist{qw}); %#ok<AGROW>
      end
      zs = whos('matrixDist');
      if zs.bytes>(2^31)
        fprintf('This is fucking disgusting, even using single precision, the var is over 2GB!!!!');
      else
        fprintf('Sacrifying precision, the variable fits in the only decent file format that matlab has available\n');
      end
    end
  end
  save([basedir filesep prefix  files{z} '.mat'], files{z});
  clear(files{z});
  res = rmfield(res, fl);
end

fprintf('DELETING FILES!!!!!!!!\n');
for z=1:numel(allfiles)
  delete(allfiles{z});
end
