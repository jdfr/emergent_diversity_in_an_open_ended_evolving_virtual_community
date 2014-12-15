function numberOfSpecies(basedir,GenI,GenF,parallel)
if not(exist('parallel', 'var'))
  parallel = [];
end

if isempty(GenF)
  generations = GenI(:)';
else
  generations = (GenI:GenF);
end

if isempty(parallel)
  for i = generation%GenI : GenF
    doAGeneration(basedir, i);
  end
elseif islogical(parallel)
  error('''Parallel'' argument must be the number of jobs!!!!!');
else
  jobArgs = {'Tag', 'numberOfSpecies', 'PathDependencies',  {pwd}};
  gs = generations';%(GenI:GenF)';
  [ranges rangesizes] = calculateRanges(numel(gs), parallel);
  gens = mat2cell(gs, rangesizes, 1);
  bass = repmat({basedir}, size(gens));
  fprintf('GOING PARALLEL...\n');
  dummy = my_dfeval(@doSeveralGenerations, bass, gens, jobArgs{:});
  fprintf('RECEIVING PARALLEL...\n');
end

function dummy = doSeveralGenerations(basedir, gens)
for i=gens(1):gens(end)
  doAGeneration(basedir, i);
end
dummy = gens;


function doAGeneration(basedir, i)
    tosave1 = ['individualsSpecies' num2str(i, '%03g')];
    if exist([basedir filesep tosave1], 'file')
      return;
    end
    i
    load(sprintf('%s/Z%03d.mat',basedir,i));
    Z = eval(sprintf('Z%03d', i));
    load(sprintf('%s/matrixDist%03d.mat',basedir,i));
    matrixDist = eval(sprintf('matrixDist%03d', i));
    load(sprintf('%s/indicesClasificados%03d.mat',basedir,i));
    indicesClasificados = eval(sprintf('indicesClasificados%03d', i));

    [individualsSpecies proto offsetproto sim numGroups] = numberGroup(basedir,i,Z,matrixDist,indicesClasificados);

    
    tosave2 = ['proto' num2str(i, '%03g')];
    tosave3 = ['offsetproto' num2str(i, '%03g')];
    tosave4 = ['sim' num2str(i, '%03g')];
    tosave5 = ['numGroups' num2str(i, '%03g')];

    eval([tosave1 ' = individualsSpecies;']);
    eval([tosave2 ' = proto;']);
    eval([tosave3 ' = offsetproto;']);
    eval([tosave4 ' = sim;']);
    eval([tosave5 ' = numGroups;']);

    save([basedir filesep tosave1],tosave1);
    save('-v7.3',[basedir filesep tosave2],tosave2)
    save([basedir filesep tosave3],tosave3);
    save([basedir filesep tosave4],tosave4);
    save([basedir filesep tosave5],tosave5);

    clear('individualsSpecies','proto','offsetproto','sim','numGroups', tosave1,tosave2,tosave3,tosave4,tosave5);
