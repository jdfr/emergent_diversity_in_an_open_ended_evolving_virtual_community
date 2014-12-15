function dummy = scriptgeneratedataset(prefix, newdir, option, useRandom, numensayos, subdirs)
%script to launch a simulation with some parameters, to a directory
%reflecting its alpha value

%cd rene/src2 ; matlab -nodisplay
% cd rene/src
% matlab -nodisplay
% addpath(pwd); cd ..; scriptgeneratedataset('uno', 1);
% addpath(pwd); cd ..; scriptgeneratedataset('uno', 2);
% addpath(pwd); cd ..; scriptgeneratedataset('uno', 3);
% addpath(pwd); cd ..; scriptgeneratedataset('uno', 4);
% addpath(pwd); cd ..; scriptgeneratedataset('uno', 5);
% addpath(pwd); cd ..; scriptgeneratedataset('uno', 6);
% addpath('lsystemdani\src'); scriptgeneratedataset('uno', 7);
% addpath('lsystemdani\src'); scriptgeneratedataset('uno', 8);
% addpath('lsystemdani\src'); scriptgeneratedataset('uno', 9);

%scriptgeneratedataset('uno', {1});
% sched = findResource('scheduler', 'type', 'torque'); j = findJob(sched)
% j = findJob(sched)

if isempty(prefix)
  prefix = '';
else
  prefix = [prefix filesep];
end

dummy = [];

if isempty(newdir)
  newdir = getNameByTime;
elseif iscell(newdir)
  newdir = [getNameByTime newdir{1}];
end

justSpecies = iscell(option);
if justSpecies
  option = option{1};
end
if option<0
  newbranches = true;
  option = -option;
end

if not(justSpecies)
  baseparams = {'ensayos', mat2str(numensayos), ...
                'prob', '[0 0.05 0.05 0.05 0.05 0.05 0.05 0]', ...
                'generaciones', '500', ...
                'tam_poblacion', '1', ...
                'stopIfTooTall', '20000', ...
                'tooLongEvTime', '15000', ...
                'initialPos', '5000', ...
                'PRESIONES', '0.001', ...
                'terclus.ranges', '[1 1]', ...
                'initialGeneration', '0', ...
                'initialPopulation', '{''G''}', ...
                'plotDendroTree', 'true'};

  if useRandom
    baseparams = [baseparams, {'seedRandom', {sum(clock)*1e5}}];
  end
end
  fprintf('hola1\n');
  if ispc
    path = '';
    origdir = 'genomas\porhacer\G0.5F1\2009_Jul_29_12_16_39\3';
  else
    path = 'rene/src';
    origdir = 'rene/dataset/G1F1/uno';%'rene/src/2009_Jul_29_12_16_39/3';
  end
  switch option
    case 1
      subdir = [prefix 'dataset' filesep 'G1F1'];
      addparams = {'alphaG', '1', ...
                   'alphaF', '1'};
    case 2
      subdir = [prefix 'dataset' filesep 'G0.7F1'];
      addparams = {'alphaG', '0.7', ...
                   'alphaF', '1'};
    case 3
      subdir = [prefix 'dataset' filesep 'G0.5F1'];
      addparams = {'alphaG', '0.5', ...
                   'alphaF', '1'};
    case 4
      subdir = [prefix 'dataset' filesep 'G1F0.5'];
      addparams = {'alphaG', '1', ...
                   'alphaF', '0.5'};
    case 5
      subdir = [prefix 'dataset' filesep 'G0.7F0.5'];
      addparams = {'alphaG', '0.7', ...
                   'alphaF', '0.5'};
    case 6
      subdir = [prefix 'dataset' filesep 'G0.5F0.5'];
      addparams = {'alphaG', '0.5', ...
                   'alphaF', '0.5'};
    case 7
      subdir = [prefix 'dataset' filesep 'G1F0.1'];
      addparams = {'alphaG', '1', ...
                   'alphaF', '0.1'};
    case 8
      subdir = [prefix 'dataset' filesep 'G0.7F0.1'];
      addparams = {'alphaG', '0.7', ...
                   'alphaF', '0.1'};
    case 9
      subdir = [prefix 'dataset' filesep 'G0.5F0.1'];
      addparams = {'alphaG', '0.5', ...
                   'alphaF', '0.1'};
    case 10
      subdir = [prefix 'dataset' filesep 'G0.8F1'];
      addparams = {'alphaG', '0.8', ...
                   'alphaF', '1'};
    case 11
      subdir = [prefix 'dataset' filesep 'G0.8F0.5'];
      addparams = {'alphaG', '0.8', ...
                   'alphaF', '0.5'};
    case 12
      subdir = [prefix 'dataset' filesep 'G0.8F0.1'];
      addparams = {'alphaG', '0.8', ...
                   'alphaF', '0.1'};
    case 13
      subdir = [prefix 'dataset' filesep 'G0.75F1'];
      addparams = {'alphaG', '0.75', ...
                   'alphaF', '1'};
    case 14
      subdir = [prefix 'dataset' filesep 'G0.75F0.5'];
      addparams = {'alphaG', '0.75', ...
                   'alphaF', '0.5'};
    case 15
      subdir = [prefix 'dataset' filesep 'G0.75F0.1'];
      addparams = {'alphaG', '0.75', ...
                   'alphaF', '0.1'};
    case 16
      subdir = [prefix 'dataset' filesep 'PROBOLDG0.75F1'];
      addparams = {'alphaG', '0.75', ...
                   'alphaF', '1', ...
                   'prob', '[0 0.05 0.05 0 0 0.05 0.05 0]'};
    case 17
      subdir = [prefix 'dataset' filesep 'PROBOLDG0.5F1'];
      addparams = {'alphaG', '0.5', ...
                   'alphaF', '1', ...
                   'prob', '[0 0.05 0.05 0 0 0.05 0.05 0]'};
    case 18
      subdir = [prefix 'dataset' filesep 'PROBOLDG1F1'];
      addparams = {'alphaG', '1', ...
                   'alphaF', '1', ...
                   'prob', '[0 0.05 0.05 0 0 0.05 0.05 0]'};
    case 19
      subdir = [prefix 'dataset' filesep 'PROBOLDG0.8F1'];
      addparams = {'alphaG', '0.8', ...
                   'alphaF', '1', ...
                   'prob', '[0 0.05 0.05 0 0 0.05 0.05 0]'};
    otherwise
      error('jarlll!');
  end
  
if newbranches
  subdir = [subdir 'NB'];
  addparams = [addparams, {'T',              'setfield(ACIParams.T, ''sameLevelRule'', false)', ...
                           'plotDendroTree', 'false'} ...
              ];
end

if not(justSpecies)
  fprintf('hola2\n');
  main('embedded', {subdir, 'replay'}, {origdir}, newdir, baseparams{:}, addparams{:});
end
if justSpecies
  fprintf('hola3\n');
  newd = [pwd filesep subdir filesep newdir];
  fprintf('about to make speciations for %s\n', newd);
  d = dir(newd);
  d = d([d.isdir]);
  d = {d.name}';
  fprintf('subdirs: %s\n', any2str(d));
  d = d(cellfun(@(x)x(1)~='.', d));
  if justSpecies && exist('subdirs', 'var')
    d = d(subdirs);
  end
  for k=1:numel(d)
    basedir = [newd filesep d{k}];
    makeAllSpeciations(basedir, 1, [], path, true);
  end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ahora = getNameByTime

ahora = datestr(now, 'yyyy_mmm_dd_HH_MM_SS');
