function disruptStatsScript(datafile, whattoshow, dst, names, types, rndseeds, gens, times, numrnd)

%disruptStatsScript('lsystemdani\dataset\T2dataset_weno_jaccard_151all.mat', [], 'MP', {'lsystemdani\dataset\G1F1\uno\1=0.001_3', 'lsystemdani\dataset\G0.75F1\uno\1=0.001_1'}, {2 2}, [], {151 151}, [], []);   

%dst: 'MP' for JACCARD, 'LOGBN' for abs(log10(nbranches1)-log10(nbranches2))   

%disruptStatsScript('T2_G1F1_3.mat', [], {'rene/dataset/G1F1/uno/1=0.001_3'}, 2);
%disruptStatsScript({'T2dataset.mat' []}, [], 'defaultT2');
%disruptStatsScript('lsystemdani\dataset\T2dataset_weno_LOGBN.mat', [], 'LOGBN');

%main('local', {'', 'replay'}, 'genomas\porhacer\G0.8F1\2009_Jul_29_12_21_07\2', 'resultados\2009_Jul_29_12_21_07_harsh_T3',     'ensayos', '1', 'mode_op_insercion', '''one_level''', 'prob', '[0 0.05 0 0.05 0 0.05 0.05 0]', 'generaciones', '500', 'tam_poblacion', '1', 'initialPopulation', '{''[G]''}', 'initialGeneration', '0');
%main('local', {'', 'replay'}, 'genomas\porhacer\G1F1\2009_Jul_29_14_22_12\5',   'resultados\2009_Jul_29_14_22_12_veryharsh_T3', 'ensayos', '1', 'mode_op_insercion', '''one_level''', 'prob', '[0 0.05 0 0.05 0 0.05 0.05 0]', 'generaciones', '500', 'tam_poblacion', '1', 'initialPopulation', '{''[G]''}', 'initialGeneration', '0');
%main('local', {'', 'replay'}, 'genomas\porhacer\G0.5F1\2009_Jul_29_12_16_39\3', 'resultados\2009_Jul_29_12_16_39_mild_T3',      'ensayos', '1', 'mode_op_insercion', '''one_level''', 'prob', '[0 0.05 0 0.05 0 0.05 0.05 0]', 'generaciones', '500', 'tam_poblacion', '1', 'initialPopulation', '{''[G]''}', 'initialGeneration', '0');

%main('embedded', {'', 'replay'}, 'rene/dataset/G0.8F1/uno', 'dataset/G0.8F1/T3uno', 'ensayos', '1', 'mode_op_insercion', '''one_level''', 'prob', '[0 0.05 0 0.05 0 0.05 0.05 0]', 'generaciones', '500', 'tam_poblacion', '1', 'initialPopulation', '{''[G]''}', 'initialGeneration', '0', 'stopIfTooTall', '20000', 'tooLongEvTime', '15000', 'plotDendroTree', 'true', 'seedRandom', {sum(clock)*1e5});
%main('embedded', {'', 'replay'}, 'rene/dataset/G1F1/uno',   'dataset/G1F1/T3uno',   'ensayos', '1', 'mode_op_insercion', '''one_level''', 'prob', '[0 0.05 0 0.05 0 0.05 0.05 0]', 'generaciones', '500', 'tam_poblacion', '1', 'initialPopulation', '{''[G]''}', 'initialGeneration', '0', 'stopIfTooTall', '20000', 'tooLongEvTime', '15000', 'plotDendroTree', 'true', 'seedRandom', {sum(clock)*1e5});
%main('embedded', {'', 'replay'}, 'rene/dataset/G0.5F1/dos', 'dataset/G0.5F1/T3dos', 'ensayos', '1', 'mode_op_insercion', '''one_level''', 'prob', '[0 0.05 0 0.05 0 0.05 0.05 0]', 'generaciones', '500', 'tam_poblacion', '1', 'initialPopulation', '{''[G]''}', 'initialGeneration', '0', 'stopIfTooTall', '20000', 'tooLongEvTime', '15000', 'plotDendroTree', 'true', 'seedRandom', {sum(clock)*1e5});

%cd rene/src ; matlab -nodisplay
%addpath(pwd); cd ..
%main('embedded', {'', 'replay'}, 'rene/exampleG0.9F1/dummy', 'exampleG0.9F1/G0.9F1/newsim1', 'ensayos', '1', 'tam_poblacion', '1', 'initialGeneration', '0', 'seedRandom', {sum(clock)*1e5});
%main('embedded', {'', 'replay'}, 'rene/exampleG0.9F1/dummy', 'exampleG0.9F1/G0.9F1/newsim2', 'ensayos', '1', 'tam_poblacion', '1', 'initialGeneration', '0', 'seedRandom', {sum(clock)*1e5});
%main('embedded', {'', 'replay'}, 'rene/exampleG0.9F1/dummy', 'exampleG0.9F1/G0.9F1/newsim3', 'ensayos', '1', 'tam_poblacion', '1', 'initialGeneration', '0', 'seedRandom', {sum(clock)*1e5});
%main('embedded', {'', 'replay'}, 'rene/exampleG0.9F1/dummy', 'exampleG0.9F1/G0.9F1/newsim4', 'ensayos', '1', 'tam_poblacion', '1', 'initialGeneration', '0', 'seedRandom', {sum(clock)*1e5});
%main('embedded', {'', 'replay'}, 'rene/exampleG0.9F1/dummy', 'exampleG0.9F1/G0.9F1/newsim5', 'ensayos', '1', 'tam_poblacion', '1', 'initialGeneration', '0', 'seedRandom', {sum(clock)*1e5});




% ACIParams.prob                                               = [0 0.05 0 0.05 0 0.05 0.05 0];
% ACIParams.prob                                               = [1 0 0 0 0.05 0.05 0.05 0];
% ACIParams.prob                                               = [0 0.05 0.05 0 0 0.05 0.05 0];

% disruptStatsScript({'T2new.mat', []}, [], 'defaultT2');
% disruptStatsScript({'T3new.mat', []}, [], 'defaultT3');

if ispc
  veryharsh = 'genomas\porhacer\G1F1\2009_Jul_29_14_22_12\5\1=0.001_1';
  harsh     = 'genomas\porhacer\G0.8F1\2009_Jul_29_12_21_07\2\1=0.001_1';
  mild      = 'genomas\porhacer\G0.5F1\2009_Jul_29_12_16_39\3\1=0.001_1';
  veryharshT3 = 'genomas\porhacer\G1F1\2009_Jul_29_14_22_12_veryharsh_T3\1=0.0001_1';
  harshT3     = 'genomas\porhacer\G0.8F1\2009_Jul_29_12_21_07_harsh_T3\1=0.001_1';
  mildT3      = 'genomas\porhacer\G0.5F1\2009_Jul_29_12_16_39_mild_T3\1=0.001_2';
else
%   veryharsh = './2009_Jul_29_14_22_12/5/1=0.001_1';
%   harsh     = './2009_Jul_29_12_21_07/2/1=0.001_1';
%   mild      = './2009_Jul_29_12_16_39/3/1=0.001_1';
%   veryharshT3 = './2009_Jul_29_14_22_12_veryharsh_T3/1=0.0001_1';
%   harshT3     = './2009_Jul_29_12_21_07_harsh_T3/1=0.001_1';
%   mildT3      = './2009_Jul_29_12_16_39_mild_T3/1=0.001_2';
  veryharsh = 'rene/dataset/G1F1/uno/1=0.001_1';
  harsh     = 'rene/dataset/G0.75F1/uno/1=0.001_1';
  mild      = 'rene/dataset/G0.5F1/tres/1=0.001_1';
  veryharshT3 = 'rene/dataset/G1F1/T3unobis/1=0.001_1';
  harshT3     = 'rene/dataset/G0.75F1/T3uno/1=0.001_1';
  mildT3      = 'rene/dataset/G0.5F1/T3dosbis/1=0.001_1';
end
if iscell(datafile)
  numjobs = datafile{2};
  datafile = datafile{1};
else
  numjobs = 54;
end

if exist(datafile, 'file')
  S=load(datafile);
  data = S.data;
  clear S;
  if nargin>3
    if ischar(names)
      if strcmp(names, 'defaultT2')
        names = {veryharsh, harsh, mild}';
      elseif strcmp(names, 'defaultT3')
        names = {veryharshT3, harshT3, mildT3}';
      end
    end
    for k=1:numel(names)
      data(k).name = names{k};
    end
  end
else
  narg = 4;
  if nargin<narg
    error('AT LEAST PROVIDE TE FIRST OPTION');
    %names = {veryharsh, harsh, mild}';
    %names = {veryharshT3, harshT3, mildT3}';
  else
    if ischar(names)
      if strcmp(names, 'defaultT2')
        names = {veryharsh, harsh, mild}';
      elseif strcmp(names, 'defaultT3')
        names = {veryharshT3, harshT3, mildT3}';
      end
    end
    if not(iscell(names))
      names = {names};
    end
    names = names(:);
  end
  narg = narg+1;
  if nargin<narg || isempty(types)
    types = array2cell([2 2 2]');
  else
    if not(iscell(types))
      types = {types};
    end
    types = types(:);
  end
  narg = narg+1;
  if nargin<narg || isempty(rndseeds)
    rndseeds = array2cell(rand(size(names))*1e9);
  else
    if not(iscell(rndseeds))
      rndseeds = array2cell(rndseeds);
    end
    if numel(rndseeds)==1
      rndseeds = repmat(rndseeds, size(names));
    end
  end
  narg = narg+1;
  if nargin<narg || isempty(gens)
    gens = cell(size(names));
  else
    if not(iscell(gens))
      gens = array2cell(gens);
    end
    gens = gens(:);
    if numel(gens)==1
      gens = repmat(gens, size(names));
    end
  end
  narg = narg+1;
  if nargin<narg||isempty(times)
    times = repmat({100}, size(names));
  else
    if not(iscell(times))
      times = array2cell(times);
    end
    if numel(times)==1
      times = repmat(times, size(names));
    end
  end
  narg = narg+1;
  if nargin<narg||isempty(numrnd)
    numrnd = cell(size(names));
  else
    if not(iscell(numrnd))
      numrnd = array2cell(numrnd);
    end
    if numel(numrnd)==1
      numrnd = repmat(numrnd, size(names));
    end
  end
  
  dsts = repmat({dst}, size(names));
  
  data = struct('name', names, 'types', types, 'numgen', gens, 'popall', [], ...
                'pop', [], 'popstats', [], 'numrnd', numrnd, 'randompop', [], ...
                'randompopstats', [], 'ntimes', times, 'rndseed', rndseeds, ...
                'seedspop', [], 'seedsrnd', [], ...
                'disrup', [], 'disrupMean', [], 'disrupR', [], 'disrupMeanR', [], ...
                'plotLabels', [], 'done', array2cell(false(size(names))), ...
                'dst', dsts ...
                );
  save(datafile, 'data');
end
  
if ispc
  doparallel = [];
else
  doparallel = struct('numJobs', numjobs, ...%54, ...
                      'jobArgs', {{'Tag', 'tree_disruption', 'PathDependencies',  {pwd}}});%, 'MaximumNumberOfWorkers', 40};);
end

for k=1:numel(data)
  if data(k).done
    continue;
  end
  fprintf('For %s:\n', data(k).name);
  if isempty(data(k).seedsrnd)
    fprintf('PREPARING DATA\n');%: loading %s\n', data(k).name);
    clear poph; poph = loadSimulation(data(k).name, false, false, false);
    g = data(k).numgen;
    if isempty(g)
      g = max(poph.generation);
    end
    toGet  = find(g==poph.generation);
    toGet  = toGet(poph.nPixelsInSoil(toGet)==0);
    fnames = fieldnames(poph);
    values = struct2cell(poph);
    fnames = fnames(1:poph.numFieldsPop); 
    values = values(1:poph.numFieldsPop);
    for z=1:numel(values)
      values{z} = values{z}(toGet);
    end
    nv = [fnames(:)'; array2cell(values(:)')];
    data(k).popall = struct(nv{:});
    data(k).pop = data(k).popall.genome;
    [data(k).pop poppidxs] = unique(data(k).pop);
    for za=1:numel(fnames)
      data(k).popall.(fnames{za}) = data(k).popall.(fnames{za})(poppidxs);
    end
    data(k).popall = rmfield(data(k).popall, 'genome');

    fprintf('  get statistics...\n');
    data(k).popstats = getGenomeStatistics(data(k).pop);

    nr = data(k).numrnd;
    if isempty(nr)
      nr = numel(data(k).pop);
    end
    rs = data(k).rndseed;
    if not(isempty(rs))
      rand('twister', rs);
      randn('state', rs);
    end
    fprintf('  generate random genomes...\n');
    switch data(k).types
      case 2
        data(k).randompop = generateGenomeRTipo2(nr,data(k).popstats);
      case 3
        if any(data(k).popstats.numAnid>1)
          error('non-regular genomes detected!!!!');
        end
        data(k).randompop = generateGenomeRTipo3(nr,data(k).popstats);
    end

    fprintf('  generate random statistics...\n');
    data(k).randompopstats = getGenomeStatistics(data(k).randompop);

    fprintf('  generate random seeds...\n');
    data(k).seedspop = rand(size(data(k).pop))*1e9;
    data(k).seedsrnd = rand(size(data(k).randompop))*1e9;
    save(datafile, 'data');
  end
  
%   save('TTTTT.mat', 'data');
%   return
% 
  if isempty(data(k).disrup)
    fprintf('  calculate disruption for evolved genomes...\n');
    if not(isfield(data(k), 'dst'))
      error('data struct must have ''dst'' field. The default value (for old structures) should be ''MP'', the original distance used');
    end
    [data(k).disrupMean  data(k).plotLabels data(k).disrup]  = calculateDisruption(data(k).pop, data(k).ntimes, data(k).seedspop, data(k).dst, doparallel);
    save(datafile, 'data');
  end
  if isempty(data(k).disrupR)
    fprintf('  calculate disruption for random genomes...\n');
    if not(isfield(data(k), 'dst'))
      error('data struct must have ''dst'' field. The default value (for old structures) should be ''MP'', the original distance used');
    end
    [data(k).disrupMeanR data(k).plotLabels data(k).disrupR] = calculateDisruption(data(k).randompop, data(k).ntimes, data(k).seedspop, data(k).dst, doparallel);
    save(datafile, 'data');
  end
  fprintf('DONE!!!!!!\n');

  data(k).done = true;

  save(datafile, 'data');
end

if iscell(whattoshow)&&isempty(whattoshow)
  whattoshow = struct(...
    'idxs', [1 2], ...
    'points', {cell(1,2)}, ...
    'marker', {{'o' 'o'}}, ...
    'mkcolf', {{'none' 'none'}}, ...
    'mkcole', {{'r' 'b'}}, ...
    'mksiz', [10 10], ...
    'legend', {{'very harsh T2', 'harsh T2'}} ...
    );
end

if isempty(whattoshow)
  return
end
showrobustez(data, whattoshow);
