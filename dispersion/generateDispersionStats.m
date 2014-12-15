function generateDispersionStats(mode, basedir, path, recursive, torquedata)
%for several simulations, generate statistics of dispersion

if ~exist('recursive', 'var')
  recursive = false;
end
if ~exist('path', 'var')
  path = '';
end
if ~exist('torquedata', 'var')
  torquedata = [];
end
reorderfun = @(x)numel(x):-1:1;
doForSeveralSimulations(@generateDispersionStatsAux, basedir, torquedata, {mode, path}, recursive, reorderfun);


function torquedata = generateDispersionStatsAux(basedir, torquedata, mode, path)

filename = 'dispersionsss.mat';

if exist([basedir filesep filename], 'file')
  fprintf('Already calculate: this file exits: %s\n', [basedir filesep filename]);
  return
end

if isempty(torquedata) && strcmp(mode, 'torque')
  torquedata = struct(...
      'jobArgs', {{'Tag', 'manySims', 'PathDependencies',  genpath(path)}}, ...
      'jobMgr',  {findResource('scheduler', 'type', 'torque')} ...
      );
end

switch mode
  case 'embedded'
    fprintf('   LOADING SIMULATION...\n');
    poph = loadSimulation(basedir);
    fprintf('   LOADING SIMULATION DONE!!!...\n');
    morphoOnTheFly('embedded', poph, (min(poph.generation):max(poph.generation))', struct('basedir', basedir, 'filename', filename), @compileStats);
  case 'torque'
    job = my_dfevalasyncJM(...
        torquedata.jobMgr, ...
        @morphoOnTheFly, ...
        1, ...
        {'embedded'}, {basedir}, {[]}, {struct('basedir', basedir, 'filename', filename)}, {@compileStats}, ...
        torquedata.jobArgs{:} ...
        ); %#ok<NASGU>
    %error('mode not implemented yet: %s', mode);
  otherwise
    error('mode not understood: %s', mode);
end




function data = compileStats(stage, data, poph, gens, mode, thisGen, P, gen, k, Pant, nomutated, recomputed, idxInP)

switch stage(1)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'i' %'init
    x = zeros(numel(gens), 1);
    x = struct('max', x, 'var', x, 'mean', x, 'median', x); 
    names = {'pht' 'gnm' 'mng'};
    x = [names; repmat({x}, size(names))];
    data = struct(...
      'gens', {gens}, ...
      'names', {names}, ...
      'dispersion', {struct(x{:})}, ...
      'dists', {[]}, ...
      'old', {data} ...
      );
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'g' %'gen'
    if isempty(Pant)
      %first generation
      [newdists indicesClasificados] = generateDistanceMatricesFromP(P, data.names);
      if numel(indicesClasificados) ~= numel(P.genome)
        error('In the initial generation, we are not prepared (but we can be changed) to deal with plants with underground branches!!!!');
      end
    else
      %first, generate new distance matrices
      olddists = data.dists;
      newdists = cell(size(olddists));
      for z=1:numel(newdists)
        newdists{z} = permuteSQM(idxInP, olddists{z});
      end
      %next, re-generate entries in the matrix for changed individuals
      [newdists indicesClasificados] = updateDistanceMatricesFromP(P, newdists, find(not(nomutated)), recomputed, data.names);
    end
    %save full distance matrix (the fact that entries for plants with
    %underground branches are not updated is not important: they will not
    %breed, thus they will not affect the next generation)
    data.dists = newdists; 
    for z=1:numel(data.names)
      %statistics are compiled without taking into account plants with
      %underground branches
      if not(isempty(newdists{z}))
        newdists{z} = newdists{z}(indicesClasificados, indicesClasificados);
      end
      data = updateStats(data, data.names{z}, k, newdists{z});
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'e' %'end'
    old = data.old;
    data = rmfield(data, {'old', 'dists'});
    save([old.basedir filesep old.filename], 'data');
end



function data = updateStats(data, fld, k, dists)
if isempty(dists)
  v = 0;
else
  v = dists(:);
end
data.dispersion.(fld).mean(k)   = mean(v);
data.dispersion.(fld).var(k)    = var(v);
data.dispersion.(fld).median(k) = median(v);
data.dispersion.(fld).max(k)    = max(v);

