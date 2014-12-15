function subnomdir = ACI(params)
startTime = clock;

%default values for several parameters
params = putDefaultParameters(params);

temp = fieldnames(params.save.disp);
params.save.disp.all = cellfun(@(x)x, struct2cell(params.save.disp));
params.save.disp.any = any(params.save.disp.all);
params.save.disp.nms = temp(params.save.disp.all);
clear temp;

params.fitnessFunctionH = str2func(params.fitnessFunction);
if any(params.variableNumSonsMode=='#')
  sharp = find(params.variableNumSonsMode=='#');
  params.variableNumSonsInterval = eval(params.variableNumSonsMode(sharp+1:end));
  params.variableNumSonsMode     = params.variableNumSonsMode(1:sharp-1);
else
  params.variableNumSonsInterval = [];
end

terclus = params.terclus;
params  = rmfield(params, 'terclus');

params.nomdir = strrep(strrep(params.nomdir, '\', filesep), '\', filesep); %#ok<NODEF> %make sure that file separators are OK
params.nomdir = [params.nomdir,filesep,num2str(params.num_iteracion),'=',num2str(params.presion),'_',num2str(params.ensayo)];
mkdir(params.nomdir);
subnomdir = params.nomdir;
params.nomdir = [params.nomdir,filesep];

terclus.elapsed = 0;

params.fpob = fopen([params.nomdir,'poblacion.txt'],'a');
params.farb = fopen([params.nomdir,'arbol.txt'],'a');
params.fres = fopen([params.nomdir,'resultado.txt'],'a');
params.ftst = fopen([params.nomdir,'timestats.txt'],'a');
% fale = fopen([nomdir,'aleatorio.txt'],'a');
params.allfs = {params.fpob; params.farb; params.fres; params.ftst};
if params.save.disp.any
  params.fdisp = fopen([params.nomdir,'disps.txt'],'a');
  params.allfs{end+1} = params.fdisp;
end

params.evaluarOptions = {params.T};

%this struct holds a destacated first field: genome, that identifies the
%individual, and other fields, that are the individual's attributes
params.indivStruct = struct('genome', {{[]}}, 'fitness', {nan}, 'randN',{nan}, 'raster', {{[]}}, 'dimensions', {[nan nan]}, 'offset', {[nan nan]}, 'maxiy', {{[]}}, 'isLeaf', {{[]}}, 'counts', {[nan nan nan]}, 'mingnm',{{[]}});
params.numFieldsToSkip = 3; %skip first fields 'genome', 'fitness' and 'randN' when reaping results
%it is convenient to pre-fabricate this structure
jobStruct = {'G' 'indexNewG' 'job' 'rangeid' 'generation'};
% legend: G==individuals of this range
%         indexNewG== indexes in G of truly new individuals (only they will be evaluated)
%         job== job handle in charge of evaluating this range
%         rangeid== the id if this range
%         generation== the generation
jobStruct = reshape([jobStruct; cell(size(jobStruct))], 1, []);
params.jobStruct = struct(jobStruct{:}); clear jobStruct;

%here, we will store the best individual of the simulation. It cannot be
%safely stored in the population, since it might be overwritten
mejorEnConserva = struct('gen', [], 'rind', [], 'num', [], 'data',  {params.indivStruct});
jobsActive = cell(0,1);
try
    maxActive  = terclus.concJobs; %#ok<NODEF>
    terclus.rangeSizes = diff(terclus.ranges,1,2)+1;

    P = repmatStruct(params.indivStruct, [params.tam_poblacion,1]);
    if isfield(params, 'initialPopulation') && (~isempty(params.initialPopulation))
      P.genome(1:end) = params.initialPopulation(1:end);
      params          = rmfield(params, 'initialPopulation');
    else
      for i=1:params.tam_poblacion
          P.genome{i} = generar_cadena(params.long_cadena,params.n_anidamientos);  % individuos aleatorios
      end
    end

    fprintf('PLEASE WAIT WHILE INITIAL POPULATION IS EVALUATED...\n');
    [P terclus params] = startPopulation(params, terclus, P, params.initialGeneration);
    if params.plotDendroTree
%       HierarchicalClustering(params.nomdir,P,params.initialGeneration)
        tosave = ['P' num2str(params.initialGeneration, '%03g')];
        eval([tosave ' = P;']);
        save([params.nomdir filesep tosave],tosave);
    end
    i=1;
    for na=1:numel(terclus.rangeSizes)
        for nb=1:terclus.rangeSizes(na)
            fprintf(params.farb,'%d %d %d = 0 0 0\n', params.initialGeneration, na, i);     % escribir genealogía de población inicial
            %fprintf(params.farb,'%d %d %d = 0 0 0 0 0\n', params.initialGeneration, na, i);     % escribir genealogía de población inicial
            printIndividual(params,params.initialGeneration, na, i, P, [], []);
%             fprintf(fale,
            i=i+1;
        end
    end
    params = saveDisps([], P, [], params, 0);
    %for each range, the generation it was last updated
    generationInRange = repmat(params.initialGeneration, (size(terclus.rangeSizes)));
    %find the best in the population
    mejorEnConserva = findMejor(P, mejorEnConserva, terclus.ranges, generationInRange);
    fprintf('TIME ELAPSED: %g\n', etime(clock, startTime));

    doTermination = false;
    lastTime = clock;
    for cont=(params.initialGeneration+1):(params.initialGeneration+params.generaciones) %for each generation
        for rind=1:size(terclus.ranges,1) %for each range

            %if max amount of concurrent jobs is reached, wait for some to end,
            %and store it (or them, if several end at once)
            if numel(jobsActive)>=maxActive
                decisionFun = @any; %wait for some job to end
            else
                decisionFun = @alwaysTrue; %reap jobs if they are finished, but do not wait for them
            end
            [ended states] = waitForJobs(params.ftst, jobsActive, decisionFun, terclus);
            %if any job ended: reaping time!
            tooSlowEvaluation = etime(clock, lastTime)>params.tooLongEvTime;
            lastTime = clock;
            if any(ended)
                endedI = find(ended);
                %reap them
                [jobsActive terclus] = reapJobs(params, jobsActive, endedI, states(ended), terclus);
                reapedids = cellfun(@(x)sprintf('G=%gR=%g, ', x.generation, x.rangeid), jobsActive(ended), 'uniformoutput', false);
                reapedids = horzcat(reapedids{:});
                %store results in population
                [jobsActive, P, generationInRange, mejorEnConserva, params] = updateGA(params, terclus, endedI, jobsActive, P, generationInRange, mejorEnConserva, parentdata, startTime);
                znb = whos('P');

                stopIfTooTall     = (~isempty(params.stopIfTooTall)) && any((P.dimensions(:,1)-P.offset(:,1))>params.stopIfTooTall);
                if stopIfTooTall || tooSlowEvaluation
                  mprintf(params.ftst, 'We stopped prematurely!!! stopIfTooTall=%s, tooSlowEvaluation=%s\n, BYTESIZE=%s,', any2str(stopIfTooTall), any2str(tooSlowEvaluation), any2str(znb.bytes));
                  clear znb;
                  doTermination=true;
                  break;
                end
                
                mprintf(params.ftst, 'REAPED  JOBS FOR %s\n   TIME ELAPSED: %f, TIME IN CLUSTER: %f, BYTESIZE=%f,\n', reapedids, etime(clock, startTime), terclus.elapsed, any2str(znb.bytes));
                clear znb;
                % MANDAR AL CLUSTER PARA PINTAR GRAFOS PIXEL
            end
            %launch a new job!!!!
            nj = params.jobStruct;
            [nj.G nj.indexNewG terclus P parentdata]  = evolucion(params, P,cont,rind,generationInRange,terclus,mejorEnConserva); % proxima generacion
            if isempty(nj.G.genome)
              mprintf(params.ftst, 'SHUTTING DOWN THE SIMULATION, SINCE NO NEW POPULATION HAS BEEN GENERATED!!!\n');
%               cellfun(@fclose, params.allfs);
%               fprintf('closed all files...\n');
              jobsActive          = carefulDestroy(jobsActive); %#ok<NASGU>
%               fprintf('destroyed all jobs...\n');
              doTermination=true;
              break;
            end
            nj.generation = cont;
            nj.rangeid    = rind;
            nj.job        = evaluarRemoteLite(terclus, nj.G.genome(nj.indexNewG), params);
            [longs, numGs] = cellfun(@(x)deal(numel(x), sum(x=='G')), nj.G.genome(nj.indexNewG));
            jobsActive{end+1} = nj; %#ok<AGROW>
            if params.variablePopSize
              popsize = numel(nj.G.genome);
            else
              popsize = numel(P.genome);
            end
            znb = whos('P');
            mprintf(params.ftst, 'SPAWNED GEN=%g, RANGE=%g, POPSIZE=%d Nº OF STRINGS: %g\n   TIME ELAPSED=%f, TIME IN CLUSTER=%f, BYTESIZE=%f,\n', cont, rind, popsize, numel(nj.indexNewG), etime(clock, startTime), terclus.elapsed, any2str(znb.bytes));
            clear znb;
            mprintf(params.ftst, '   NEW STRING LENGTH   STATS: MIN=%g, MAX=%g, MEAN=%g, MEDIAN=%g\n', min(longs), max(longs), nanmean(longs), nanmedian(longs));
            mprintf(params.ftst, '   NEW STRING G NUMBER STATS: MIN=%g, MAX=%g, MEAN=%g, MEDIAN=%g\n', min(numGs), max(numGs), nanmean(numGs), nanmedian(numGs));
        end
        if doTermination
          break;
        end
    end

    if ~doTermination
        mprintf(params.ftst, 'BEFORE FINAL REAPING:\n   TIME ELAPSED=%f, TIME IN CLUSTER=%f\n', etime(clock, startTime), terclus.elapsed);
        %wait for all extant jobs
        [nevermind states] = waitForJobs(params.ftst, jobsActive, @all, terclus); %#ok<ASGLU>
        %reap all extant jobs
        [jobsActive terclus] = reapJobs(params, jobsActive, 1:numel(jobsActive), states, terclus);
        reapedids = cellfun(@(x)sprintf('G=%gR=%g, ', x.generation, x.rangeid), jobsActive, 'uniformoutput', false);
        reapedids = horzcat(reapedids{:});
        mprintf(params.ftst, 'AT LAST, REAPED JOBS FOR %s\n   TIME ELAPSED: %g, TIME IN CLUSTER: %g\n', reapedids, etime(clock, startTime), terclus.elapsed);
        %store results in population
        [jobsActive, P, generationInRange, mejorEnConserva, params] = updateGA(params, terclus, 1:numel(jobsActive), jobsActive, P, generationInRange, mejorEnConserva, parentdata, startTime); %#ok<ASGLU>
        writeMejor(params.fres, mejorEnConserva); %write best individual for this generation
        fprintf('Best in this run: GEN=%g RANGE=%g FIT=%g IDX=%g <%s>\n', mejorEnConserva.gen, mejorEnConserva.rind, mejorEnConserva.data.fitness, mejorEnConserva.num, mejorEnConserva.data.genome{1});
    end

catch ME
    cellfun(@fclose, params.allfs);
    fprintf('closed all files...\n');
    %  destroy(findJob(terclus.jobMgr, 'Tag', terclus.tag));
    jobsActive          = carefulDestroy(jobsActive); %#ok<NASGU>
    fprintf('destroyed all jobs...\n');
    rethrow(ME);
end
cellfun(@fclose, params.allfs);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initializes population
function   [P terclus params] = startPopulation(params, terclus, P, Generation)
%get fitnesses for all

ranges          = terclus.ranges;
cj              = terclus.concJobs;
totalj          = size(ranges,1);
numconc         = 0;
jobs            = cell(1,totalj);
activeJobs      = zeros(0,1);
lP              = 0;
while (lP<totalj)
    %launch jobs
    while (lP<totalj) && (numconc<=cj)
        lP            = lP+1;
        numconc       = numconc+1;
        nj            = params.jobStruct;
        nj.indexNewG  = 1:(diff(ranges(lP,:))+1);%ranges(lP,1):ranges(lP,2);
        nj.G          = getSubStruct(P, ranges(lP,:), 'intensive');
        
        nj.rangeid    = lP;
        nj.generation = Generation;
        nj.job        = evaluarRemoteLite(terclus, nj.G.genome, params);
        jobs{lP}      = nj;
        activeJobs = [activeJobs; lP]; %#ok<AGROW>
    end

    %wait for any job to finish
    [ended states]          = waitForJobs(params.ftst, jobs(activeJobs), @any, terclus);
    %reap finished jobs
    [jobs terclus] = reapJobs(params, jobs, activeJobs(ended), states(ended), terclus);

    endedI = find(ended);
    for k=1:numel(endedI) %dispose non-needed fields
        jobs{activeJobs(endedI(k))}.indexNewG = [];
        jobs{activeJobs(endedI(k))}.job       = [];
    end
    numconc                 = numconc - sum(ended);
    activeJobs              = activeJobs(~ended);
    if size(activeJobs,1)>size(activeJobs,2) %this may happen sometimes
        activeJobs            = activeJobs'; %just to be safe
    end
end

%wait for all extant jobs to finish
[nevermind states]      = waitForJobs(params.ftst, jobs(activeJobs), @all, terclus); %#ok<ASGLU>
%reap finished jobs
[jobs terclus] = reapJobs(params, jobs, activeJobs, states, terclus);
for k=1:numel(activeJobs) %dispose non-needed fields
    jobs{activeJobs(k)}.indexNewG = [];
    jobs{activeJobs(k)}.job       = [];
end

%update population
for k=1:numel(jobs)
  P = assignSubStruct(P, ranges(k,:), 'intensive', jobs{k}.G);
end

% [xmax,xmin,ymax,ymin] = rango(P.ramas);
if isfield(params, 'initialPos') && (~isempty(params.initialPos))
  [P colores params] = params.fitnessFunctionH(P,params,Generation, [], false, params.initialPos);
  params      = rmfield(params, 'initialPos');
else
  [P colores params] = params.fitnessFunctionH(P,params,Generation, []);
end
% f = landscape2([ramas{:}],xmax,xmin,ymax,ymin,nomdir,Generation);
% % f = fitnessParalelizado(P.ramas,xmax,xmin,ymax,ymin,terclus);


% A = pixelGraph(P.ramas,xmin,xmax,ymin,ymax);
% sparseM = landscape(A);
% f = fitness(sparseM,size([ramas{:}],2));
clear jobs;

% boundingbox = [xmin,xmax,ymin,ymax];
% ranGEMA     = [xmax,xmin,ymax,ymin];
if(params.plotting)
    drawtrees(P, colores,params, Generation);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [jobstructs terclus] = reapJobs(params, jobstructs, toReap, states, terclus)
for z=1:numel(toReap)
    e = toReap(z);
    jrid = jobstructs{e}.rangeid;
    switch states{z}
        case {'failed', 'destroyed', 'unavailable', 'notjob'}
            mprintf(params.ftst, 'JOB %s FOR RANGE %g = [%g ... %g] !!!!!\n', upper(states{e}), jrid, terclus.ranges(jrid,1), terclus.ranges(jrid,2));
            %terclus = addElapsedTime(terclus, nan);
        case 'finished'
            job         = jobstructs{e}.job;
            limitsByW   = job.UserData{1};
            ntasks      = numel(job.Tasks);
            indexes     = jobstructs{e}.indexNewG;
            for k=1:ntasks
                if isempty(job.Tasks(k).ErrorMessage)
                    goon = true;
                    output = job.Tasks(k).OutputArguments;
                    %there seems to be a race condition: sometimes output can be
                    %empty, supposedly because the job has not yet transferred the
                    %output. The workaround is simple, but it's stupid to have to do
                    %it
                    if numel(output)<4
                        hayprob = true;
                        pause(terclus.pauseTime);
                        output = job.Tasks(k).OutputArguments;
                        if numel(output)<4
                            mprintf(params.ftst, 'ERROR!!! THERE SEEMS TO BE A CONDITION PREVENTING TO REAP GENERATION %g RANGE %g TASK %g!!!\nERROR: %s\nOUTPUT: %s\n', jobstructs{e}.generation, jobstructs{e}.rangeid, k, showError(job.Tasks(k).Error), any2str(job.Tasks(k).OutputArguments));
                            goon = false;
                        end
                    else
                        hayprob = false;
                    end
                    if goon
                        if hayprob
                            mprintf(params.ftst, 'WE HAVE BEEN ABLE TO RECOVER FROM THE CONDITION :-)\n');
                        end
                        %fill results
                        fnames = fieldnames(jobstructs{e}.G);
                        numFieldsToSkip = params.numFieldsToSkip;
                        for zz=(1+numFieldsToSkip):numel(fnames) 
                          %output{z-1} corresponds to the field fnames{z}, in
                          %the positions indexes(limitsByW(k,1):limitsByW(k,2))  
                          jobstructs{e}.G.(fnames{zz})(indexes(limitsByW(k,1):limitsByW(k,2)),:) = output{zz-(numFieldsToSkip-1)};
                        end
                        terclus = addElapsedTime(terclus, output{1});
                    else
                        if hayprob
                            mprintf(params.ftst, 'WE HAVEN''T BEEN ABLE TO RECOVER FROM THE CONDITION :-(\n');
                        end
                    end
                else
                    err = job.Tasks(k).Error;
                    mprintf(params.ftst, 'JOB %g HAS TASK %g FINISHED WITH ERROR!!! %s', jrid, k, showError(err));
                end
            end

    end
    if (~strcmp('notjob', states{z})) && (isa(jobstructs{e}.job, 'distcomp.job') || isa(jobstructs{e}.job, 'distcomp.simplejob'))
        destroy(jobstructs{e}.job);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function terclus = addElapsedTime(terclus, elapsed)
%if ~any(isnan(elapsed))
terclus.elapsed = terclus.elapsed+sum(elapsed);
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mejorEnConserva = findMejor(P, mejorEnConserva, ranges, generationInRange)
[maxi, mejor] = max(P.fitness);               %#ok<ASGLU> % mejor individuo
% mejor = noplanos(mejor);
if (~isempty(mejor)) %&& ( isempty(mejorEnConserva.data.fitness) || (maxi>=mejorEnConserva.data.fitness) )
    mejorEnConserva.rind  = findRange(ranges, mejor);
    mejorEnConserva.gen   = generationInRange(mejorEnConserva.rind);
    mejorEnConserva.num   = mejor;
    mejorEnConserva.data  = getSubStruct(P, mejor, 'extensive');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeMejor(fres, mejorEnConserva)
fprintf(fres,'%g %g %g %1.6f %s\n',mejorEnConserva.gen, mejorEnConserva.rind, mejorEnConserva.num, mejorEnConserva.data.fitness, mejorEnConserva.data.genome{1}); % actualizamos el fichero de resultados

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t = alwaysTrue(varargin)
t = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = findRange(ranges, idx)
solvfun = @(ix) find((ix>=ranges(:,1)) & (ix<=ranges(:,2)));
if numel(idx)==1
    r = solvfun(idx);
else
    r = arrayfun(solvfun, idx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%update variables after a job (some jobs) is (are) reaped
function [jobsActive, P, generationInRange, mejorEnConserva, params] = updateGA(params, terclus, toUpdate, jobsActive, P, generationInRange, mejorEnConserva, parentdata, startTime)
ranges = terclus.ranges;
for k=1:numel(toUpdate)
    ei    = toUpdate(k);
    rpri  = jobsActive{ei}.rangeid;
    gen   = jobsActive{ei}.generation;
    if params.variablePopSize
      newP   = jobsActive{ei}.G;
    else
      newP   = assignSubStruct(P, ranges(rpri,:), 'intensive', jobsActive{ei}.G);
    end
    params = saveDisps(P, newP, parentdata, params, gen);
    generationInRange(rpri) = gen;
%     [xmax,xmin,ymax,ymin] = rango(newP.ramas);
%     boundingbox = [xmin,xmax,ymin,ymax];
%     ranGEMA     = [xmax,xmin,ymax,ymin];
    mprintf('      NOW, THE FITNESS FUNCTION... (TIME ELAPSED: %s)\n', mat2str(etime(clock, startTime)));
    [newP colores params] = params.fitnessFunctionH(newP,params,gen,parentdata);
    mprintf('      THE FITNESS FUNCTION IS FINISHED... (TIME ELAPSED: %s)\n', mat2str(etime(clock, startTime)));
    if(params.plotting)
        drawtrees(newP, colores,params, gen);
    end
        
    %recalculate best individual
    mejorEnConserva = findMejor(newP, mejorEnConserva, ranges, generationInRange); %#ok<ASGLU>
    %output stage
    if rpri==1
        writeMejor(params.fres, mejorEnConserva); %write best individual for this generation
    end
    for z=ranges(rpri,1):ranges(rpri,2)
      printIndividual(params,gen, rpri, z, newP, P, parentdata);
    end
    P = newP;
    clear newP;
    if params.plotDendroTree
%       HierarchicalClustering(params.nomdir,newP,gen);
        tosave = ['P' num2str(gen, '%03g')];
        eval([tosave ' = P;']);
        save([params.nomdir filesep tosave],tosave);
    end
    fprintf(params.ftst, '   POPULATION GENERATIONS UPDATE: GENERATION=%g, RANGE=%g, ALL=%s    (TIME ELAPSED: %s)\n', gen, rpri, mat2str(generationInRange), mat2str(etime(clock, startTime)));
end
%wipe out reaped jobs
jobsActive(toUpdate) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = saveDisps(P, newP, parentdata, params, gen)

if params.save.disp.any
  fprintf(params.fdisp, 'GENERATION %d\n', gen);
  if isempty(P)
    %first generation
    [newdists indicesClasificados] = generateDistanceMatricesFromP(newP, params.save.disp.nms);
    if numel(indicesClasificados) ~= numel(newP.genome)
      save([params.nomdir 'errorraro.mat'], 'newdists', 'indicesClasificados', 'P', 'newP', 'parentdata', 'params', 'gen');
      error('In the initial generation, we are not prepared (but we can be changed) to deal with plants with underground branches!!!!\n%s %s\n', mat2str(numel(indicesClasificados)), mat2str(numel(newP.genome)));
    end
  else
    %first, generate new distance matrices
    olddists = params.save.disp.dists;
    newdists = cell(size(olddists));
    for z=1:numel(newdists)
      newdists{z} = permuteSQM(parentdata.idxInP, olddists{z});
    end
    %next, re-generate entries in the matrix for changed individuals
    [newdists indicesClasificados] = updateDistanceMatricesFromP(newP, newdists, find(parentdata.mutated), find(parentdata.maskNew), params.save.disp.nms);
  end
  %save full distance matrix (the fact that entries for plants with
  %underground branches are not updated is not important: they will not
  %breed, thus they will not affect the next generation)
  params.save.disp.dists = newdists; 
  for k=1:numel(newdists)
    %statistics are compiled without taking into account plants with
    %underground branches
    if not(isempty(newdists{k}))
      newdists{k} = newdists{k}(indicesClasificados, indicesClasificados);
    end
    if isempty(newdists{k})
      v = 0;
    else
      v = newdists{k}(:);
    end
    st = struct(...
      'max',    max(v), ...
      'var',    var(v), ...
      'mean',   mean(v), ...
      'median', median(v) ...
      );
    fprintf(params.fdisp, [params.save.disp.nms{k} ' ' any2str(st) '\n']);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%destroy distributed jobs carefully
function jobsActive = carefulDestroy(jobsActive)
if (~isempty(jobsActive)) && iscell(jobsActive)
    marked = false(size(jobsActive));
    for k=1:numel(jobsActive)
        if isstruct(jobsActive{k}) && isfield(jobsActive{k}, 'job') && ( isa(jobsActive{k}.job, 'distcomp.job') || isa(jobsActive{k}.job, 'distcomp.simplejob')) && (~isempty(jobsActive{k}.job))
            for m=1:numel(jobsActive{k}.job)
                %if ~strcmp(jobsActive{k}.job(m).State, 'destroyed')
                destroy(jobsActive{k}.job(m));
                %end
            end
            marked(k) = true;
        end
    end
    jobsActive(marked) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = showError(error)
if isempty(error)
    str='<EMPTY ERROR>';
elseif (~isstruct(error)) || (~all(ismember(lower(fieldnames(error)), {'message', 'identifier', 'stack'})))
    str=['<NO STRUCT ERROR: ' any2str(error) '>'];
else
    lines = cell(size(error.stack));
    for m=1:numel(error.stack)
        lines{m} = sprintf('    File: %s\n  Name: %s\n  Line: %g\n', error.stack(m).file, error.stack(m).name, error.stack(m).line);
    end
    str = horzcat(...
        sprintf('ERROR:\n  Message:    <<%s>>\n  Identifier: <<%s>>\n  Stack:\n', error.message, error.identifier), ...
        lines{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [timeSpent] = drawtrees(P, colores, params, Generation)%#ok<INUSD> 
if (~params.plottingSample)
  timeSpent = 0;
  return;
end

drawing.figsize = [1000 1000]; %#ok<UNRCH>
drawing.FontSize = 4;
drawing.rowcols = [6 6];

switch params.plotSampleMode
  case 'best'
    [nevermind idx] = sort(P.fitness, 'descend'); %#ok<ASGLU>
    clear nevermind;
    idx = idx(1:min(numel(P.raster), prod(drawing.rowcols)));
  case 'firstIndexes'
    idx =     1:min(numel(P.raster), prod(drawing.rowcols));
  otherwise
    error('plotSampleMode not understood <%s>!!!', any2str(params.plotSampleMode));
end

%draw a panel of trees
timeStart = clock;
fig = figure('Visible', 'off', 'Units', 'pixels', 'PaperUnits', 'points', 'PaperPosition', [0 0 drawing.figsize]);
axis equal; %#ok<UNRCH>
for z=1:numel(idx)
  i = idx(z);
  if (numel(P.raster{i})>0) && (numel(P.raster{i}{1})>0) %~any(isnan(P.ramas{i}))
    subplot(drawing.rowcols(1),drawing.rowcols(2),z,'align');
    switch params.fitnessFunction
      case 'fitnessSimilarity'
        params = drawTwoTreesTogether(params, P, i);
      otherwise
        drawtree(P.raster{i},  P.dimensions(i,:), P.offset(i,:), 0, i, colores, 'fit');
    end
    title(sprintf('%s\n%g=>%g',P.genome{i},i,P.fitness(i)),'FontSize',drawing.FontSize);%4);
  elseif (numel(P.raster{i})>0) && (numel(P.raster{i}{1})==0)
    title(sprintf('THIS GENOME PRODUCES NO TREE: %s', P.genome{i}));
  else
    title(sprintf('ERROR IN EVALUATION FOR %s', P.genome{i}));
  end
end
drawnow;
saveas(fig,[params.nomdir,'sample' num2str(Generation, '%03g')],'png');
close(fig);


timeSpent = etime(clock, timeStart);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newjob = evaluarRemoteLite(terclus, G, params)%, drawing)

nargs = 8;
if isempty(terclus.jobMgr) || isempty(G)
  outputs = cell(1,nargs);
  limitsByW = [1 numel(G)];
  if ~params.plotting 
      returnIndividuals = 0;
  else
      returnIndividuals = numel(G); %return all individuals
  end
  [outputs{:}] = evaluar(G,returnIndividuals,params.evaluarOptions{:});
  newjob = struct('Tasks', {struct('OutputArguments', {outputs}, 'ErrorMessage', {''})}, ...
                  'State', {'finished'}, ...
                  'UserData', {{limitsByW, returnIndividuals}}, ...
                  'fakeJob', {true});
else

  concW       = terclus.concurrentWorkers; %number of concurrent workers

  [limitsByW indsByW] = calculateRanges(numel(G), concW);
  if ~params.plotting 
      returnIndividuals = zeros(size(indsByW));
  else
      returnIndividuals = indsByW; %return all individuals
  end

  tempseed = rand(params.randmethod);
  newjob  = createJob(terclus.jobMgr,terclus.jobArgs{:});
  rand(params.randmethod, tempseed);
  newjob.UserData = {limitsByW, returnIndividuals};
  for k=1:numel(indsByW)
    createTask(newjob,@evaluar, nargs, {G(limitsByW(k,1):limitsByW(k,2)), returnIndividuals(k), params.evaluarOptions{:}});
  end
  fprintf('He creado %d tasks\n', numel(indsByW));
  submit(newjob);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printIndividual(params, generation, rangeid, z, P, Pant, parentdata) %z is the individual
%fprintf(params.fpob,'%d %d %d %.10f %d %s\n',generation,rangeid,z,P.fitness(z),P.randN(z), P.genome{z});
height = max(P.maxiy{z}); if isempty(height); height = int32(0); end
if params.correc.zeroHeightOK
  cond = P.maxiy{z}>=0;
else
  cond = P.maxiy{z}>0;
end
if isempty(Pant)
  gnm = P.genome{z};
  mng = P.mingnm{z};
else
  if parentdata.mutated(z)
    gnm = P.genome{z};
  else
    gnm = 'I';
  end
  if parentdata.maskNew(z)
    mng = P.mingnm{z};
  else
    mng = 'I';
  end
end
fprintf(params.fpob,'%d %d %d %.10f %d %d %d %d %d %d %d %s %s\n',generation,rangeid,z,P.fitness(z),P.randN(z), height, P.dimensions(z,2), P.counts(z,1), P.counts(z,2), P.counts(z,3), sum(P.isLeaf{z}(cond)), mng, gnm); %P.mingnm{z}, P.genome{z});
%fprintf(params.fpob,'%d %d %d %.10f %10d %10d %10d %10d %10d %10d %10d %s %s\n',generation,rangeid,z,P.fitness(z),P.randN(z), height, P.dimensions(z,2), P.counts(z,1), P.counts(z,2), P.counts(z,3), sum(P.isLeaf{z}(P.maxiy{z}>0)), mng, gnm); %P.mingnm{z}, P.genome{z});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = drawTwoTreesTogether(params, P, numTree)

rasterTree = P.raster{numTree};
offsetTree = P.offset(numTree,:);

if isfield(params, 'rastermaster')
  rastermaster = params.rastermaster;
else
  rastermaster = zeros(params.treeToMatch{5}, 'uint8');
  rastermaster(sub2ind(size(rastermaster), params.treeToMatch{2}{1}, params.treeToMatch{2}{2})) = 2;
  params.rastermaster = rastermaster;
end

offsetMaster = params.treeToMatch{3};

ys = rasterTree{1}-offsetTree(1)+offsetMaster(1);
xs = rasterTree{2}-offsetTree(2)+offsetMaster(2);

minY = min(ys);
minX = min(xs);

if minY<=0
  rastermaster = [zeros(1-minY, size(rastermaster,2), 'uint8'); rastermaster];
  ys = ys+1-minY;
end
if minX<=0
  rastermaster = [zeros(size(rastermaster,1), 1-minX, 'uint8'), rastermaster];
  xs = xs+1-minX;
end

maxY = max(ys);
maxX = max(xs);

if maxY>size(rastermaster,1)
  rastermaster = [rastermaster; zeros(maxY-size(rastermaster,1), size(rastermaster,2), 'uint8')];
end
if maxX>size(rastermaster,2)
  rastermaster = [rastermaster, zeros(size(rastermaster,1), maxX-size(rastermaster,2), 'uint8')];
end

rastermaster(sub2ind(size(rastermaster), ys, xs)) = 3;

imshow(rastermaster, [1 1 1; 0 0 0; 1 0 0; 0 0 1], 'InitialMagnification', 'fit');
set(gca, 'YDir', 'normal');
