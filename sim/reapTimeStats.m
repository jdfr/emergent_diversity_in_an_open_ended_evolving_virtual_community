function info = reapTimeStats(runstr,gens,ranges,show)
%reap the time statistics for a given simulation. Runstr is the directory
%of the simulation

if ~exist('show', 'var')
  show = true;
end

stfields = {'minl', 'ming', 'maxl', 'maxg', 'meanl', 'meang', 'medianl', 'mediang', 'numStrings', ...
            'timeSpawned', 'timeReaped', 'totalClusterSpawned', 'totalClusterReaped', ...
            'bestFitness', ...
            };
statdim  = cell(gens, ranges);
args = reshape([stfields; repmat({statdim}, size(stfields))], 1, []);
info = struct(args{:});
clear statdim args;

nomdir = [runstr filesep];

% fpob = fopen([nomdir,'poblacion.txt'],'r');
fres = fopen([nomdir,'resultado.txt'],'r');
ftst = fopen([nomdir,'timestats.txt'],'r');

ntok = '([0-9\.eE]+)';

AUTOMATED_SPAWNED = {...
  {'numStrings',          'STRINGS: ', true}, ...
  {'timeSpawned',         'ELAPSED=',  false}, ...
  {'totalClusterSpawned', 'CLUSTER=',  true}, ...
  {'minl',                'MIN=',      false}, ...
  {'maxl',                'MAX=',      false}, ...
  {'meanl',               'MEAN=',     false}, ...
  {'medianl',             'MEDIAN=',   true}, ...
  {'ming',                'MIN=',      false}, ...
  {'maxg',                'MAX=',      false}, ...
  {'meang',               'MEAN=',     false}, ...
  {'mediang',             'MEDIAN=',   false}, ...
  };

while ~feof(ftst)
  line = fgetl(ftst);
  if ~isempty(strmatch('SPAWNED', line))
    tok   = regexp(line, ['GEN=' ntok], 'tokens', 'once');
    gen   = str2double(tok{1});
    tok   = regexp(line, ['RANGE=' ntok], 'tokens', 'once');
    range = str2double(tok{1});
    if (gen>gens) || (range>ranges)
      for k=1:numel(AUTOMATED_SPAWNED)
        if AUTOMATED_SPAWNED{k}{3}
          fgetl(ftst);
        end
      end
      continue
    end
    for k=1:numel(AUTOMATED_SPAWNED)
      FIELD   = AUTOMATED_SPAWNED{k}{1};
      LABEL   = AUTOMATED_SPAWNED{k}{2};
      READNEW = AUTOMATED_SPAWNED{k}{3};
      tok   = regexp(line, [LABEL ntok], 'tokens', 'once');
      if not(isempty(tok))
        val   = str2double(tok{1});
        info(gen,range).(FIELD) = val;
      end
      if READNEW
        line = fgetl(ftst);
      end
    end
  elseif (~isempty(strmatch('REAPED', line))) || (~isempty(strmatch('AT LAST, REAPED', line)))
    tok   = regexp(line, ['G=' ntok], 'tokens');
    ges   = cellfun(@(x)str2double(x{1}), tok);
    tok   = regexp(line, ['R=' ntok], 'tokens');
    ras   = cellfun(@(x)str2double(x{1}), tok);
    line  = fgetl(ftst);
    tok   = regexp(line, ['ELAPSED: ' ntok], 'tokens', 'once');
    elaps = str2double(tok{1});
    tok   = regexp(line, ['CLUSTER: ' ntok], 'tokens', 'once');
    clust = str2double(tok{1});
    for z=1:numel(ges)
      if (ges(z)>gens) || (ras(z)>ranges)
        continue;
      end
      info(ges(z),ras(z)).timeReaped         = elaps;
      info(ges(z),ras(z)).totalClusterReaped = clust;
    end
  elseif ~isempty(strmatch('   POPULATION', line))
    %do nothing
  elseif (~isempty(strmatch('BEFORE FINAL REAPING', line))) || ...
         (~isempty(strmatch('BEFORE FINAL DRAWING', line))) || ...
         (~isempty(strmatch('AFTER FINAL DRAWING', line)))
    fgetl(ftst);
    %do nothing
  else
    fprintf('I CANNOT UNDERSTAND THIS <%s>\n', line);
  end
end

bestFitness = zeros(1,gens);
for k=1:gens
  if feof(fres)
    break;
  end
  line = fgetl(fres);
  %info(k,1).bestFitness = sscanf(line, '%f');
  info(k,1).bestFitness = sscanf(line, '%*d %*d %*d %f');
end

fclose(ftst);
fclose(fres);

if show
  r1ts    = diff([0, [info(:,1).timeSpawned]]);
  r1tr    = diff([0, [info(:,1).timeReaped]]);
  tg      = [info(:,2).timeReaped] - [info(:,1).timeSpawned];
  r1tsc   = diff([0, [info(:,1).totalClusterSpawned]]);
  r1trc   = diff([0, [info(:,1).totalClusterReaped]]);
  numstrs = [info(:,1).numStrings];
  meanl   = [info(:,1).meanl];
  maxl    = [info(:,1).maxl];
  totlen  = meanl.*numstrs;

  tims = 1:gens;
  plot(tims, r1tsc, 'm', tims, tg, 'c', tims, meanl, 'r', tims, maxl, 'b', tims, totlen, 'k', tims, bestFitness, 'g'); grid on;
  legend({'t. cluster', 'tot. gen. time for gen.', 'mean ind. len.', 'max ind. len.', 'tot. ind. len.', 'best fitness'});
end