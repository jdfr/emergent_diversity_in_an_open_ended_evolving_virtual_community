function dispersionsTXT2MAT(basedir)

doForSeveralSimulations(@doit, basedir, [], {}, true, []);

function d = doit(basedir, d)
tic;


if exist([basedir filesep 'dispersionsss.mat'], 'file')
  fprintf('This simulation already has the processed file!!!\n');
  return;
end
if not(exist([basedir filesep 'disps.txt'], 'file'))
  fprintf('This simulation has not generated the appropriated RAW file!!!!');
  return;
end

f = fopen([basedir filesep 'disps.txt']);

text = textscan(f, '%s', 1, 'BufSize', 1024*1024*10, 'delimiter', char(29), 'endOfLine', char(29));

fclose(f);

fprintf('disps.txt LOADED (%s)...\n', mat2str(toc));

%tokens=regexp(t{1}{1}, 'GENERATION ([0-9]+)[ \n\r]*([a-zA-Z0-9]+)[ \n\r]*(struct\([a-zA-Z0-9 \{\}\[\]\''\.\,\;]+\))', 'tokens');
tokens=regexp(text{1}{1}, '(GENERATION ([0-9]+)|[a-zA-Z0-9]+ struct\([a-zA-Z0-9 \{\}\[\]\''\.\,\;]+\))', 'tokens');

gens=regexp(text{1}{1}, 'GENERATION ([0-9]+)', 'tokens');

fprintf('text PARSED (%s)...\n', mat2str(toc));

initgen = str2double(gens{1}{1});
lastgen = str2double(gens{end}{1});

gens = (initgen:lastgen)';

init = nan(size(gens));
names = {};
dispersion = struct;
maxg = struct;

strg = 'GENERATION ';

idx=0;
g=nan;
for k=1:numel(tokens)
  t = tokens{k}{1};
  if isempty(strmatch(strg, t))
    p = find(t==' ', 1, 'first');
    fld = t(1:p-1);
    val = eval(t(p+1:end));
    maxg.(fld) = g;
    if not(isfield(dispersion, fld))
      names{end+1} = fld; %#ok<AGROW>
      dispersion.(fld) = val;
      fnames = fieldnames(val);
      for z=1:numel(fnames)
        dispersion.(fld).(fnames{z}) = init;
      end
    end
    fnames = fieldnames(val);
    for z=1:numel(fnames)
      dispersion.(fld).(fnames{z})(idx) = val.(fnames{z});
    end
  else
    idx=idx+1;
    lastg = g;
    g = str2double(t(numel(strg)+1:end));
    if not(isnan(lastg)) && ((g-lastg)~=1)
      error('generations are not in order!!!!');
    end
  end
end

maxg = struct2cell(maxg);
maxg = max([maxg{:}]);

mask = gens<=maxg;
gens = gens(mask);
f1 = fieldnames(dispersion);
for k=1:numel(f1)
  f2 = fieldnames(dispersion.(f1{k}));
  for z=1:numel(f2)
    v = dispersion.(f1{k}).(f2{z});
    v = v(mask);
    dispersion.(f1{k}).(f2{z}) = v;
  end
end

data = struct('gens', {gens}, 'names', {names}, 'dispersion', {dispersion}); %#ok<NASGU>

fprintf('text PROCESSED (%s)...\n', mat2str(toc));

save([basedir filesep 'dispersionsss.mat'], 'data');

fprintf('data SAVED (%s)...\n', mat2str(toc));
