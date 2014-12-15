function [total, tiempos, numFracasados] = sumarTiemposDeComputoNoWorkers(basedir)

d = dir(basedir);

total = 0;
tiempos = [];
numFracasados = 0;

for k=1:numel(d)
  if d(k).isdir
    if d(k).name(1)~='.'
      [n1 n2 nf] = sumarTiemposDeComputoNoWorkers([basedir filesep d(k).name]);
      total = total + n1;
      numFracasados = numFracasados + nf;
      tiempos = [tiempos; n2]; %#ok<AGROW>
    end
  elseif strcmp(d(k).name, 'timestats.txt');
    fprintf('getting %s...\n', basedir);
    [nuevo isFracasado] = getTiempo([basedir filesep d(k).name]);
    numFracasados = numFracasados+isFracasado;
    total = total + nuevo;
    tiempos(end+1,1) = nuevo; %#ok<AGROW>
  end
end

function [tiempo isFracasado] = getTiempo(file)

numPopsizes = 0;
numPopsizes1 = 0;

f = fopen(file, 'r');

tiempo = -111;

while ~feof(f)
  l = fgetl(f);
  trimmed = l(~isspace(l));
  if (numel(trimmed)>12) && all(trimmed(1:12)=='TIMEELAPSED=')
    zz=find(trimmed==',', 1);
    if isempty(zz)
      error('this should never happen!!!!');
    end
    tiempo = str2double(trimmed(13:zz));
    if isnan(tiempo)
      error('error 1 processing <%s>!!!!', l);
    end
  end
  if (numel(trimmed)>8)
    pos = strfind(trimmed, 'POPSIZE=');
    if ~isempty(pos)
      pos = pos+8;
      tr = trimmed(pos:end);
      pos2 = find(upper(tr)=='N', 1, 'first');
      nn = tr(1:(pos2-1));
      popsize = str2double(nn);
      if isnan(popsize)
        error('error 2 processing <%s>, <%s>!!!!', l, nn);
      end
      numPopsizes  = numPopsizes  + 1;
      numPopsizes1 = numPopsizes1 + (popsize<2);
    end
  end
end

if tiempo==-111
  error('This should never happen!!!!!');
end

ratio = numPopsizes1/numPopsizes;

isFracasado = (ratio) > 0.95;