function fits = loadFitness(basedirs)

fits = cell(size(basedirs));

for k=1:numel(basedirs)
  fprintf('loading %03d/%03d %s...\n', k, numel(basedirs), basedirs{k});
  f = fopen([basedirs{k} filesep 'poblacion.txt']);
  fits{k} = textscan(f, '%f %*f %*f %f %*[^\n]', 'BufSize', 8192*2-5);
  fclose(f);
end
