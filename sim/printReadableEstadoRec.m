function printReadableEstadoRec(basedir)
fprintf('processing dir %s\n', basedir);
dirr        = dir(basedir);
subdirs     = arrayfun(@(x) x.isdir && (x.name(1)~='.'), dirr);
anyEstado   = any(arrayfun(@(x)strcmp(x.name, 'estado.mat'), dirr));
notReadable = ~any(arrayfun(@(x)strcmp(x.name, 'readableEstado.txt'), dirr));
if anyEstado && notReadable
  load([basedir filesep 'estado.mat']);
  printReadableEstado([basedir filesep 'readableEstado.txt'], ACIParams, 'ACIParams');
  clear ACIParams;
  fprintf('PRINTED STATE  %s\n', basedir);
end
if any(subdirs)
  subdirs = dirr(subdirs);
  for k=1:numel(subdirs)
    printReadableEstadoRec([basedir filesep subdirs(k).name]);
  end
end
    