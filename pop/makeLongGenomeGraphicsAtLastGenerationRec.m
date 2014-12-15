function makeLongGenomeGraphicsAtLastGenerationRec(basedir, fileToTest)
fprintf('Probing dir <%s>...\n', basedir);
dirr        = dir(basedir);
subdirs     = arrayfun(@(x) x.isdir && (x.name(1)~='.'), dirr);
if     exist([basedir filesep 'arbol.txt'], 'file') && ...
       exist([basedir filesep 'poblacion.txt'], 'file') && ...
       exist([basedir filesep '..' filesep 'estado.mat'], 'file') && ...
       (~exist([basedir filesep fileToTest], 'file'))
  fprintf('Generating graphs for dir <%s>...\n', basedir);
  try
    makeLongGenomeGraphics(basedir, inf, true);
    fprintf('DONE! Generated graphs for dir <%s>...\n', basedir);
  catch ME
    fprintf('HOW SAD!!!! An Error has been generated:\n');
    fprintf(showError(ME));
  end
end
if any(subdirs)
  subdirs = dirr(subdirs);
  for k=1:numel(subdirs)
    makeLongGenomeGraphicsAtLastGenerationRec([basedir filesep subdirs(k).name], fileToTest);
  end
end

