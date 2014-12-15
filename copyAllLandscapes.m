function idx = copyAllLandscapes(basedir, putdir, prefix, copysimname, copyidx, modes, thelast, idx)
%recursively searchs simulations, copying the last landscape of each one to
%another directory

%b = 'rene/dataset/test'; bd = cellfunc(@(x) [b x], {'110NB_old', 'LONG', 'LONGG', 'LONGOLD', 'LONGOLD2', 'OKPRFX'}); for k=1:numel(bd); copyAllLandscapes(bd{k}, [bd{k} filesep 'results'], '', true, false, {'scale'}); end

%b = 'rene/dataset/test'; bd = cellfunc(@(x) [b x], {'110NB_old', '110NB', 'LONG', 'LONGG', 'LONGOLD', 'LONGOLD2', 'OKPRFX'}); for k=1:numel(bd); copyAllLandscapes(bd{k}, [bd{k} filesep 'results'], '', true, false, cellfunc(@(x)['#generations_meanMatrixDistWithSTD_' x '.png'], {'gnm', 'mng', 'pht'})); end   

if not(exist('idx', 'var'))
  idx = 0;
end
if not(exist('thelast', 'var'))
  thelast = nan;
end

d = dir(basedir);

nm = {d.name}';
isd = [d.isdir]';
ds = cellfun(@(x)all(x(1:min(numel(x), 2))~='.'), nm);
isds = isd&ds;
dirs = nm(isds);

for k=numel(dirs):-1:1
  idx = copyAllLandscapes([basedir filesep dirs{k}], putdir, prefix, copysimname, copyidx, modes, thelast, idx);
end

fprintf('Checking %s...\n', basedir);

if all(ismember({'poblacion.txt', 'arbol.txt', 'timestats.txt'}, nm))
  fprintf('Seems to be a simulation in %s...\n', basedir);
  nms = cellfun('prodofsize', nm);
  arepng = nms>(13+4);
  arepng(arepng) = cellfun(@(x)all('.png'==x(end-3:end)), nm(arepng));
  arepng(arepng) = cellfun(@(x)all('entornoraster'==lower(x(1:13))), nm(arepng));
  nmpng = nm(arepng);
  [ens ensidx] = sort(cellfun(@(x)str2double(x(14:end-4)), nmpng));
  if not(isempty(ens))
    fprintf('Going (%d) for sim %s\n', thelast, basedir);
    if copysimname
      params = load([basedir filesep '..' filesep 'estado.mat']);
      seps = find(basedir==filesep);
      subs = ['_G' num2str(params.ACIParams.alphaG, '%.2f') '_' basedir(seps(end-1)+1:seps(end)-1) '_'];
      clear params;
    else
      subs = '';
    end
    if copyidx
      idxstr = ['_' mat2str(idx) '.png'];
    else
      idxstr = '';
    end
    if iscell(prefix)
      prefix1 = prefix{1};
      prefix2 = prefix{2};
    else
      prefix1 = prefix;
      prefix2 = '';
    end
    if isnan(thelast) || (numel(ens)<thelast)
      ilast = numel(ens);
    else
      ilast = thelast;
    end
    fprintf('Selected ilast: %s (thelast=%s, numel(ens)=%s)\n', mat2str(ilast), mat2str(thelast), mat2str(numel(ens)));
    last = mat2str(ens(ilast));
    fnm = [basedir filesep nmpng{ensidx(ilast)}];
    for k=1:numel(modes)
      switch modes{k}
        case 'copy'
          fprintf('copying for %s \n', fnm);
          copyfile(fnm, [putdir filesep prefix1 subs last prefix2 idxstr '.png']);
        case 'scale'
          fprintf('loading for %s \n', fnm);
          [img imgmap] = imread(fnm, 'png');
          sz = size(img);
          if sz(1)>sz(2)
            if sz(1)>1000
              fprintf('scaling V for %s \n', fnm);
              [img imgmap] = imresize(img, imgmap, [1000 nan]);
            end
          else
            if sz(2)>1000
              fprintf('scaling H for %s \n', fnm);
              [img imgmap] = imresize(img, imgmap, [nan 1000]);
            end
          end
          fprintf('writing scaled for %s \n', fnm);
          pngname = [putdir filesep prefix1 subs last '_sc' mat2str(sz(2)) 'X' mat2str(sz(1))  idxstr  '.png'];
          imwrite(img, imgmap, pngname, 'png');
          fprintf('wrote scaled to %s \n', pngname);
        otherwise %assume that the string is the name of a file which we must copy
          fil = modes{k};
          doit = true;
          if modes{k}(1)=='#'
            fil = fil(2:end);
            doit = exist([basedir filesep fil], 'file');
          end
          if doit
            fprintf('copying file %s\n', [prefix1 subs last prefix2 fil]);
            copyfile([basedir filesep fil], [putdir filesep prefix1 subs last prefix2 fil]);
          else
            fprintf('file does not exists: %s\n', [basedir filesep fil]);
          end
      end
    end
    idx = idx+1;
  end
end
  
