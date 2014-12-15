function ztrees = loadSampleTreesFromSimulation(basedir, nsamples, generations)

if iscell(basedir)
  poph = basedir{1};
  basedir = basedir{2};
else
  poph = loadSimulation(basedir, false, false, false);
end

ztrees = {'genome', 'g. minimized', 'height', 'width', 'nbranches', 'image', 'generation', 'individual'};

for k=1:numel(generations)
  g = generations(k);
  tg = find(poph.generation==g);
  tg = tg(poph.nPixelsInSoil(tg)==0);
  ns = min(nsamples, numel(tg));
  sel = randperm(numel(tg));
  sel = tg(sel(1:ns));
  gnm = poph.genome(sel);
  gmin = cellfun(@sanitizeString, gnm, 'uniformoutput', false);
  inds = poph.individual(sel);
  v = sprintf('P%03d', g);
  file = [basedir filesep v '.mat'];
  images = cell(size(gnm));
  if exist(file, 'file')
   load(file, v);
   if exist(v, 'var')
     P = eval(v);
     clear(v);
     images = cellfun(@myDraw, P.raster(inds), makecell(P.dimensions(inds,:)), 'uniformoutput', false);
   end
  end
  
  ztrees = [ztrees;
    gnm, gmin, ...
    array2cellMX([poph.height(sel), poph.width(sel), poph.nbranches(sel)]), ...
    images, array2cellMX([poph.generation(sel), inds]) ];
end

function x=makecell(x)
x = mat2cell(x, ones(size(x,1),1), size(x,2));

function array = myDraw(raster, dimensions)
ys = raster{1}; if size(ys,2)>1; ys = ys(:); end;
xs = raster{2}; if size(xs,2)>1; xs = xs(:); end;
%array = accumarray([ys,xs],uint16(indexInColors+2), uint16(dimensions));
array = false(dimensions);
array(sub2ind(dimensions, ys, xs)) = true;
array = array(end:-1:1,:);
if size(array,2)==1
  array(end,2)=false;
end
