function res = makePacoGenomeGraph(poph, res)
if nargin<2
  tic;
  genomeLengths = cellfun('prodofsize', poph.genome);
  toc, tic;
  minimized = cell(size(genomeLengths));
  genome = poph.genome;
  nchars = 0;
  nm = numel(minimized);
  for k=1:nm
    minimized{k}       = sanitizeString(genome{k}, false);
    if (mod(k, 100)==0)
      fprintf(repmat('\b', 1, nchars));
      nchars = fprintf('max: %06d, k==%06d, etime: %g', nm, k, toc);
    end
  end
  fprintf(repmat('\b', 1, nchars));
  toc, tic;

  minimizedLength = cellfun('prodofsize', minimized);
  
  ming = min(poph.generation);
  maxg = max(poph.generation);

  generations = (ming:maxg)';

  compressmedian  = zeros(size(generations));
  compressmean    = zeros(size(generations));
  compressstd     = zeros(size(generations));
  compressmin     = zeros(size(generations));
  compressmax     = zeros(size(generations));
  medians  = zeros(size(generations));
  means    = zeros(size(generations));
  stds     = zeros(size(generations));
  mins     = zeros(size(generations));
  maxs     = zeros(size(generations));
  mediansm = zeros(size(generations));
  meansm   = zeros(size(generations));
  stdsm    = zeros(size(generations));
  minsm    = zeros(size(generations));
  maxsm    = zeros(size(generations));
  toc, tic;

  for k=1:numel(generations);
    gs          = poph.generation==generations(k);
    gl          = genomeLengths(gs);
    gm          = minimizedLength(gs);
    compress    = gm./gl;
    compressmedian(k) = median(compress);
    compressmean(k)   = mean(compress);
    compressstd(k)    = std(compress);
    compressmin(k)    = min(compress);
    compressmax(k)    = max(compress);
    medians(k)  = median(gl);
    means(k)    = mean(gl);
    stds(k)     = std(gl);
    mins(k)     = min(gl);
    maxs(k)     = max(gl);
    mediansm(k) = median(gm);
    meansm(k)   = mean(gm);
    stdsm(k)    = std(gm);
    minsm(k)    = min(gm);
    maxsm(k)    = max(gm);
  end
  toc

  res = struct('generations', generations, 'minimizedGenomes', {minimized}, ... 
               'compression', struct('mean', compressmean, 'median', compressmedian, 'std', compressstd, 'min', compressmin, 'max', compressmax), ...
               'nominimized', struct('mean', means,        'median', medians,        'std', stds,        'min', mins,        'max', maxs), ...
               'minimized',   struct('mean', meansm,       'median', mediansm,       'std', stdsm,       'min', minsm,       'max', maxsm));

  
else
  nani       = nan(size(res.generations, 2), size(res.generations, 1));
  barsx      = reshape([res.generations';                            res.generations';                            nani], [], 1);
  barsy      = reshape([(res.nominimized.mean+res.nominimized.std)'; (res.nominimized.mean-res.nominimized.std)'; nani], [], 1);
  barsym     = reshape([(res.minimized.mean  +res.minimized.std)';   (res.minimized.mean  -res.minimized.std)';   nani], [], 1);
  barsyc     = reshape([(res.compression.mean+res.compression.std)'; (res.compression.mean-res.compression.std)'; nani], [], 1);
  barargs    = {'Color', 'r', 'LineWidth', 2};
  plotargs   = {'Color', 'b', 'LineWidth', 1};
  minargs    = {'Color', 'c', 'LineWidth', 1};%, 'LineStyle', '--'};
  maxargs    = {'Color', 'm', 'LineWidth', 1};%, 'LineStyle', '--'};
  medianargs = {'Color', 'k', 'LineWidth', 1};%, 'LineStyle', '--'};
  
  figure('Color', 'w');
  hold on;
%   axis off;
  stairs(barsx, barsy, 'Color', repmat(0.9, 1, 3), 'LineWidth', 2);
  stairs(res.generations, res.nominimized.mean,  'Color', repmat(0.1, 1, 3), 'LineWidth', 1);
%   stairs(res.generations, res.minimized.mean,  'MarkerSize', 1, 'Color', 'k', 'LineStyle', 'none', 'Marker', '.');
  plot(res.generations, res.minimized.mean,  'k:');
  legend({'mean +/- standard deviation', 'mean', 'minimized''s mean'});

  
%   showplot(barsx, barsy,  barargs, res.generations, res.nominimized.mean, plotargs, res.nominimized.min, minargs, res.nominimized.max, maxargs, res.nominimized.median, medianargs);
%   title('genome length mean and spread (mean +/- std)');
%   
%   showplot(barsx, barsym, barargs, res.generations, res.minimized.mean,   plotargs, res.minimized.min,   minargs, res.minimized.max,   maxargs, res.minimized.median,   medianargs);
%   title('minimized genome length mean and spread (mean +/- std)');
%   
%   showplot(barsx, barsyc, barargs, res.generations, res.compression.mean, plotargs, res.compression.min, minargs, res.compression.max, maxargs, res.compression.median, medianargs);
%   title('ratio (minimized/non-minimized) mean and spread (mean +/- std)');
%   
%   figure; hold on;
%     stairs(res.generations, res.nominimized.mean,  'Color', 'r');
%     stairs(res.generations, res.minimized.mean,    'Color', 'b');
%     legend({'mean genome length', 'mean minimized'});
%   hold off; grid on;
  
%   len   = cellfun('prodofsize', poph.genome);
%   lenm  = cellfun('prodofsize', res.minimizedGenomes);
%   ming  = min(poph.generation);
%   maxg  = max(poph.generation);
%   maxl  = max(len);
%   maxlm = max(lenm);
%   
%   rangesx  = [ming maxg];
%   rangesy  = [0 maxl];
%   rangesym = [0 maxlm];%maxl];
%   rangesyh = [0 500];%max(poph.height)];%maxm];%maxl];
%   rangesyb = [0 1000];%max(poph.nbranches)];%maxm];%maxl];
%   
%   histo  = hist3([poph.generation, len],  {rangesx(1):rangesx(2), rangesy(1):rangesy(2)});
%   histom = hist3([poph.generation, lenm], {rangesx(1):rangesx(2), rangesym(1):rangesym(2)});
%   histoh  = hist3([poph.generation, poph.height],     {rangesx(1):rangesx(2), rangesyh(1):rangesyh(2)});
%   histob  = hist3([poph.generation, poph.nbranches],  {rangesx(1):rangesx(2), rangesyb(1):rangesyb(2)});
%   histomh  = hist3([lenm, poph.height],  {rangesym(1):rangesym(2), rangesyh(1):rangesyh(2)});
%   histomb  = hist3([lenm, poph.nbranches],  {rangesym(1):rangesym(2), rangesyb(1):rangesyb(2)});
%   scatterargs = {'Color', 'k', 'LineWidth', 2};
% 
%   allrangex = rangesx(1):rangesx(2);
%   allrangem = rangesym(1):rangesym(2);
%   allrangeh = rangesyh(1):rangesyh(2);
%   allrangeb = rangesyb(1):rangesyb(2);
%   
%   limitm = 101; limith = 101; limitb = 100;
%   
%   scatterplot(histoh,  rangesx,  rangesyh, res.generations, [], 'generations',             'tree height',        'amount of individuals', scatterargs);
%   title('amount of individuals by generation and tree height'); 
% %   scatterplot(histob,  rangesx,  rangesyb, res.generations, [], 'generations',             'number of branches', 'amount of individuals', scatterargs);
%   scatterplot(histomh(1:limitm, 1:limith), allrangem([1 limitm]), allrangeh([1 limith]), res.generations, [], 'minimized genome length', 'tree height',        'amount of individuals', scatterargs);
%   title('amount of individuals by minimized genome length and tree height'); 
%   scatterplot(histomb(1:limitm, 1:limitb), allrangem([1 limitm]), allrangeb([1 limitb]), res.generations, [], 'minimized genome length', 'number of branches', 'amount of individuals', scatterargs);
%   title('amount of individuals by minimized genome length and number of branches'); 
%   
%   scatterplot(histo,  rangesx, rangesy,  res.generations, res.nominimized.mean, 'generations', 'genome length',           'amount of individuals', scatterargs);
%   title('amount of individuals by generation and genome length, and mean genome length');
%   scatterplot(histom, rangesx, rangesym, res.generations, res.minimized.mean,   'generations', 'minimized genome length', 'amount of individuals', scatterargs);
%   title('amount of individuals by generation and genome length, and mean minimized genome length');
%   
%   scatterplot(log10(histoh),  rangesx,  rangesyh, res.generations, [], 'generations',             'tree height',        'log10( amount of individuals )', scatterargs);
%   title('log10( amount of individuals ) by generation and tree height'); 
% %   scatterplot(histob,  rangesx,  rangesyb, res.generations, [], 'generations',             'number of branches', 'amount of individuals', scatterargs);
%   scatterplot(log10(histomh(1:limitm, 1:limith)), allrangem([1 limitm]), allrangeh([1 limith]), res.generations, [], 'minimized genome length', 'tree height',        'log10( amount of individuals )', scatterargs);
%   title('log10( amount of individuals ) by minimized genome length and tree height'); 
%   scatterplot(log10(histomb(1:limitm, 1:limitb)), allrangem([1 limitm]), allrangeb([1 limitb]), res.generations, [], 'minimized genome length', 'number of branches', 'amount of individuals', scatterargs);
%   title('log10( amount of individuals ) by minimized genome length and number of branches'); 
%   
%   scatterplot(log10(histo),  rangesx, rangesy,  res.generations, res.nominimized.mean, 'generations', 'genome length',           'log10( amount of individuals )', scatterargs);
%   title('log10( amount of individuals ) by generation and genome length, and mean genome length');
%   scatterplot(log10(histom), rangesx, rangesym, res.generations, res.minimized.mean,   'generations', 'minimized genome length', 'log10( amount of individuals )', scatterargs);
%   title('log10( amount of individuals ) by generation and minimized genome length, and mean minimized genome length');
%   
%   
end

function scatterplot(histo, rx, ry, x, y, labelx, labely, labelc, scatterargs)
figure;
imagesc(rx, ry, histo');
axis xy;
hold on;
if (~isempty(x)) && (~isempty(y))
  line(x, y, scatterargs{:});
end
xlabel(labelx);
ylabel(labely);
cb = colorbar;
ylabel(cb, labelc);

function showplot(barx, bary, barargs, x, y, plotargs, miny, minargs, maxy, maxargs, mediany, medianargs)
figure;
stairs(barx, bary, barargs{:});
hold on;
stairs(x, mediany, medianargs{:});
stairs(x, miny,    minargs{:});
stairs(x, maxy,    maxargs{:});
stairs(x, y,       plotargs{:});
hold off;
grid on;
legend({'mean +/- standard deviation', 'median', 'min', 'max', 'mean'});
%legend({'standard deviation', 'min', 'max', 'mean'});
