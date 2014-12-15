function res = analyzeGenomeReduced(poph, res, doShow)

if nargin<2 || isempty(res)
  if ischar(poph)
    if exist([poph filesep 'poph.mat'], 'file')
      poph = load([poph filesep 'poph.mat']);
      poph = poph.poph;
    else
      poph = loadSimulation(poph);
    end
  end

  gs = poph.generation;
  gn = poph.genome;

  if not(isfield(poph, 'minGenome'))
    fprintf('CALCULATING MINGENOME...\n');
    mn = cell(size(gn));
    g = nan;
    for k=1:numel(mn)
      if gs(k)~=g
        g = gs(k);
        if mod(g, 10)==0
          fprintf('GENERATION %d\n', g);
        end
      end
      mn{k} = sanitizeString(gn{k});
    end
    poph.minGenome = mn;
    fprintf('MINGENOME CALCULATED...\n');
  else
    mn = poph.minGenome;
  end

  fprintf('CALCULATING T3...\n');
  GT3  = false(size(gn));
  MGT3 = false(size(gn));
  for k=1:numel(mn)
    GT3(k)  = strIsT3(gn{k});
    MGT3(k) = strIsT3(mn{k});
  end
  fprintf('T3 CALCULATED\n');
  fprintf('CALCULATING SIZES...\n');
  Gs  = cellfun('prodofsize', gn);
  MGs = cellfun('prodofsize', mn);
  ratio = MGs./Gs;
  fprintf('SIZES CALCULATED\n');

  ming      = min(gs);
  maxg      = max(gs);
  gens      = (ming:maxg)';
  meanC     = zeros(size(gens));
  medianC   = zeros(size(gens));
  minC      = zeros(size(gens));
  maxC      = zeros(size(gens));
  stdC      = zeros(size(gens));
  meanG     = zeros(size(gens));
  medianG   = zeros(size(gens));
  minG      = zeros(size(gens));
  maxG      = zeros(size(gens));
  stdG      = zeros(size(gens));
  meanMG    = zeros(size(gens));
  medianMG  = zeros(size(gens));
  minMG     = zeros(size(gens));
  maxMG     = zeros(size(gens));
  stdMG     = zeros(size(gens));
  ratioGT3  = zeros(size(gens));
  ratioMGT3 = zeros(size(gens));

  for k=1:numel(gens)
    g     = gens(k);
    thisG = find(gs==g);
    Gsz   =  Gs(thisG);
    MGsz  = MGs(thisG);
    T3G   =  GT3(thisG);
    T3MG  = MGT3(thisG);
    rt    = ratio(thisG);
    ratioGT3(k)  = sum(T3G)/numel(T3G);
    ratioMGT3(k) = sum(T3MG)/numel(T3MG);
    minC(k)      =  min(rt);
    maxC(k)      =  max(rt);
    meanC(k)     = mean(rt);
    medianC(k)   = median(rt);
    stdC(k)      =  std(rt);
    minG(k)      =  min(Gsz);
    maxG(k)      =  max(Gsz);
    meanG(k)     = mean(Gsz);
    medianG(k)   = median(Gsz);
    stdG(k)      =  std(Gsz);
    minMG(k)     =  min(MGsz);
    maxMG(k)     =  max(MGsz);
    meanMG(k)    = mean(MGsz);
    medianMG(k)  = median(MGsz);
    stdMG(k)     =  std(MGsz);
  end

  res = struct('generations', gens, ... 
               'nonminT3', ratioGT3, 'minT3', ratioMGT3, ...
               'compression', struct('mean', meanC,        'median', medianC,        'std', stdC,        'min', minC,        'max', maxC), ...
               'nominimized', struct('mean', meanG,        'median', medianG,        'std', stdG,        'min', minG,        'max', maxG), ...
               'minimized',   struct('mean', meanMG,       'median', medianMG,       'std', stdMG,       'min', minMG,       'max', maxMG));
end

if nargin<3 || not(doShow)
  return
end
  
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
x=axes;
set(x,'FontName','times','FontSize',12)
hold on
  stairs(barsx, barsy, 'Color', repmat(0.9, 1, 3), 'LineWidth', 2);
  stairs(res.generations, res.nominimized.mean,  'Color', repmat(0.1, 1, 3), 'LineWidth', 1);
%   stairs(res.generations, res.minimized.mean,  'MarkerSize', 1, 'Color', 'k', 'LineStyle', 'none', 'Marker', '.');
  plot(res.generations, res.minimized.mean,  'k:');
  legend({'std of genome length','average genome length', 'average useful genome length'},'Location','NorthWest');
  xlabel('generations')
  ylabel('genome length')

  
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
