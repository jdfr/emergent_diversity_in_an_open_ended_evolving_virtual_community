function allres = showFitnessPerGenRec(imgcap, overwrite, besafe, base)
%b = 'rene/dataset/test'; bd = cellfunc(@(x) [b x], {'110NB', '110NB_old', 'LONG', 'LONGG', 'LONGOLD', 'LONGOLD2', 'OKPRFX'}); allres = showFitnessPerGenRec(false, true, false, bd{k}); end      

%b = 'rene/dataset/test'; bd = cellfunc(@(x) [b x], {'110NB', '110NB_old'}); doForSeveralSimulationsParallel('torque', makedummy({@(x,varargin)makeAllGraphics(x, true), @(x, varargin)showFitnessPerGenRec(true, true, false, x(1:find(x==filesep, 1, 'last')-1))}), {}, bd, 'rene/src', true, []);  

allres = [];

d = dir(base);
dn = {d.name}';
subdirs = ([d.isdir]') & not(strcmp('.', dn)) & not(strcmp('..', dn));
alldirs = cellfunc(@(x) [base filesep x], dn(subdirs));
if any(strcmpi('estado.mat', dn))
  for k=1:numel(alldirs)
    dd = dir(alldirs{k});
    dd = {dd.name}';
    if overwrite || not(any(strcmp('generations_popsize.png', dd)))
      if not(besafe) || any(strcmp('genomeLengthVSwidth_AtLastGen.png', dd))
        poph = loadSimulation(alldirs{k});
        allres = showFitnessPerGen({poph}, imgcap, alldirs{k});
        clear poph;
      end
    end
  end
else
  alls = cell(numel(alldirs),1);
  for k=1:numel(alldirs)
    alls{k} = showFitnessPerGenRec(imgcap, overwrite, besafe, alldirs{k});
  end
  allres = vertcat(allres, alls{:});
end
  


function res = showFitnessPerGen(pophs, imgcap, base)

allstr = {'potential fitness',       'fitnesspot',@(x)x.fitnesspot,                       false; ...
          'competitive pressure',    'pressure',  @(x)x.pressure,                         false; ...
          'fitness',                 'fitness',   @(x)x.fitness,                          false; ...
          'genotype length',         'gnmlen',    @(x)cellfun('prodofsize', x.genome),    true; ...
          'reduced genotype length', 'gnmlenred', @(x)cellfun('prodofsize', x.minGenome), true; ...
          'number of branches',      'nbranches', @(x)x.nbranches,                        true; ...
          'height',                  'height',    @(x)max(x.height, 0.1),                           true; ...
          };

%figures showing fitness in three well defined simulations

%showFitnessPerGen(pophs, names, true, 'lsystemdani\src');
names = {['\alpha=' num2str(pophs{1}.ACIParams.alphaG)]};
labs = {''};

legargs = { {'Location', 'Best'} };
ylim = {'Layer', 'top'};%, 'XLim', [min(pophs{1}.generation) max(pophs{1}.generation)]};

res = struct;

res.alphaG = pophs{1}.ACIParams.alphaG;
res.dir    = base;

widths = {2};
colors = {repmat(0, 1, 3)}';
lineargs = cellfunc(@(x, y) {'LineWidth', x, 'Color', y}, widths, colors);
largs = cellfunc(@(x, y) {'LineWidth', x, 'Color', y}, {2 2}', {repmat(0, 1, 3) repmat(0.7, 1, 3)}');

for q=1:size(allstr,1)
  for k=1:numel(pophs)
    poph = pophs{k};
    poph.fitnesspot = calculateCeilFitness(poph.ACIParams, poph.height, poph.nleafsOnTop, poph.nbranches, poph.nleafs, poph.nPixelsInSoil);
    poph.pressure   = 1 - poph.fitness ./ poph.fitnesspot;
    ming = min(poph.generation);
    maxg = max(poph.generation);
    gens = (ming:maxg)';
    meanf  = zeros(size(gens));
    stdf  = zeros(size(gens));
    popsizes = zeros(size(gens));
    biomass  = zeros(size(gens));
    variables = allstr{q,3}(poph);
    for z=1:numel(gens)
      g = gens(z);
      thisg = find(poph.generation==g);
      thisg = thisg(poph.nPixelsInSoil(thisg)==0);
      fs = variables(thisg);
      meanf(z) = mean(fs);
      stdf(z) = std(fs);
      popsizes(z) = numel(thisg);
      biomass(z) = sum(poph.nbranches(thisg));
    end
    
    if q==1
      res(k).gens = gens;
      res(k).popsize = popsizes;
      res(k).biomass = biomass;
      %if imgcap; img = imcapture(h, 'all', 300); imwrite(img, [base filesep 'generations_popsize.png']); close(h); end
      if imgcap;
        h = plotVals({gens(:)}, {popsizes(:)}, largs, {}, 'generations', 'population size', '', {});
        saveas(h, [base filesep 'generations_popsize.png'], 'png');
        close(h);
      end
      %if imgcap; img = imcapture(h, 'all', 300); imwrite(img, [base filesep 'generations_biomass.png']); close(h); end
      if imgcap; 
        h = plotVals({gens(:)}, {biomass(:)}, largs, {}, 'generations', 'biomass (total amount of branches)', '', {});
        saveas(h, [base filesep 'generations_biomass.png'], 'png');
        close(h);
      end
    end
    
    res(k).([allstr{q,2} 'MEAN']) = meanf(:);
    res(k).([allstr{q,2} 'STD'])  = stdf(:);

    vsa = variables(poph.nPixelsInSoil==0);
    %if imgcap; img = imcapture(h, 'all', 300); imwrite(img, [base filesep 'generations_' allstr{q,2} '_' labs{k} '.png']); close(h); end
    if imgcap;
      h = showLogHist3(poph.generation(poph.nPixelsInSoil==0), 'generation', gens, vsa, allstr{q,1}, 100, [allstr{q,1} ' in ' names{k} ' environment'], [1 1 1; jet(255)], allstr{q,4});
      if not(isempty(ylim))
        set(gca, ylim{:});
      end
      saveas(h, [base filesep 'generations_' allstr{q,2} '_' labs{k} '.png'], 'png');
      close(h);
    end

  %   ratios = popsizes(2:end)./popsizes(1:end-1);
  % 
    ylab = ['mean ' allstr{q,1}];
    vals = stdf;
    leg1 = '\sigma (standard deviation)';
    l1 = [gens(:), meanf(:)];
    l2x = [gens(:)'; gens(:)'; nan(size(gens(:)'))];
    l2y = [meanf(:)'+vals(:)';meanf(:)'-vals(:)'; nan(size(gens(:)'))];
    %if imgcap; img = imcapture(h, 'all', 300); imwrite(img, [base filesep 'generations_mean' allstr{q,2} 'WithSTD_' labs{k} '.png']); close(h);  end
    if imgcap;
      h = figure;
      hstd = line(l2x(:), l2y(:),  largs{2}{:});
      hmn  = line(l1(:,1), l1(:,2),  largs{1}{:});
      xl = xlabel('generations');
      yl = ylabel(ylab);
      tit = title(sprintf('mean %s in %s environment', allstr{q,1}, names{k}));
      doForAll(xl, yl, tit);
      if not(isempty(ylim))
        set(gca, ylim{:});
      end
      %set(gca, 'XLim', [ming maxg]);
      legend([hmn hstd], {['mean ' allstr{q,1}], leg1}, legargs{k}{:});
      saveas(h, [base filesep 'generations_mean' allstr{q,2} 'WithSTD_' labs{k} '.png'], 'png');
      close(h);
    end
  %   figure;
  %   plot(gens(2:end), ratios);
  %   set(gca, 'YLim', [0 2], 'XLim', [1, gens(end)]);
  %   tit = title(sprintf('mean offspring amount in %s environment', names{k}));
  end
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = plotVals(vxs, vys, lineargs, names, xlab, ylab, titl, legargs)
if not(exist('legargs', 'var'))
  legargs = {};
end
h = figure;
for k=1:numel(vxs)
  line(vxs{k}, vys{k}, lineargs{k}{:});
end
if not(isempty(names))
  legend(names, legargs{:});
end
if not(isempty(xlab))
  xl = xlabel(xlab);
  doForAll(xl);
end
if not(isempty(ylab))
  yl = ylabel(ylab);
  doForAll(yl);
end
if not(isempty(titl))
  tit = title(titl);
  doForAll(tit);
end
%doForAll(xl, yl, tit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doForAll(varargin)
% fs = 10;%14;
% set(gca, 'FontSize', fs, 'fontname', 'times');
% %grid on; %set(gca, 'XGrid', 'on'); 
% %set(gca, 'XGrid', 'on'); 
% for k=1:numel(varargin)
%   set(varargin{k},  'FontSize', fs, 'fontname', 'times');
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ceilFitness = calculateCeilFitness(params, height, nleafsOnTop, nbranches, nleafs, negY)
haveNotTooManyNeg  = negY<=params.negYThreshold;
if (~params.tooTallArePlaced) && (~isempty(params.heightT))
  areNotTooTall    = height<=params.heightT;
else
  areNotTooTall    = true(size(haveNotTooManyNeg));
end
haveSomeFitness    = haveNotTooManyNeg & areNotTooTall;
ceilFitness        = zeros(size(haveNotTooManyNeg));
ceilFitness(haveSomeFitness) = nleafsOnTop(haveSomeFitness);
switch params.divideFitness
  case 'byLeafs'
    % dividiendo por el numero de hojas del arbol
    do       = nleafs>0; % |G|= 10
    ceilFitness(do)  = (ceilFitness(do)-negY(do))./((nleafs(do)*params.factorG).^params.alphaG);
    ceilFitness(~do) = 0;
  case 'byBranches'
    % dividiendo por el numero de ramas distinguibles del arbol *|G|
    do               = nbranches>0; % |G|= 10
    ceilFitness(do)  = (ceilFitness(do)-negY(do))./((nbranches(do)*params.factorG).^params.alphaG);
    ceilFitness(~do) = 0;
end
ceilFitness(ceilFitness<0) = -1;
