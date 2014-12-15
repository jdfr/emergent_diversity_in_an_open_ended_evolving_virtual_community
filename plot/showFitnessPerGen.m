function showFitnessPerGen(pophs, imgcap, base)
%figures showing fitness in three well defined simulations

%showFitnessPerGen(pophs, names, true, 'lsystemdani\src');
names = {'\alpha = 1      (very harsh)', '\alpha = 0.75 (harsh)', '\alpha = 0.5   (mild)'};
%names = {'very harsh', 'harsh', 'mild'};
labs = {'1VH', '2H', '3M'};

if not(exist('legargs', 'var'))
  legargs = { {'Location', 'SouthEast'}, ...
              {'Location', 'SouthEast'}, ...
              {'Location', 'SouthEast'} };
end
if not(exist('ylim', 'var')) %|| isempty(ylim)
  ylim = {...
            {'Layer', 'top', 'YLim', [0 1.0+.4], 'XLim', [0 500]}; ...
            {'Layer', 'top', 'YLim', [0 1.4], 'XLim', [0 500]}; ...
            {'Layer', 'top', 'YLim', [0 1.4], 'XLim', [0 500]}; ...
         };
%   ylim = {};
else
  ylim = {'YLim', ylim};
end


widths = repmat({2}, numel(pophs), 1);
colors = {repmat(0.8, 1, 3), repmat(0.5, 1, 3), repmat(0, 1, 3)}';
lineargs = cellfunc(@(x, y) {'LineWidth', x, 'Color', y}, widths, colors);
lidx = [1 3]; largs = cellfunc(@(x, y) {'LineWidth', x, 'Color', y}, {3 2}', colors(lidx));

genss = cell(size(pophs));
popss = cell(size(pophs));
bioms = cell(size(pophs));

for k=1:numel(pophs)
  poph = pophs{k};
  ming = min(poph.generation);
  maxg = max(poph.generation);
  gens = (ming:maxg)';
  genss{k} = gens;
  meanf  = zeros(size(gens));
  stdf  = zeros(size(gens));
  popsizes = zeros(size(gens));
  biomasses  = zeros(size(gens));
  for z=1:numel(gens)
    g = gens(z);
    thisg = find(poph.generation==g);
    thisg = thisg(poph.nPixelsInSoil(thisg)==0);
    fs = poph.fitness(thisg);
    meanf(z) = mean(fs);
    stdf(z) = std(fs);
    popsizes(z) = numel(thisg);
    biomasses(z) = sum(poph.nbranches(thisg));
  end
  popss{k} = popsizes;
  bioms{k} = biomasses;
  
%   h = showLogHist3(poph.generation(poph.fitness>=0), 'generation', gens, poph.fitness(poph.fitness>=0), 'fitness', linspace(0, 3, 100), ['fitness in ' names{k} ' environment'], [1 1 1; jet(255)], false);
%   if imgcap; img = imcapture(h, 'all', 600); imwrite(img, [base filesep 'generations_fitness_' labs{k} '.png']); end
  continue
  
  ratios = popsizes(2:end)./popsizes(1:end-1);

  ylab = 'mean fitness';
  vals = stdf;
  leg1 = 'standard deviation';
  h = figure;
  l1 = [gens(:), meanf(:)];
  l2x = [gens(:)'; gens(:)'; nan(size(gens(:)'))];
  if strcmp(labs{k}, '1VH')
    fprintf('capping method 1...\n');
    l2y = [min(meanf(:)'+vals(:)', 1);meanf(:)'-vals(:)'; nan(size(gens(:)'))];
  else
    l2y = [meanf(:)'+vals(:)';meanf(:)'-vals(:)'; nan(size(gens(:)'))];
  end
  hstd = line(l2x(:), l2y(:),  largs{1}{:});
  hmn  = line(l1(:,1), l1(:,2),  largs{2}{:});
  xl = xlabel('generation');
  yl = ylabel(ylab);
  txtbx = annotation('textbox', [0.92 0.9 0 0], 'HorizontalAlignment', 'right', 'String', regexprep(names{k}, '[ ]+', ' '), 'LineStyle', 'none', 'FitBoxToText', 'on');
  %tit = title(sprintf('mean fitness in %s environment', names{k}));
  doForAll(xl, yl, txtbx);%, tit);
  if not(isempty(ylim))
    if iscell(ylim{1})
      fprintf('capping method 2...\n');
      set(gca, ylim{k}{:});
    else
      set(gca, ylim{:});
    end
  end
  %set(gca, 'XLim', [ming maxg]);
  hl = legend([hmn hstd], {'mean fitness', leg1}, legargs{k}{:});
  doForAll(hl);
%  if imgcap; img = imcapture(h, 'all', 600); imwrite(img, [base filesep 'generations_meanFitnessWithSTD_' labs{k} '.png']); end
  if imgcap;
    z=get(h, 'papersize');
    za=[z(2)*1.5 z(1)*0.75]; set(h, 'papersize', za, 'PaperOrientation', 'portrait', 'PaperPosition', [0 0 za]);
    img = imcapture(h, 'all', 400);
    %figure; imshow(img);
    imwrite(img, [base filesep 'generations_meanFitnessWithSTD_' labs{k} '.png']);
  end
%   figure;
%   plot(gens(2:end), ratios);
%   set(gca, 'YLim', [0 2], 'XLim', [1, gens(end)]);
%   tit = title(sprintf('mean offspring amount in %s environment', names{k}));
end

widths = repmat({4}, numel(pophs), 1);
lineargs = cellfunc(@(x, y) {'LineWidth', x, 'Color', y}, widths, colors);
ylab = 'population size';
xlab = 'generation';
titl = '';
legargs{1} = {'Location', 'SouthEast'};
ylim = {'XLim' [0 500]};% 'layer' 'top'};% 'YLim' [0 1] };
h = plotVals(genss, popss, lineargs, names, xlab, ylab, titl, legargs{1}, ylim);
if imgcap;
  z=get(h, 'papersize');
  za=[z(2)*1.5 z(1)*1]; set(h, 'papersize', za, 'PaperOrientation', 'portrait', 'PaperPosition', [0 0 za]);
  img = imcapture(h, 'all', 300);
  %figure; imshow(img);
  imwrite(img, [base filesep 'generations_zpopsize.png']);
close(h);
end

ylab = 'biomass';
ylim = {'XLim' [0 500] 'layer' 'top' 'YLim' [0 860000] 'YTick' (0:4:86)*1e4};
%ylim = {'XLim' [0 500] 'layer' 'top' 'YLim' [0 860000]  'YTick' [1 10 100 1e3 1e4 1e5] 'YScale', 'log'};
%ylim = {'XLim' [0 500] 'layer' 'top' 'YLim' [0 1.2e5]  };
%legargs{1} = {'Location', 'SouthEast'};
legargs{1} = {'Location', 'NorthEast'};
h = plotVals(genss, bioms, lineargs, names, xlab, ylab, titl, legargs{1}, ylim);
if imgcap;
  z=get(h, 'papersize');
  za=[z(2)*1.5 z(1)*1*5];
  %za=[z(2)*1.5 z(1)*1];
  set(h, 'papersize', za, 'PaperOrientation', 'portrait', 'PaperPosition', [0 0 za]);
  img = imcapture(h, 'all', 300);
  %figure; imshow(img);
  imwrite(img, [base filesep 'generations_zbiomass.png']);
close(h);
end

ylab = 'biomass'; xlab = 'population size';
ylim = {'layer' 'top' 'YLim' [0 860000] 'YTick' (0:4:86)*1e4};
%ylim = {'layer' 'top' 'YLim' [0 860000] 'YTick' [1 10 100 1e3 1e4 1e5] 'YScale', 'log'};
%ylim = {'XLim' [0 500] 'layer' 'top' 'YLim' [0 1.2e5]  };
%legargs{1} = {'Location', 'SouthEast'};
legargs{1} = {'Location', 'NorthEast'};
h = plotVals(popss, bioms, lineargs, names, xlab, ylab, titl, legargs{1}, ylim);
if imgcap;
  z=get(h, 'papersize');
  za=[z(2)*1.5 z(1)*1*5];
  %za=[z(2)*1.5 z(1)*1];
  set(h, 'papersize', za, 'PaperOrientation', 'portrait', 'PaperPosition', [0 0 za]);
  img = imcapture(h, 'all', 300);
  %figure; imshow(img);
  imwrite(img, [base filesep 'gpopsize_zbiomass.png']);
close(h);
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function h = plotVals(vxs, vys, lineargs, names, xlab, ylab, titl, legargs)
% if not(exist('legargs', 'var'))
%   legargs = {};
% end
% h = figure;
% for k=1:numel(vxs)
%   line(vxs{k}, vys{k}, lineargs{k}{:});
% end
% legend(names, legargs{:});
% xl = xlabel(xlab);
% yl = ylabel(ylab);
% tit = title(titl);
% doForAll(xl, yl, tit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = plotVals(vxs, vys, lineargs, names, xlab, ylab, titl, legargs, ylim)
if not(exist('legargs', 'var'))
  legargs = {};
end
h = figure;
for k=1:numel(vxs)
  line(vxs{k}, vys{k}, lineargs{k}{:});
end
hl = legend(names, legargs{:});
xl = xlabel(xlab);
yl = ylabel(ylab);
labs = {xl, yl, hl};
if not(isempty(titl))
  labs{end+1} = title(titl);
end
doForAll(labs{:});
if not(isempty(ylim))
  set(gca, ylim{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doForAll(varargin)
fs = 30;
set(gca, 'FontSize', fs, 'fontname', 'times');
%grid on; %set(gca, 'XGrid', 'on'); 
%set(gca, 'XGrid', 'on'); 
for k=1:numel(varargin)
  set(varargin{k},  'FontSize', fs, 'fontname', 'times');
end
