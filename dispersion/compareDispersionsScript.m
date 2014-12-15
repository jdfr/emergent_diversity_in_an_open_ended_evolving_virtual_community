function compareDispersionsScript(prefix, measurename, pophs, ylim, imgcap, legargs)
%a script to compare the evolution of the mean of the distance matrix over
%time, for a triplet of simulations: 1st very harsh, 2nd harsh, 3rd mild

%compareDispersionsScript('MP_', 'phenotypic diversity', [], {'YLim' [0 1] 'XLim' [0 500] 'layer' 'top'}, true, {'Location', 'SouthEast'})
%compareDispersionsScript('MP_', 'phenotypic diversity', [], {'YLim' [0 1] 'XLim' [0 500] 'layer' 'top'}, true, {'Location', [0.7 0.74310126 0.106770833 0.085443037]})

if not(exist('imgcap', 'var'))
  imgcap = false;
end
if not(exist('legargs', 'var'))
  legargs = {'Location', 'SouthEast'};
end

if isempty(ylim)
  ylim = {};
% else
%   ylim = {'YLim', ylim};
end

%compareDispersionsScript('MP_', 'JACCARD DISTANCE'); compareDispersionsScript('E1_', 'GENOMIC DISTANCE'); compareDispersionsScript('E2_', 'REDUCED GENOMIC DISTANCE');
if ispc
  base = 'lsystemdani\dataset';
  subs = {'G1F1\uno\1=0.001_3', 'G0.75F1\uno\1=0.001_1', 'G0.5F1\dos\1=0.001_1'};
else
  base = 'rene/dataset';
  subs = {'G1F1/uno/1=0.001_3', 'G0.75F1/uno/1=0.001_1', 'G0.5F1/dos/1=0.001_1'};
end
dirs = cellfunc(@(x)[base filesep x], subs);
labels = cellfunc(@(x)x(1:find(x==filesep, 1)), subs);
names = {'\alpha = 1      (very harsh)', '\alpha = 0.75 (harsh)', '\alpha = 0.5   (mild)'};
labs = {'1VH', '2H', '3M'};
widths = repmat({2}, numel(dirs), 1);
colors = {repmat(0.8, 1, 3), repmat(0.5, 1, 3), repmat(0, 1, 3)}';
lineargs = cellfunc(@(x, y) {'LineWidth', x, 'Color', y}, widths, colors);
lidx = [1 3]; largs = cellfunc(@(x, y) {'LineWidth', x, 'Color', y}, {3 2}', colors(lidx));

ylab=measurename;%'mean over all the values of the distance matrix';

tit1 = '';%sprintf('FAST TRACK DIVERSITY FOR %s (daniel''s measure)', measurename);
tit2 = '';%sprintf('FAST TRACK DIVERSITY FOR %s (rene''s measure)', measurename);
baseimg = 'lsystemdani\src\';
base = 'lsystemdani\figurasrene\newer\';
load([base prefix 'matrixDist_simpleMeasures.mat']); 
% [gens daniels renes variances] = getCurves(dirs, prefix);
% save([prefix 'matrixDist_simpleMeasures.mat'], 'gens', 'daniels', 'renes', 'variances');

[h a] = plotVals(gens, daniels, lineargs, names, 'generation', ylab, tit1, legargs, ylim);
if imgcap
  makesnapshot(h, a, 300, [baseimg prefix 'generations_meanMatrixDist']);
%  makesnapshotPDF(h, a, 300, [baseimg prefix 'generations_meanMatrixDist.pdf']);
end

stds = cellfunc(@sqrt, variances);


vals = stds;
leg1 = 'standard deviation';
for k=1:numel(gens)
  h = figure;
  l1 = [gens{k}(:), daniels{k}(:)];
  l2x = [gens{k}(:)'; gens{k}(:)'; nan(size(gens{k}(:)'))];
  l2y = [daniels{k}(:)'+vals{k}(:)';daniels{k}(:)'-vals{k}(:)'; nan(size(gens{k}(:)'))];
  h1 = line(l2x(:), l2y(:),  largs{1}{:});
  h2 = line(l1(:,1), l1(:,2),  largs{2}{:});
  xl = xlabel('generation');
  yl = ylabel(ylab);
  %txtbx = annotation('textbox', [0.65 0.34310126 0.106770833 0.085443037], 'HorizontalAlignment', 'left', 'String', regexprep(names{k}, '[ ]+', ' '), 'LineStyle', 'none', 'FitBoxToText', 'on');
  %zz = [0.13 0.94 0 0];
  zz = [0.23 1 0 0];
  txtbx = annotation('textbox', zz, 'HorizontalAlignment', 'left', 'String', regexprep(names{k}, '[ ]+', ' '), 'LineStyle', 'none', 'FitBoxToText', 'on');
  %tit = title(sprintf('%s in %s environment', measurename, names{k}));
  doForAll(xl, yl, txtbx);
  set(gca, ylim{:});
  a = gca;
  legend([h2 h1], {measurename, leg1}, legargs{:});
  if imgcap
%     makesnapshotPDF(h, a, 300, [baseimg prefix 'generations_meanMatrixDistWithSTD' labs{k} '.pdf']);
    makesnapshot(h, a, 300, [baseimg prefix 'generations_meanMatrixDistWithSTD' labs{k}]);
  end
end
return

% vals = variances;
% leg1 = 'VAR of dissimilarity';
% for k=1:numel(gens)
%   h = figure;
%   l1 = [gens{k}(:), daniels{k}(:)];
%   l2x = [gens{k}(:)'; gens{k}(:)'; nan(size(gens{k}(:)'))];
%   l2y = [daniels{k}(:)'+vals{k}(:)';daniels{k}(:)'-vals{k}(:)'; nan(size(gens{k}(:)'))];
%   line(l2x(:), l2y(:),  largs{1}{:});
%   line(l1(:,1), l1(:,2),  largs{2}{:});
%   xl = xlabel('generations');
%   yl = ylabel(ylab);
%   tit = title(sprintf('%s in %s environment', measurename, names{k}));
%   doForAll(xl, yl, tit);
%   legend({leg1, 'mean dissimilarity'}, 'Location', 'SouthEast');
%   if imgcap; img = imcapture(h, 'all', 600); imwrite(img, [base prefix 'generations_meanMatrixDistWithVAR.png']); end
% end

% h = plotVals(gens, renes, lineargs, names, 'generations', 'mean of the sum of the row with minimal sum of values', tit2);
h = plotVals(gens, variances, lineargs, names, 'generation', 'variance over all the values of the matrix', tit1);
if imgcap; img = imcapture(h, 'all', 600); imwrite(img, [base prefix 'generations_varianceMatrixDist.png']); end
% h = plotVals(gens, cellfunc(@sqrt, variances), lineargs, names, 'generations', 'standard deviation over all the values of the matrix', tit1);

% h = plotVals(daniels, variances, lineargs, names, 'mean over all the values of the matrix', 'variance over all the values of the matrix', tit1);
% h = plotVals(daniels, cellfunc(@sqrt, variances), lineargs, names, 'mean over all the values of the matrix', 'std over all the values of the matrix', tit1);

%[popsizes biomass] = cellfunc(@(x)cellfun(@(y)deal(sum(x.generation==y), sum(x.nbranches(x.generation==y))), array2cell(1:max(x.generation))'), pophs);

% h = plotVals(daniels, popsizes, lineargs, names, 'diversity', 'population size', tit1);
% h = plotVals(daniels, biomass, lineargs, names, 'diversity', 'biomass', tit1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gens daniels renes variances] = getCurves(dirs, prefix)

gens = cell(3,1);
daniels = cell(3,1);
renes = cell(3,1);
variances = cell(3,1);

for k=1:numel(dirs)
  fprintf('Going for %s, %s...\n', prefix, dirs{k});
  genI = 1;

  lastl = getLastLine([dirs{k} filesep 'poblacion.txt']);
  if isempty(lastl)
    error('Mu raro, no hay líneas con caracteres en el fichero!!!!!!');
  end
  genF = sscanf(lastl, '%f', 1);
  if isempty(genF)
    error('Mu raro, mira linea: <%s>!!!!!!', lastl);
  end
  clear lastl;
  gens{k} = genI:genF;
  [daniels{k} renes{k} variances{k}] = superSimpleDispersion(dirs{k}, prefix);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h a] = plotVals(vxs, vys, lineargs, names, xlab, ylab, titl, legargs, ylim)
if not(exist('legargs', 'var'))
  legargs = {};
end
h = figure;
for k=1:numel(vxs)
  line(vxs{k}, vys{k}, lineargs{k}{:});
end
legend(names, legargs{:});
xl = xlabel(xlab);
yl = ylabel(ylab);
a = gca;
labs = {xl, yl};
if not(isempty(titl))
  labs{end+1} = title(titl);
end
doForAll(labs{:});
set(gca, ylim{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doForAll(varargin)
fs = 30;
set(gca, 'FontSize', fs, 'fontname', 'times');
%grid on; %set(gca, 'XGrid', 'on'); 
%set(gca, 'XGrid', 'on'); 
for k=1:numel(varargin)
  set(varargin{k},  'FontSize', fs, 'fontname', 'times');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makesnapshot(h, a, dpi, file)
%makesnapshotIMG(h, a, dpi, file);
makesnapshotPDF(h, a, dpi, file)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makesnapshotIMG(h, a, dpi, file)
z=get(h, 'papersize');
za=[z(2)*1.2 z(1)*0.75]; set(h, 'papersize', za, 'PaperOrientation', 'portrait', 'PaperPosition', [0 0 za]);
img = imcapture(h, 'all', dpi);
imwrite(img, [file '.png']);
close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makesnapshotPDF(h, a, dpi, file)
z=get(h, 'papersize');
za=[z(2)*1.2 z(1)*0.75];
set(h, ...
    'papersize', za, ...
    'PaperOrientation', 'portrait', ...
    'PaperPosition', [0 0 za]);
p = get(a, 'Position');
shifts = [-0.11 -0.05 0.17 0.07]; %[left bottom width height]
set(a, 'Position', p+shifts);
% xl = get(a,'XLabel');
% set(xl, 'units', 'normalized');
% p = get(xl, 'Position');
% p(2) = p(2)+0.1;
% set(xl, 'Position', p);
saveas(h, [file '.pdf']);
close(h);
