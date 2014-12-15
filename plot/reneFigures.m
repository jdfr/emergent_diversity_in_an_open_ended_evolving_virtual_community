function reneFigures
%script to gather data and generate figures for Rene. Not to be executed as
%presented, but to execute the chunks of lines I need at each moment

base = 'lsystemdani\dataset';
subs = {'G1F1\uno\1=0.001_3', 'G0.75F1\uno\1=0.001_1', 'G0.5F1\dos\1=0.001_1'};
labels = {'very harsh', 'harsh', 'mild'};
labs = {'1VH', '2H', '3M'};
dirs = cellfunc(@(x)[base filesep x], subs);
pophs = cellfunc(@(x)loadSimulation(x, false), dirs);
pophs = getPixels(dirs, pophs);
res = load('lsystemdani\src\MP_matrixDist_ALL.mat');
res = res.res;
dists_mdscale = load('lsystemdani\src\MP_monster_mdscale.mat');
dists_mdscale = dists_mdscale.dists_mdscale;
pophs = interpretMonsterDists(res, dists_mdscale.sstress_mdscaleY, pophs); close all hidden;

lastGens = cellfun(@(x)max(x.generation), pophs);
lastPNames = arrayfunc(@(x)sprintf('P%03g', x), lastGens);
lastPs = cellfunc(@(x, y)load([x filesep y '.mat']), basedirs, lastPNames);

% img = showTreesAndComplexity(pophs{2}, lastGens(2), lastPs{2}, pophs{2}.compression9, false);
% img = showTreesAndComplexity(pophs{2}, lastGens(2), lastPs{2}, pophs{2}.nbranches, false);
% img = showTreesAndComplexity(pophs{2}, lastGens(2), lastPs{2}, pophs{2}.npixels, false);

for k=1:numel(pophs); h = showLogHist3(pophs{k}.generation(pophs{k}.nPixelsInSoil==0), 'generation', [], pophs{k}.mdscale_sstress(pophs{k}.nPixelsInSoil==0), 'MDS metric', 100, labels{k}, [1 1 1; jet(255)], false); postpros(h); img = imcapture(h, 'all', 600); imwrite(img, ['generation_MDSmetric_' labs{k} '.png']); end
for k=1:numel(pophs); h = showLogHist3(pophs{k}.generation(pophs{k}.nPixelsInSoil==0), 'generation', [], pophs{k}.compression9(pophs{k}.nPixelsInSoil==0), 'compression', 100, labels{k}, [1 1 1; jet(255)], true); postpros(h); img = imcapture(h, 'all', 600); imwrite(img, ['generation_compression_' labs{k} '.png']); end
for k=1:numel(pophs); h = showLogHist3(pophs{k}.nbranches(pophs{k}.nPixelsInSoil==0), 'nbranches', 100, pophs{k}.compression9(pophs{k}.nPixelsInSoil==0), 'compression', 100, labels{k}, [1 1 1; jet(255)], false); postpros(h); img = imcapture(h, 'all', 600); imwrite(img, ['nbranches_compression_' labs{k} '.png']); end
for k=1:numel(pophs); h = showLogHist3(pophs{k}.nbranches(pophs{k}.nPixelsInSoil==0), 'nbranches', 100, pophs{k}.pixel(pophs{k}.nPixelsInSoil==0), 'npixels', 100, labels{k}, [1 1 1; jet(255)], false); postpros(h); img = imcapture(h, 'all', 600); imwrite(img, ['nbranches_npixels_' labs{k} '.png']); end
for k=1:numel(pophs); h = showLogHist3(pophs{k}.pixel(pophs{k}.nPixelsInSoil==0), 'npixels', 100, pophs{k}.compression9(pophs{k}.nPixelsInSoil==0), 'compression', 100, labels{k}, [1 1 1; jet(255)], false); postpros(h); img = imcapture(h, 'all', 600); imwrite(img, ['npixels_compression_' labs{k} '.png']); end
for k=1:numel(pophs); h = showLogHist3(pophs{k}.generation(pophs{k}.nPixelsInSoil==0), 'generation', [], pophs{k}.pixel(pophs{k}.nPixelsInSoil==0), 'npixels', 100, labels{k}, [1 1 1; jet(255)], true); postpros(h); img = imcapture(h, 'all', 600); imwrite(img, ['generation_npixels_' labs{k} '.png']); end
for k=1:numel(pophs); h = showLogHist3(pophs{k}.generation(pophs{k}.nPixelsInSoil==0), 'generation', [], pophs{k}.nbranches(pophs{k}.nPixelsInSoil==0), 'nbranches', 100, labels{k}, [1 1 1; jet(255)], true); postpros(h); img = imcapture(h, 'all', 600); imwrite(img, ['generation_nbranches_' labs{k} '.png']); end


imgnames = {...
  'lsystemdani\dataset\G0.5F1\dos\1=0.001_1\entornoRaster100.png';...
  'lsystemdani\dataset\G0.75F1\uno\1=0.001_1\entornoRaster100.png';...
  'lsystemdani\dataset\G1F1\uno\1=0.001_3\entornoRaster100.png';...
  };
imgwrite = 'lsystemdani\src\stack100.png';
sep = 100;
stackEnvironmentsTrueScale(imgnames, imgwrite, sep)

imgnames = {...
  'lsystemdani\dataset\G0.75F1\uno\1=0.001_1\entornoRaster340.png';...
  'lsystemdani\dataset\G1F1\uno\1=0.001_3\entornoRaster340.png';...
  };
imgwrite = 'lsystemdani\src\stack340.png';
sep = 100;
stackEnvironmentsTrueScale(imgnames, imgwrite, sep)

compareDispersionsScript('E2_', 'genotypic diversity', [], {'YLim' [0 60] 'XLim' [0 500] 'layer' 'top'}, true, {'Location', [0.65 0.57310126 0.106770833 0.085443037]});
compareDispersionsScript('MP_', 'phenotypic diversity', [], {'YLim' [0 1] 'XLim' [0 500] 'layer' 'top'}, true, {'Location', 'SouthEast'});


% function compareDispersionsScript(idx)
% %compareDispersionsScript('MP_', 'JACCARD DISTANCE'); compareDispersionsScript('E1_', 'GENOMIC DISTANCE'); compareDispersionsScript('E2_', 'REDUCED GENOMIC DISTANCE');
% if ispc
%   base = 'lsystemdani\dataset';
%   subs = {'G1F1\uno\1=0.001_3', 'G0.75F1\uno\1=0.001_1', 'G0.5F1\dos\1=0.001_1'};
% else
%   base = 'rene/dataset';
%   subs = {'G1F1/uno/1=0.001_3', 'G0.75F1/uno/1=0.001_1', 'G0.5F1/dos/1=0.001_1'};
% end
% dirs = cellfunc(@(x)[base filesep x], subs);
% labels = cellfunc(@(x)x(1:find(x==filesep, 1)), subs);
% legendnames = {'very harsh', 'harsh', 'mild'};
% legendpos = {'SouthEast', 'SouthEast', 'NorthEast'};
% widths = repmat({2}, numel(dirs), 1);
% colors = {repmat(0.8, 1, 3), repmat(0.5, 1, 3), repmat(0, 1, 3)}';
% lineargs = cellfunc(@(x, y) {'LineWidth', x, 'Color', y}, widths, colors);
% 
% pophs = cellfunc(@(x)loadSimulation([basedir filesep 'poph.mat'], false), dirs);
% pophs = cellfunc(@(x)x.poph, pophs);
% gens = cell(size(pophs));
% popsizes = cell(size(pophs));
% biomass = cell(size(pophs));
% for k=1:numel(pophs)
%   gens{k} = (0:max(pophs{k}.generation))';
%   popsizes{k} = zeros(size(gens{k}));
%   biomass{k} = zeros(size(gens{k}));
%   for z=1:numel(gens{k})
%     g = gens{k}(z);
%     thisg = pophs{k}.generation==g;
%     thisg(thisg) = pophs{k}.nPixelsInSoil(thisg)==0;
%     popsizes{k}(z) = sum(thisg);
%     biomass{k}(z) = sum(pophs{k}.nbranches(thisg);
%   end
% end
%     
% measures = {'MP_', 'JACCARD DISTANCE'; 'E1_', 'GENOMIC DISTANCE'; 'E2_', 'REDUCED GENOMIC DISTANCE'};
% prefixes = measures(:,1);
% measurenames = measures(:,2);
% 
% visible = false;
% [gens daniels renes] = getCurves(dirs, prefixes{idx});
% h = plotVals(visible, gens, daniels, lineargs, legendnames, legendpos{idx}, 'generations', 'mean over all the values of the distance matrix');
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [gens daniels renes] = getCurves(dirs, prefix)
% 
% gens = cell(3,1);
% daniels = cell(3,1);
% renes = cell(3,1);
% 
% for k=1:numel(dirs)
%   fprintf('Going for %s, %s...\n', prefix, dirs{k});
%   genI = 1;
% 
%   lastl = getLastLine([dirs{k} filesep 'poblacion.txt']);
%   if isempty(lastl)
%     error('Mu raro, no hay líneas con caracteres en el fichero!!!!!!');
%   end
%   genF = sscanf(lastl, '%f', 1);
%   if isempty(genF)
%     error('Mu raro, mira linea: <%s>!!!!!!', lastl);
%   end
%   clear lastl;
%   gens{k} = genI:genF;
%   [daniels{k} renes{k}] = superSimpleDispersion(dirs{k}, prefix);
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function h = plotVals(visible, vxs, vys, lineargs, legendnames, legendpos, xlab, ylab, titl)
% h = figure('Visible', visible);
% for k=1:numel(vxs)
%   line(vxs{k}, vys{k}, lineargs{k}{:});
% end
% legend(legendnames, legendpos);
% xl = xlabel(xlab);
% yl = ylabel(ylab);
% tit = title(titl);
% doForAll(xl, yl, tit);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function doForAll(varargin)
% fs = 14;
% set(gca, 'FontSize', fs, 'fontname', 'times');
% %grid on; %set(gca, 'XGrid', 'on'); 
% %set(gca, 'XGrid', 'on'); 
% for k=1:numel(varargin)
%   set(varargin{k},  'FontSize', fs);
% end
