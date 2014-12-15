function generateDispersionPNG(basedir, recursive)
%b = 'rene/dataset/test'; bd = cellfunc(@(x) [b x], {'110NB', '110NB_old', 'LONG', 'LONGG', 'LONGOLD', 'LONGOLD2', 'OKPRFX'}); for k=1:numel(bd); generateDispersionPNG(bd{k}, true); end   

if ~exist('recursive', 'var')
  recursive = false;
end

imgcap = true;
ylim = struct(...
          'pht', {{'layer' 'top'}}, ...
          'gnm', {{'layer' 'top'}}, ...
          'mng', {{'layer' 'top'}} ...
          );
% ylim = struct(...
%           'pht', {{'YLim', [0 1],   'XLim', [0 1000], 'layer' 'top'}}, ...
%           'gnm', {{'YLim', [0 500], 'XLim', [0 1000], 'layer' 'top'}}, ...
%           'mng', {{'YLim', [0 100], 'XLim', [0 1000], 'layer' 'top'}} ...
%           );
measurename = struct(...
          'pht', 'phenotypic diversity', ...
          'gnm', 'genotypic diversity (full)', ...
          'mng', 'genotypic diversity (red.)' ...
          );
legargs = {'Location', 'Best'};%'SouthEast'};

doForSeveralSimulations(@generateDispersionStatsAux, basedir, [], {imgcap, ylim, measurename, legargs}, recursive);

function data = generateDispersionStatsAux(basedir, data, imgcap, ylims, measurenames, legargs)

filename = 'dispersionsss.mat';

if not(exist([basedir filesep filename], 'file'))
  fprintf('   DISPERSIONS.MAT NOT FOUND!!!\n');
  return
else
  fprintf('   PROCESSING...\n');
end

f2 = ['..' filesep 'estado.mat'];

d = load([basedir filesep filename]); 
d = d.data;

p = load([basedir filesep f2]);
p = p.ACIParams;

r3 = @(x) [x x x];

largs = cellfunc(@(x, y) {'LineWidth', x, 'Color', y}, {3 2}', {r3(0.8), r3(0)}');

name = ['\alpha = ' mat2str(p.alphaG)];

gens = d.gens;

types = fieldnames(d.dispersion);

for z=1:numel(types)
  means = d.dispersion.(types{z}).mean;
  stds  = realsqrt(d.dispersion.(types{z}).var);
  ylab  = measurenames.(types{z});
  ylim  = ylims.(types{z});

  vals = stds;
  legs = {'mean', 'standard deviation'};
    h = figure;
    l1 = [gens(:), means(:)];
    l2x = [gens(:)'; gens(:)'; nan(size(gens(:)'))];
    l2y = [means(:)'+vals(:)';means(:)'-vals(:)'; nan(size(gens(:)'))];
    h1 = line(l2x(:), l2y(:),  largs{1}{:});
    h2 = line(l1(:,1), l1(:,2),  largs{2}{:});
    xl = xlabel('generation');
    yl = ylabel(ylab);
    txtbx = title(name);
    doForAll(xl, yl, txtbx);
    set(gca, ylim{:});
    legend([h2 h1], legs, legargs{:});
    if imgcap; saveas(h, [basedir filesep 'generations_meanMatrixDistWithSTD_' types{z} '.png'], 'png'); end
%     if imgcap
%       img = makesnapshot(h, 150);
%       imwrite(img, [basedir filesep 'generations_meanMatrixDistWithSTD_' types{z} '.png']);
%       close(h);
%     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doForAll(varargin)
% fs = 14;
% set(gca, 'FontSize', fs, 'fontname', 'times');
% %grid on; %set(gca, 'XGrid', 'on'); 
% %set(gca, 'XGrid', 'on'); 
% for k=1:numel(varargin)
%   set(varargin{k},  'FontSize', fs, 'fontname', 'times');
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function img = makesnapshot(h, dpi)
z=get(h, 'papersize');
za=[z(2)*1.2 z(1)*0.75]; set(h, 'papersize', za, 'PaperOrientation', 'portrait', 'PaperPosition', [0 0 za]);
img = imcapture(h, 'all', dpi);
