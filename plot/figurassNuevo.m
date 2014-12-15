function hs=figurassNuevo(figuras, toSelect, clustType, useDefaultColors, shows, ejes, ticks, fname, wpos)
%ejes: {[minPop maxPop] [minSpecs maxSpecs] [minBio maxBio]}  
%ticks: lo mesmo con ticks

%load(figurass_6sims)
%figurass(figuras, [6, 2, 5], true, {[0 800], [0 150], [0 3.6e4]}, {0:100:800, 0:20:120, 0:0.5e4:3.5e4}, 'PHASEPLOT_%s', [1000 750]);   


if not(exist('fname', 'var'))
  fname = '';
end
if not(exist('wpos', 'var'))
  wpos = '';
end
if not(exist('ejes', 'var'))
  ejes = {[],[],[]};
end
if not(exist('ticks', 'var'))
  ticks = {[],[],[]};
end

ejes = struct('pop', ejes{1}, 'specs', ejes{2}, 'bio', ejes{3});
ticks = struct('pop', ticks{1}, 'specs', ticks{2}, 'bio', ticks{3});

fprintf('population size == number of individuals in each generation\n');
fprintf('biomass         == sum of amounts of branches of individuals in each generation\n');

location = 'Best';%'SouthEast';

figuras.tipos = figuras.tipos(toSelect);
figuras.tabla = figuras.tabla(:,toSelect);

nmi = 1;
coli = 2;
idxi = 3;
legendi = 4;
prmGi = 5;
prmFi = 6;
geni = 7;
popi = 8;
bioi = 9;
genci = 10;
if isnumeric(clustType)
  if numel(clustType)==1
    clustType = repmat(clustType, numel(figuras.tipos), 1);
  end
  divis = clustType;
else
  if ischar(clustType)
    clustType = repmat({clustType}, numel(figuras.tipos), 1);
  end
  divis = zeros(size(clustType));
  for z=1:numel(divis)
    switch clustType{z}
      case 'morpho'; divis(z) = 11;
      case 'ed';     divis(z) = 12;
      case 'edred';  divis(z) = 13;
      case 'mut';    divis(z) = 14;
    end
  end
end

if nargin>2
  if islogical(useDefaultColors) && useDefaultColors
    colors = {repmat(0.8, 1, 3), repmat(0.5, 1, 3), repmat(0, 1, 3)};
    figuras.tabla(coli,:) = colors;
  else
    if not(islogical(useDefaultColors))
      figuras.tabla(coli,:) = useDefaultColors;
    end
  end
end

hs=[];

% shows(1) = true;
% shows(2)=false;
% shows(3)=true;
% shows(4)=true;
% shows(5)=true;

if shows(1)
  %figura 7 (diversity vs population size)
  hs(end+1,1)=figure;
  for k=1:numel(figuras.tipos);
    kk=k;
    divi = divis(kk);
    zi = ismember(figuras.tabla{geni,kk}, figuras.tabla{genci,kk});
    line(figuras.tabla{divi, kk}, figuras.tabla{popi,kk}(zi), ...
         'LineWidth', 2.5, 'Color', figuras.tabla{coli,kk});
  end;
  xl = xlabel('diversity (number of clusters)'); yl = ylabel('population size'); %zl = zlabel('generations');
  legend(figuras.tabla(legendi, :), 'Location', location); 
  doForAll(xl, yl);%, zl);
  if not(isempty(ejes.specs))
    set(gca, 'XLim', ejes.specs);
  end
  if not(isempty(ejes.pop))
    set(gca, 'YLim', ejes.pop);
  end
  if not(isempty(ticks.specs))
    set(gca, 'XTick', ticks.specs);
  end
  if not(isempty(ticks.pop))
    set(gca, 'YTick', ticks.pop);
  end
  dowrite(fname, wpos, 'DP');

end


if shows(2)
  %diversity vs biomass
  hs(end+1,1)=figure;
  for k=1:numel(figuras.tipos); 
    kk=k; 
    divi = divis(z);
    zi = ismember(figuras.tabla{geni,kk}, figuras.tabla{genci,kk}); 
    line(figuras.tabla{divi, kk}, figuras.tabla{bioi,kk}(zi), ...
         'LineWidth', 2.5, 'Color', figuras.tabla{coli,kk}); 
  end; 
  xl = xlabel('diversity (number of species)'); yl = ylabel('biomass'); %zl = zlabel('generations');
  legend(figuras.tabla(legendi, :), 'Location', location); 
  doForAll(xl, yl);%, zl);
  if not(isempty(ejes.specs))
    set(gca, 'XLim', ejes.specs);
  end
  if not(isempty(ejes.bio))
    set(gca, 'YLim', ejes.bio);
  end
  if not(isempty(ticks.specs))
    set(gca, 'XTick', ticks.specs);
  end
  if not(isempty(ticks.bio))
    set(gca, 'YTick', ticks.bio);
  end
  %axis([0 120 0 36000]);
  dowrite(fname, wpos, 'DB');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    idxs = [genci divis(1); geni popi; geni bioi];
    ylabs = {'diversity', 'population size', 'biomass'};
  for mm=1:numel(ylabs)
    if not(shows(mm+2))
      continue;
    end
    %figura 7 (diversity vs population size)
    hs(end+1,1)=figure;
    for k=1:numel(figuras.tipos);
      kk=k;
  %     divi = divis(kk);
      line(figuras.tabla{idxs(mm,1),kk}, figuras.tabla{idxs(mm,2),kk}, ...
           'LineWidth', 2.5, 'Color', figuras.tabla{coli,kk});
    end;
    xl = xlabel('generations'); yl = ylabel(ylabs{mm}); %zl = zlabel('generations');
    legend(figuras.tabla(legendi, :), 'Location', location); 
    doForAll(xl, yl);%, zl);
    % if not(isempty(ejes.specs))
    %   set(gca, 'XLim', ejes.specs);
    % end
    % if not(isempty(ejes.pop))
    %   set(gca, 'YLim', ejes.pop);
    % end
    % if not(isempty(ticks.specs))
    %   set(gca, 'XTick', ticks.specs);
    % end
    % if not(isempty(ticks.pop))
    %   set(gca, 'YTick', ticks.pop);
    % end
    % dowrite(fname, wpos, 'DP');
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function doForAll(varargin)
fs = 14;
set(gca, 'FontSize', fs, 'fontname', 'times');
%grid on; %set(gca, 'XGrid', 'on'); 
%set(gca, 'XGrid', 'on'); 
for k=1:numel(varargin)
  set(varargin{k},  'FontSize', fs);
end


function dowrite(fname, wpos, typ)
if not(isempty(fname))
  name  = sprintf(fname, typ);
  if not(isempty(wpos))
    set(gcf, 'Visible', 'off');
    set(gcf, 'PaperUnit', 'points', 'PaperPosition', [0 0 wpos]);
  end
  saveas(gcf, name, 'png');
  if not(isempty(wpos))
    set(gcf, 'Visible', 'on');
  end
  close(gcf);
end
