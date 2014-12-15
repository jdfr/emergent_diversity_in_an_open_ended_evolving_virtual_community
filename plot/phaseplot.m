function phaseplot(figuras, toSelect, useDefaultColors, ejes, ticks, fname, wpos)
%ejes: {[minPop maxPop] [minSpecs maxSpecs] [minBio maxBio]}  
%ticks: lo mesmo con ticks

%load(figurass_6sims)
%figurass(figuras, [6, 2, 5], true, {[0 800], [0 150], [0 3.6e4]}, {0:100:800, 0:20:120, 0:0.5e4:3.5e4}, 'PHASEPLOT_%s', [1000 750]);   

if nargin<6
  fname = '';
  wpos=[];
end
if nargin<7
  wpos=[];
end
% if (nargin<1) || isempty(figuras)
%   if nargin<2
%     useall=false;
%   end
%   load figurass;
%   if ~useall
%     figuras = figuras3;
%   end
%   clear figuras3;
% end

if nargin<4
  ejes= {[],[],[]};
end
ejes = struct('pop', ejes{1}, 'specs', ejes{2}, 'bio', ejes{3});
if nargin<5
  ticks= {[],[],[]};
end
ticks = struct('pop', ticks{1}, 'specs', ticks{2}, 'bio', ticks{3});

fprintf('population size == number of individuals in each generation\n');
fprintf('biomass         == sum of amounts of branches of individuals in each generation\n');

%kk=find(cell2mat(figuras.tabla(9,:))==k);

location = 'SouthEast';

figuras.tipos = figuras.tipos(toSelect);
figuras.tabla = figuras.tabla(:,toSelect);
if nargin>2
  if islogical(useDefaultColors) && useDefaultColors
    colors = {repmat(0.8, 1, 3), repmat(0.5, 1, 3), repmat(0, 1, 3)};
    figuras.tabla(8,:) = colors;
  else
    figuras.tabla(8,:) = useDefaultColors;
  end
end

% %figura 6 (history vs population size)
% figure;
% for k=1:numel(figuras.tipos); kk=k; line(figuras.tabla{2, kk}, figuras.tabla{3,kk}, 'LineWidth', 2.5, 'Color', figuras.tabla{8,kk}); end; legend(figuras.tabla(7, :)); xl = xlabel('generations'); yl = ylabel('population size');
% doForAll(xl, yl);


%figura 7 (diversity vs population size)
figure;
for k=1:numel(figuras.tipos); kk=k; zi = ismember(figuras.tabla{2,kk}, figuras.tabla{6,kk}); line(figuras.tabla{5, kk}, figuras.tabla{3,kk}(zi), 'LineWidth', 2.5, 'Color', figuras.tabla{8,kk}); end;
xl = xlabel('diversity (number of species)'); yl = ylabel('population size'); %zl = zlabel('generations');
legend(figuras.tabla(7, :), 'Location', location); 
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
%axis([0 300 0 1200]);
%%axis([0 120 0 800]);
% figuras.tabla{2,4}=figuras.tabla{2,4}(1:501);
% figuras.tabla{3,4}=figuras.tabla{3,4}(1:501);
% figuras.tabla{4,4}=figuras.tabla{4,4}(1:501);
% figuras.tabla{5,4}=figuras.tabla{5,4}(1:500);
% figuras.tabla{6,4}=figuras.tabla{6,4}(1:500);
dowrite(fname, wpos, 'DP');



% %history vs biomass
% figure;
% for k=1:numel(figuras.tipos); kk=k; line(figuras.tabla{2, kk}, figuras.tabla{4,kk}, 'LineWidth', 2.5, 'Color', figuras.tabla{8,kk}); end; legend(figuras.tabla(7, :)); xl = xlabel('generations'); yl = ylabel('biomass');
% doForAll(xl, yl);


%diversity vs biomass
figure;
for k=1:numel(figuras.tipos); kk=k; zi = ismember(figuras.tabla{2,kk}, figuras.tabla{6,kk}); line(figuras.tabla{5, kk}, figuras.tabla{4,kk}(zi), 'LineWidth', 2.5, 'Color', figuras.tabla{8,kk}); end; 
xl = xlabel('diversity (number of species)'); yl = ylabel('biomass'); %zl = zlabel('generations');
legend(figuras.tabla(7, :), 'Location', location); 
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

% %diversity vs biomass
% figure;
% for k=1:numel(figuras.tipos); kk=k; zi = ismember(figuras.tabla{2,kk}, figuras.tabla{6,kk}); line(figuras.tabla{5, kk}, figuras.tabla{4,kk}(zi), figuras.tabla{3,kk}(zi), 'LineWidth', 2.5, 'Color', figuras.tabla{8,kk}); end; 
% xl = xlabel('diversity (number of species)'); yl = ylabel('biomass'); zl = zlabel('population size');
% legend(figuras.tabla(7, :), 'Location', location); 
% doForAll(xl, yl, zl);
% if not(isempty(ejes.specs))
%   set(gca, 'XLim', ejes.specs);
% end
% if not(isempty(ejes.bio))
%   set(gca, 'YLim', ejes.bio);
% end
% if not(isempty(ejes.pop))
%   set(gca, 'ZLim', ejes.pop);
% end
% if not(isempty(ticks.specs))
%   set(gca, 'XTick', ticks.specs);
% end
% if not(isempty(ticks.bio))
%   set(gca, 'YTick', ticks.bio);
% end
% if not(isempty(ticks.pop))
%   set(gca, 'ZTick', ticks.pop);
% end
% view(3);
% %axis([0 120 0 36000 0 800]);

% %biomass vs population size
% figure;
% for k=1:numel(figuras.tipos); kk=k; line(figuras.tabla{4,kk}, figuras.tabla{3, kk}, 'LineWidth', 2.5, 'Color', figuras.tabla{8,kk}); end; legend(figuras.tabla(7, :)); xl = xlabel('biomass'); yl = ylabel('population size');
% doForAll(xl, yl);

function doForAll(varargin)
fs = 28;
set(gca, 'FontSize', fs);
%grid on; %set(gca, 'XGrid', 'on'); 
set(gca, 'XGrid', 'on'); 
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
