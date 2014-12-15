function h = showCompression(poph, compression, titname, colmap, useLogY)

puthalfticks=true;

minc = min(compression);
maxc = max(compression);
ming = min(poph.generation);
maxg = max(poph.generation);
if useLogY
  minc = log10(minc);
  maxc = log10(maxc);
  compression = log10(compression);
  cs = linspace(minc, maxc, 100);
else
  cs = linspace(minc, maxc, 100);
end
gs = ming:maxg;
generations = poph.generation;
h = myhist3(generations, 'generations', gs, [], compression, 'complexity', cs, realpow(10, cs), true, colmap);
if useLogY
  a = get(h, 'currentaxes');
  %yt = get(a, 'YTick');
  yt = ceil(minc):floor(maxc);
  if puthalfticks
    yt = [log10(realpow(10, yt')/2), yt']';
    yt = yt(:);
    yt = yt(2:end); %skip the first half
    ytceil = ceil(maxc);
    ytceil2 = log10(realpow(10, ytceil)/2);
    if ytceil2<=maxc
      yt(end+1) = ytceil2;
    end
  end
  set(a, 'YTick', yt, 'YTickLabel', realpow(10, yt));
end
if not(isempty(titname))
  title(titname);
end
% set(a, 'YScale', 'log');

% [H3 C3] = hist3([generations, compression], {gs, cs});
% 
% 
% 
% complexities = linspace(
% 
% axesargs = {'FontName','times','FontSize',18};
% lw = 2;
% lw2 = 4;
% linefun = @stairs;%@line; %@stairs;
% 
% if not(exist('filename', 'var'))
%   filename = '';
% end
% 
% figargs = {'Color', 'w'};
% 
% if not(isempty(filename))
%   figargs = [figargs {'visible', 'off'}];
% end
% 
% h = figure(figargs{:});
%   
% x=axes;
% set(x,axesargs{:}, 'YLim', [0 300]);
% hold on
% cols = {repmat(0, 1, 3), repmat(0.4, 1, 3), repmat(0.7, 1,3)};
% for k=1:numel(res)
%   linefun(res{k}.generations, res{k}.nominimized.mean,  'Color', cols{k}, 'LineWidth', lw);
%   nm{end+1,1} = nm1{k};
% end
% for k=1:numel(res)
%   linefun(res{k}.generations, res{k}.minimized.mean,  'Color', cols{k}, 'LineWidth', lw2, 'LineStyle', '-');%':');
%   nm{end+1,1} = nm2{k};
% end
% lh = legend(nm,'Location','NorthWest');
% set(lh, legendargs{:});
% xlabel('generations')
% ylabel('genome length')
% 
% if not(isempty(filename))
%   img = imcapture(h, 'all', 600);
%   imwrite(img, filename);
%   close(h);
%   clear img;
% end
%   
% nm1 = cellfun(@(x)['mean ratio between genome lengths (' x ')'],         names, 'uniformoutput', false);
% nm2 = cellfun(@(x)[x ' (reduced)'], names, 'uniformoutput', false);
% nm3 = cellfun(@(x)[x ' (nonreduced)'], names, 'uniformoutput', false);
% nm =  {};
% figure('Color', 'w');
% x=axes;
% set(x,axesargs{:});
% hold on
% % for k=1:numel(res)
% %   linefun(res{k}.generations, res{k}.compression.mean,  'Color', cols{k}, 'LineWidth', lw);
% %   nm{end+1,1} = nm1{k};
% % end
% for k=1:numel(res)
%   linefun(res{k}.generations, res{k}.nonminT3,  'Color', cols{k}, 'LineWidth', lw, 'LineStyle', '-');
%   nm{end+1,1} = nm3{k};
% end
% for k=1:numel(res)
%   linefun(res{k}.generations, res{k}.minT3,  'Color', cols{k}, 'LineWidth', lw2, 'LineStyle', '-');%':');
%   nm{end+1,1} = nm2{k};
% end
% lh = legend(nm,'Location','NorthWest');
% set(lh, legendargs{:});
% xlabel('generations')
% ylabel('ratio')
%   
