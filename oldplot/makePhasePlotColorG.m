function makePhasePlotColorG(x,y,generations) %#ok<DEFNU>

% figargs = {'Visible', 'off', 'Units', 'pixels', 'PaperUnits', 'points', 'PaperPosition', [0 0 paperSize paperSize]};
lineargs   = {'LineStyle', 'none', 'Marker', 'd'};%, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'};
cmap = 'cool';

% fig = figure(figargs{:});
cols = generations-min(generations)+1;
cmap3=feval(cmap, (max(cols)));
set(gcf, 'Colormap', cmap3);
hc = colorbar('Location', 'EastOutside');%, 'YTick', numgens, 'YTickLabel', strgens);
freezeColors(hc);

line(x,y, 'Color', 'k');
%for k=1:numel(x)
for k=numel(x):-1:1
  color = cmap3(cols(k),:);
  line(x(k),y(k), lineargs{:}, 'MarkerFaceColor', color, 'MarkerEdgeColor', color);
end
grid on; 
% xlabel(xlab); ylabel(ylab); title(titulo);
% if tofile
%   drawnow;
%   saveas(fig,[basedir filesep 'species_' filename '.png'],'png');
%   close(fig);
% end