function antecesores(basedir, generation, rangeid, individual)

g = reapIndividualGenealogy(basedir, generation, rangeid, individual);

ancestors = g.ancestor;
generation = g.gDescendant;
cambios = g.change;
sa = size(ancestors,1)-1;
drawing.rowcols = [ceil(sqrt(sa)) ceil(sqrt(sa))];
fig = figure; %('Visible', 'off', 'Units', 'pixels', 'PaperUnits', 'points');
for i = 1 : sa
    ramas = ls2(ancestors{i}, 'canonical');
    subplot(drawing.rowcols(1),drawing.rowcols(2),i); %,'align'
        title(sprintf('%s\n%g=>%s\n',ancestors{i},generation(i),cambios{i}));%4);

%     plot([1 3],[2 5]);
    drawtree(ramas, 0,'k');
%     paintTree(ancestors{i});
%     axis equal
    axis([0 1000 0 1000])

    axis off
end
drawnow
% saveas(fig,[nomdir,'sample' num2str(Generation, '%03g')],'png');
% close(fig);



   