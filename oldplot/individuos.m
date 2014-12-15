function individuos(P,finaldir)

long = length(P.raster);
dimension = P.dimensions;
maxDim = max(dimension);

for i = 1 : long
    h = figure('visible','off');
    subAxis = [[0 maxDim(2)]-((maxDim(2)-dimension(i,2))/2) 0 maxDim(1)];
    drawtree(P.raster{i},dimension(i,:),P.offset(i,:),0,1,[0 0 0], subAxis);
    saveas(h,[finaldir,'individuo' num2str(i, '%03g')],'png');
    close(h);
end

