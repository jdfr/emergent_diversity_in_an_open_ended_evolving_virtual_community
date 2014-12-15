function allIndividualsOfAnSpecie(basedir,Generacion,individualsSpecies)

numSpecies = length(individualsSpecies);
load(sprintf('%s/P%03d.mat',basedir,Generacion));
P=eval(sprintf('P%03d', Generacion));

for i = 1 : numSpecies
    i
    h = figure('visible','off');
    les = length(individualsSpecies{i});
    m = ceil(sqrt(les));
    if(m * (m-1) >= les)
        n = m-1;
    else
        n = m;
    end
    for j = 1 : les
        subplot(m,n,j)
        dimension = [max(P.raster{individualsSpecies{i}(j)}{1}) max(P.raster{individualsSpecies{i}(j)}{2})];
                drawtree(P.raster{individualsSpecies{i}(j)},dimension,P.offset(individualsSpecies{i}(j),:),0,1,[0 0 0], 'fit');
        title({['INDIVIDUAL: ',mat2str(individualsSpecies{i}(j))],['DIMENSION: ',mat2str(dimension)]},'FontSize',6);
    end
            saveas(h,[basedir,'especie' num2str(i,'%03g')],'png');
        close(h);

end
