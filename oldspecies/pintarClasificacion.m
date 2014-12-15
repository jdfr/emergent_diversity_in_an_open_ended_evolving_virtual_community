function fig = pintarClasificacion(basedir,numSpecies,individualsSpecies,protos,offsets,numGenerations)

independiente = false;

h1 = figure;
A = numSpecies;
for i = 0 : numGenerations - 1
    numSpecies{i+1}(find(numSpecies{i+1}==0))=[];
    numSpeciesperG(i+1) = length(find(numSpecies{i+1}>0));
    plot(i,numSpecies{i+1},'b*');
    hold on;
end
xlabel('generations')
ylabel('species')

h2 = figure;
plot(numSpeciesperG);
xlabel('generations')
ylabel('number of species')


esp = find(A{numGenerations}>0);
proto = protos{numGenerations}(esp);
offset = offsets{numGenerations}(esp);
eachgroup = individualsSpecies{numGenerations}(esp);

if independiente
    for i = 1 : length(proto)
        fig = figure('Visible', 'off', 'Units', 'pixels', 'PaperUnits', 'points', 'PaperPosition', [0 0 1200 1200]);
        dimension = [max(proto{i}{1}) max(proto{i}{2})];
        subAxis = [0 dimension(1,2) 0 dimension(1,1)];
        drawtree(proto{i},dimension,offset{i},0,1,[0 0 0], subAxis);
        title({mat2str(eachgroup{i}),['DIMENSION: ',mat2str(dimension)]},'FontSize',6);
        drawnow;
        saveas(fig,[basedir,'especie' num2str(i,'%03g')],'png');
        close(fig);
    end
else
fig = figure;
%     fig = figure('Visible', 'off', 'Units', 'pixels', 'PaperUnits', 'points', 'PaperPosition', [0 0 10000 10000]);
    paintPrototypes(eachgroup,proto,offset,[]);
end