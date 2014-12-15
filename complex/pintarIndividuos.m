function sizePNG=pintarIndividuos(basedir,finaldir,Generation,ind1,ind2)

poph = loadSimulation(basedir);

numIndividuals = zeros(1,Generation+1);
for j = 1 : Generation+1
    numIndividuals(j) = sum(poph.generation == j-1); % numIndividuals(1) = numero de individuos en la generacion 0
end

genomes = poph.genome(sum(numIndividuals(1:Generation))+1:sum(numIndividuals(1:Generation+1)));

for i = ind1:ind2
    fig = figure('visible','off');
    paintTree(genomes{i});
    saveas(fig,[finaldir,'ind' num2str(i, '%03g')],'png');
    close(fig);
end

tama = dir(finaldir);
for i = 3 : ind2-ind1+1+2
    sizePNG(i-2) = tama(i).bytes;
end