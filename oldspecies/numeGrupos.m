function [pop grupos] = numeGrupos(basedir,GenI,GenF)

poph = loadSimulation(basedir);
numIndividuals = zeros(1,GenF+1);
for j = 1 : GenF+1
    numIndividuals(j) = sum(poph.generation == j-1); % numIndividuals(1) = numero de individuos en la generacion 0
end
for i = GenI : GenF
    load(sprintf('%s/numGroups%03d.mat',basedir,i));
    numGroups = eval(sprintf('numGroups%03d', i));

    grupos(i) = numGroups;
end
pop = numIndividuals(2:end);
plot(pop,grupos)

xlabel('size of the population')
ylabel('number os species')
