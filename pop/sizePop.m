function numIndividuals = sizePop(basedir,numGenerations)

poph = loadSimulation(basedir);

numIndividuals = zeros(1,numGenerations+1);
for j = 1 : numGenerations+1
    numIndividuals(j) = sum(poph.generation == j-1); % numIndividuals(1) = numero de individuos en la generacion 0
end
