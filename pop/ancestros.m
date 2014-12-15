function ances = ancestros(basedir,indice,generacion)

poph = loadSimulation(basedir);

numGenerations = max(poph.generation);

numIndividuals = zeros(1,numGenerations);
for i = 1 : numGenerations+1
    numIndividuals(i) = sum(poph.generation == i-1); % numIndividuals(1) = numero de individuos en la generacion 0
end

indicesPredecesores = poph.tree.iA;

predecesores = cell(1,numGenerations);
for i = 2 : numGenerations+1  
    if i == 668
        disp('hola')
    end
    predecesores{i-1} = indicesPredecesores(sum(numIndividuals(1:i-1))+1:sum(numIndividuals(1:i)));
end

ances = ones(1,generacion);
for i = generacion:-1:1
    if i == 668
        disp('hola')
    end
    ances(i) = predecesores{i}(indice);
    indice = ances(i);
end