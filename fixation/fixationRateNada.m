function [fixationRate naturaleza] = fixationRateNada(basedir)

poph = loadSimulation(basedir);

numGenerations = max(poph.generation);

numIndividuals = zeros(1,numGenerations);
for i = 1 : numGenerations+1
    numIndividuals(i) = sum(poph.generation == i-1); % numIndividuals(1) = numero de individuos en la generacion 0
end


changes = poph.tree.change;
indicesPredecesores = poph.tree.iA;

changesG{1}= changes(1:numIndividuals(1));
indicesPredecesoresG{1}= indicesPredecesores(1:numIndividuals(1));

for i = 2 : numGenerations+1
    changesG{i} = changes(sum(numIndividuals(1:i-1))+1:sum(numIndividuals(1:i)));
    indicesPredecesoresG{i} = indicesPredecesores(sum(numIndividuals(1:i-1))+1:sum(numIndividuals(1:i)));
end
    
for i = 1 : numGenerations
    lchanges = length(changesG{i});
    duplNada = [];
    indDupl = [];
    for j = 1 : lchanges
        duplNada = find(changesG{i}{j} == '=');
        if ~isempty(duplNada)
            indDupl = [indDupl j];
        end
    end
    if isempty(indDupl)
        fixationRate(i) = 0;
    else
        lindDupl = length(indDupl);
        fij = 0;
        for j = 1 : lindDupl
            f = find(indicesPredecesoresG{i+1} == indDupl(j));
            if ~isempty(f)
                fij = fij + 1;
            end
        end
        fixationRate(i) = fij/lindDupl; %fijadas/ocurridas
    end
    naturaleza(i) = 1/(2*numIndividuals(i));
end
        
        
        
        
        