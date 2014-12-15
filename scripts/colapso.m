function cola = colapso(basedir)

for alphaG = 0.1:0.1:1
    alphaG
    for alphaF = 0.1:0.1:1
        alphaF
        direc1 = [basedir filesep 'G' mat2str(alphaG) 'F' mat2str(alphaF)];
        direc2 = dir(direc1);
        [a b direc3] = direc2.name;
        col = zeros(1,5);
        for i = 1 : 5
            direc = [basedir filesep 'G' mat2str(alphaG) 'F' mat2str(alphaF) filesep direc3 filesep '1=0.0001_' mat2str(i)];
            poph = loadSimulation(direc);
            numGenerations = max(poph.generation);

            numIndividuals = zeros(1,numGenerations+1);
            for j = 1 : numGenerations+1
                numIndividuals(j) = sum(poph.generation == j-1); % numIndividuals(1) = numero de individuos en la generacion 0
            end
            maxi = max(numIndividuals);
            if numIndividuals(numGenerations+1) < maxi-((25*maxi)/100) || numGenerations < 500
                col(i) = 1;
            end
        end
        cola(round(alphaG*10),round(alphaF*10)) = sum(col)

    end
end

