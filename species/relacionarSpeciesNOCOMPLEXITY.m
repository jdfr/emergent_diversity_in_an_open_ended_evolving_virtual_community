function [especies ]= relacionarSpeciesNOCOMPLEXITY(arg1,numGenerations)

if ischar(arg1)
  basedir = arg1;
  poph = loadSimulation(basedir);
elseif iscell(arg1)
  poph = arg1{1};
  basedir = arg1{2};
else
  error('arg1 has incorrect type!!!!');
end

numIndividuals = zeros(1,numGenerations+1);
for i = 1 : numGenerations+1
    numIndividuals(i) = sum(poph.generation == i-1); % numIndividuals(1) = numero de individuos en la generacion 0
end

indicesPredecesores = poph.tree.iA;

for i = 1 : numGenerations-1
    i
    predecesores = indicesPredecesores(sum(numIndividuals(1:i+1))+1:sum(numIndividuals(1:i+2)));

    load(sprintf('%s/numGroups%03d.mat',basedir,i));
    numGroupsA = eval(sprintf('numGroups%03d', i));

    load(sprintf('%s/individualsSpecies%03d.mat',basedir,i));
    individualsSpeciesA = eval(sprintf('individualsSpecies%03d', i));

    load(sprintf('%s/numGroups%03d.mat',basedir,i+1));
    numGroupsS = eval(sprintf('numGroups%03d', i+1));

    load(sprintf('%s/individualsSpecies%03d.mat',basedir,i+1));
    individualsSpeciesS = eval(sprintf('individualsSpecies%03d', i+1));

    for j = 1 : numGroupsA
        longEspecie = length(individualsSpeciesA{numGroupsA}{j});
        cont = cell(longEspecie,1);
        cont2 = cell(longEspecie,1);
        for k = 1 : longEspecie
            t = 1;
            f = find(predecesores == individualsSpeciesA{numGroupsA}{j}(k));
            for q = 1 : length(f)
                for p = 1 : numGroupsS
                    if ~isempty(find(individualsSpeciesS{numGroupsS}{p}==f(q), 1))
                        cont{k}{t} = p;
                        t = t + 1;
                    end
                end
            end
            if ~isempty(cont{k})
                cont2{k}=unique([cont{k}{:}]);
            else
                cont2{k}=[];
            end
        end
        especies{i}{j}=unique([cont2{:}]);
    end
end

tosave1 = ['especies' num2str(numGenerations, '%03g')];

eval([tosave1 ' = especies;']);

save([basedir filesep tosave1],tosave1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dist=distancia(X)

long = length(X);
k = 1;
for i = 1 : long-1
    for j = (i+1): long
        dist(k) = abs(X(j) - X(i));
        k = k + 1;
    end
end

