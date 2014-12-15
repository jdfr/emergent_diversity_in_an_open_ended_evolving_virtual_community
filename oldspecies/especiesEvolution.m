function especiesEvolution(basedir)

poph = loadSimulation(basedir);
numGenerations = max(poph.generation);

for i = 1 : numGenerations
    load(sprintf('%s\\P%03d.mat',basedir,i));
    [Z matrixDist indicesClasificados] = HierarchicalClustering(eval(sprintf('P%03d', i)),basedir,i);
    [eachgroup proto offsetproto sim numGroups] = numberGroup(eval(sprintf('P%03d', i)),Z,matrixDist,indicesClasificados,basedir,i);
    
    tosave1 = ['Z' num2str(i, '%03g')];
    eval([tosave1 ' = Z;']);
    save([basedir filesep tosave1],tosave1);
    
    tosave2 = ['matrixDist' num2str(i, '%03g')];
    eval([tosave2 ' = matrixDist;']);
    save([basedir filesep tosave2],tosave2);
    
    tosave3 = ['indicesClasificados' num2str(i, '%03g')];
    eval([tosave3 ' = indicesClasificados;']);
    save([basedir filesep tosave3],tosave3);
    
    tosave4 = ['eachgroup' num2str(i, '%03g')];
    eval([tosave4 ' = eachgroup;']);
    save([basedir filesep tosave4],tosave4);
    
    tosave5 = ['proto' num2str(i, '%03g')];
    eval([tosave5 ' = proto;']);
    save([basedir filesep tosave5],tosave5);
    
    tosave6 = ['sim' num2str(i, '%03g')];
    eval([tosave6 ' = sim;']);
    save([basedir filesep tosave6],tosave6);
    
    tosave7 = ['offsetproto' num2str(i, '%03g')];
    eval([tosave7 ' = offsetproto;']);
    save([basedir filesep tosave7],tosave7);
    
    tosave8 = ['numGroups' num2str(i, '%03g')];
    eval([tosave8 ' = numGroups;']);
    save([basedir filesep tosave8],tosave8);
end


