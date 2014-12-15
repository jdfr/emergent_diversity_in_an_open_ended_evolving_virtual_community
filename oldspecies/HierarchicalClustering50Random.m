function [randpos,randaxioms,dist]=HierarchicalClustering50Random(basedir,Generation,numIndividuals)


% las generaciones empiezan en 0, tal y como aparecen en poblacion.txt

pop = loadSimulation(basedir);
[axioms{1:numIndividuals}] = pop.genome{((Generation)*numIndividuals)+1 : ((Generation)*numIndividuals)+numIndividuals};

r = randperm(200);
randpos = r(1:20);

[randaxioms{1:20}] = axioms{randpos};

dist = zeros(1,(20*(20-1))/2);
k = 1;
for i = 1 : 20-1
    i
    axiom2 = randaxioms{i};
    for j = i+1 : 20
        j
        axiom1 = randaxioms{j};
        dist(k) = 1 - similitud(axiom1,axiom2);
        k = k + 1;
    end
end

Z = linkage(dist);
dendrogram(Z,0);

% [...] = dendrogram(...,'colorthreshold',t) 
% assigns a unique color to each group of nodes in the dendrogram where the linkage is less than the threshold t. 
% t is a value in the interval [0,max(Z(:,3))]. 
% Setting t to the string 'default' is the same as t = .7(max(Z(:,3))). 
% 0 is the same as not specifying 'colorthreshold'. 
% The value max(Z(:,3)) treats the entire tree as one group and colors it all one color.










