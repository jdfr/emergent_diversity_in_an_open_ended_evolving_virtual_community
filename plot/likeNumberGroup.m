function likeNumberGroup(basedir)
c=95;
gen=450;

nm = sprintf('Z%03d', gen);
z=load([basedir filesep nm '.mat']);
Z=z.(nm);
nm = sprintf('indicesClasificados%03d', gen);
z=load([basedir filesep nm '.mat']);
indicesClasificados=z.(nm);
nm = sprintf('individualsSpecies%03d', gen);
z=load([basedir filesep nm '.mat']);
individualsSpecies2=z.(nm);
nm = sprintf('matrixDist%03d', gen);
z=load([basedir filesep nm '.mat']);
matrixDist=z.(nm);

eachgroup=repmat({{}}, c, 1);
meanDisrup=repmat({{}}, c, 1);
individualsSpecies=repmat({{}}, c, 1);
for c=1:c
groups = cluster(Z,'maxclust',c,'criterion','distance');
        for i = 1 : c
            eachgroup{c}{i}= find(groups == i);
            [meanDisrup{c}{i}]= prototype2(eachgroup{c}{i},matrixDist);
            individualsSpecies{c}{i} = indicesClasificados(eachgroup{c}{i});
        end


        sim(c) = sum([meanDisrup{c}{:}])/c;
        simi = sim(c);
end
        
        z=z;
function [mini]= prototype2(eachgroup,matrixDist)

% prototype = less disruptive individual of the group
leachgroup=length(eachgroup);

for i = 1 : leachgroup
    disruptive(i) = sum(matrixDist(eachgroup(i),eachgroup));
end

[mini pos]= min(disruptive);
