function [individualsSpecies sim numGroups] = numberGroupNOPROTO(basedir,Generacion,Z,matrixDist,indicesClasificados, options)
%,basedir,Generacion)
% Z = linkage(dist);

if ~isempty(Z)
    sZ = size(Z,1);
    c = 1;
    simi = 1;
    while c<=sZ+1 && simi >0

        % groups = cluster(Z,'cutoff',c);
        groups = cluster(Z,'maxclust',c,'criterion','distance');
        for i = 1 : c
            eachgroup{c}{i}= find(groups == i);
            [meanDisrup{c}{i}]= prototype2(eachgroup{c}{i},matrixDist,options);
            individualsSpecies{c}{i} = indicesClasificados(eachgroup{c}{i});
        end


        sim(c) = sum([meanDisrup{c}{:}])/c;
        simi = sim(c);

        c = c + 1;
    end
else
    individualsSpecies{1}{1} = indicesClasificados;
    eachgroup{1}{1} = indicesClasificados;
    sim = 0;
end

 numGroups = calculateNumGroups(sim, options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sim = similitud3(pixel1,proto)

% prototype = proto

if isempty(pixel1) && isempty(proto)
    sim = 1;
elseif (isempty(pixel1) && ~isempty(proto)) || (isempty(proto) && ~isempty(pixel1))
    sim = 0;
else
    [inter rows] = intersect(proto(:,1:2),pixel1,'rows');
    suma = sum(proto(rows,3));
    sumatotal = sum(proto(:,3));

    sim = suma/sumatotal;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mini]= prototype2(eachgroup,matrixDist,options)

% prototype = less disruptive individual of the group
leachgroup=length(eachgroup);

if options.oldWayToFindPrototype
  for i = 1 : leachgroup
      disruptive(i) = sum(matrixDist(eachgroup(i),eachgroup));
  end
else
  for i = 1 : leachgroup
      disruptive(i) = sum(matrixDist(eachgroup(i),eachgroup))/leachgroup;
  end
end

[mini pos]= min(disruptive);

% meanDisrup = sum(sum(matrixDist(eachgroup,eachgroup)))/leachgroup^2;


