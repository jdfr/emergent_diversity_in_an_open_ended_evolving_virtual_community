function [individualsSpecies proto offsetproto sim numGroups] = numberGroup(basedir,Generacion,Z,matrixDist,indicesClasificados)
%,basedir,Generacion)
% Z = linkage(dist);
load(sprintf('%s/P%03d.mat',basedir,Generacion));
P=eval(sprintf('P%03d', Generacion));
raster = P.raster;
offset = P.offset;

if ~isempty(Z)
    sZ = size(Z,1);
    c = 1;
    simi = 1;
    if ~isempty(indicesClasificados)
        raster= raster(indicesClasificados,:);
        offset = offset(indicesClasificados,:);
    end
    while c<=sZ+1 && simi >0

        % groups = cluster(Z,'cutoff',c);
        groups = cluster(Z,'maxclust',c,'criterion','distance');
        for i = 1 : c
            eachgroup{c}{i}= find(groups == i);
            [proto{c}{i} offsetproto{c}{i} meanDisrup{c}{i}]= prototype2(raster,offset,eachgroup{c}{i},matrixDist);
            individualsSpecies{c}{i} = indicesClasificados(eachgroup{c}{i});
        end


        sim(c) = sum([meanDisrup{c}{:}])/c;
        simi = sim(c);

        c = c + 1;
    end
else
    individualsSpecies{1}{1} = indicesClasificados;
    eachgroup{1}{1} = indicesClasificados;
    proto{1}{1} = raster{1};
    offsetproto{1}{1} = offset(1,:);
    sim = 0;
end

if sim == 0
    derivada = 0;
else
    derivada = abs(diff(sim));
end

h1 = figure('visible','off');
bar(sim);
xlabel('number of groups');
ylabel('dispersion');
saveas(h1,[basedir filesep 'dispersion' num2str(Generacion, '%03g')],'png');
close(h1);

h2 = figure('visible','off');
bar(derivada);
xlabel('number of groups');
ylabel('derivative of the dispersion');
saveas(h2,[basedir filesep 'derivadaDispersion' num2str(Generacion, '%03g')],'png');
close(h2);

f = find(derivada <= 1);
if numel(f) == 0
    f = length(derivada)+1;
end
numGroups = f(1);

%paintPrototypes(eachgroup{f(1)},proto{f(1)},offsetproto{f(1)},indicesClasificados,basedir,Generacion)
 


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function total = prototype(raster,offset,eachgroup)

% continue prototype
leachgroup=length(eachgroup);

for i = 1 : leachgroup

    raster1 = raster{eachgroup(i)};
    pixel1 = double([raster1{:}]);
    offset1 = offset(eachgroup(i),:);
    pixel1(:,1) = pixel1(:,1) - offset1(1);
    pixel1(:,2) = pixel1(:,2) - offset1(2);
    pixel1 = [pixel1 ones(1,size(pixel1,1))'];
    if i == 1
        total = pixel1;
    else
        [diftotalsvspixel1 pos1] = setdiff(total(:,1:2),pixel1(:,1:2),'rows');
        [difpixel1vstotal pos2] = setdiff(pixel1(:,1:2),total(:,1:2),'rows');

        totalrowsarepixel1 = setdiff(1:size(total,1),pos1);
        total(totalrowsarepixel1,3) = total(totalrowsarepixel1,3) + 1;
        total = [total ; pixel1(pos2,:)];



    end

end

if ~isempty(total)
    total(:,3) = total(:,3)/max(total(:,3));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [proto offsetproto mini]= prototype2(raster,offset,eachgroup,matrixDist)

% prototype = less disruptive individual of the group
leachgroup=length(eachgroup);

for i = 1 : leachgroup
    disruptive(i) = sum(matrixDist(eachgroup(i),eachgroup));
end

[mini pos]= min(disruptive);
minidisruptive = eachgroup(pos);
proto = raster{minidisruptive};
offsetproto = offset(minidisruptive,:);

% meanDisrup = sum(sum(matrixDist(eachgroup,eachgroup)))/leachgroup^2;


