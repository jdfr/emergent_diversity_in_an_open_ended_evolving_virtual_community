function [individualsSpecies proto offsetproto sim numGroups] = numberGroup2sim(basedir1,Gen1,basedir2,Gen2,Z,matrixDist,indicesClasificados1,indicesClasificados2,finaldir)
%,basedir,Generacion)
% Z = linkage(dist);
load(sprintf('%s/P%03d.mat',basedir1,Gen1));
P1=eval(sprintf('P%03d', Gen1));
raster1 = P1.raster;
offset1 = P1.offset;

load(sprintf('%s/P%03d.mat',basedir2,Gen2));
P2=eval(sprintf('P%03d', Gen2));
raster2 = P2.raster;
offset2 = P2.offset;


if ~isempty(Z)
    sZ = size(Z,1);
    c = 1;
    simi = 1;
    if ~isempty(indicesClasificados1)
        raster1= raster1(indicesClasificados1,:);
        offset1 = offset1(indicesClasificados1,:);
    end
    if ~isempty(indicesClasificados2)
        raster2= raster2(indicesClasificados2,:);
        offset2 = offset2(indicesClasificados2,:);
    end
    lraster1 = size(raster1,1);

    raster = [raster1;raster2];
    offset = [offset1; offset2];

    while c<=sZ+1 && simi >0

        % groups = cluster(Z,'cutoff',c);
        groups = cluster(Z,'maxclust',c,'criterion','distance');
        for i = 1 : c
            eachgroup{c}{i}= find(groups == i);
            [proto{c}{i} offsetproto{c}{i} meanDisrup{c}{i}]= prototype2(raster,offset,eachgroup{c}{i},matrixDist);
            for j = 1 : length(eachgroup{c}{i})
                if eachgroup{c}{i}(j) <= lraster1
                    individualsSpecies{c}{i}(j) = indicesClasificados1(eachgroup{c}{i}(j));
                else
                    individualsSpecies{c}{i}(j) = indicesClasificados2(eachgroup{c}{i}(j)-lraster1);
                end
            end
        end
        sim(c) = sum([meanDisrup{c}{:}])/c;
        simi = sim(c);

        c = c + 1;
    end
else
    individualsSpecies{1}{1} = [indicesClasificados1 indicesClasificados2];
    eachgroup{1}{1} = [indicesClasificados1 indicesClasificados2];
    proto{1}{1} = [raster1{1};raster2{1}];
    offsetproto{1}{1} = [offset1(1,:);offset2(1,:)];
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
saveas(h1,[finaldir filesep 'dispersion'],'png');
close(h1);

h2 = figure('visible','off');
bar(derivada);
xlabel('number of groups');
ylabel('derivative of the dispersion');
saveas(h2,[finaldir filesep 'derivadaDispersion'],'png');
close(h2);

f = find(derivada <= 1);
if numel(f) == 0
    f = length(derivada)+1;
end
numGroups = f(1);

paintPrototypes2sim(individualsSpecies{f(1)},proto{f(1)},offsetproto{f(1)},finaldir)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function paintPrototypes2sim(individualsSpecies,proto,offsetproto,finaldir)

grandote = true;

leachgroup = length(individualsSpecies);
drawing.FontSize = 6;

m = ceil(sqrt(leachgroup));
if(m * (m-1) >= leachgroup)
    n = m-1;
else
    n = m;
end

dimension = zeros(leachgroup, 2);
notEmptyTrees = cellfun(@(x) (~isempty(x{1})), proto);

dimension(notEmptyTrees, :) = cell2mat(cellfun(@(x) ([max(x{1}) max(x{2})]), proto, 'UniformOutput', false)');

maxDim = max(dimension);

h = figure('visible','off');

for i = 1 : leachgroup

    if(grandote)
        subAxis = [0 dimension(i,2) 0 dimension(i,1)];
    else
        subAxis = [[0 maxDim(2)]-((maxDim(2)-dimension(i,2))/2) 0 maxDim(1)];
    end
    subplot(m,n,i);%,'align');
    drawtree(proto{i},dimension(i,:),offsetproto{i},0,1,[0 0 0], subAxis);
    title({mat2str(individualsSpecies{i}),['DIMENSION: ',mat2str(dimension(i,:))]},'FontSize',drawing.FontSize);

end
drawnow;
saveas(h,[finaldir filesep 'grupos'],'png');
close(h);



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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sim = similitud2(raster1,offset1,raster2,offset2)%, show)

if isempty([raster1{:}]) && isempty([raster2{:}])
    sim = 1;
elseif (isempty([raster1{:}]) && ~isempty([raster2{:}])) || (~isempty([raster1{:}]) && isempty([raster2{:}]))
    sim = 0;
else
    pixel1 = [raster1{:}];
    pixel2 = [raster2{:}];

    pixel1(:,1) = pixel1(:,1) - offset1(1);
    pixel1(:,2) = pixel1(:,2) - offset1(2);

    pixel2(:,1) = pixel2(:,1) - offset2(1);
    pixel2(:,2) = pixel2(:,2) - offset2(2);

    %     sraster1 = size(pixel1,1);
    %     sraster2 = size(pixel2,1);

    sim = size(intersect(pixel1,pixel2,'rows'),1)/size(union(pixel1,pixel2,'rows'),1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

function showInf(src,event,raster,dimension,offset,H,perm,indicesClasificados)
%This function shows the individual's information in the figure
ax = get(H(1),'Parent');

%get current point

%this means that the function has been called by clicking on the tree
currentPoint= get(ax, 'CurrentPoint');
plotTree    = true;
cx            = currentPoint(1,1);
cy            = currentPoint(1,2);
cx = round(cx);
fprintf('cx=%g, cy=%g\n', cx, cy);
if plotTree
    if isempty(indicesClasificados)
        if isempty(raster{perm(cx)}{1})
            sprintf('empty individual')
        else
            fig = figure;
            drawtree(raster{perm(cx)}, dimension(perm(cx),:), offset(perm(cx),:), 0, 1, [0 0 0], 'fit');
            title(sprintf('individual: %03d', perm(cx)));
        end
    else
        if isempty(raster{perm(cx)}{1})
            sprintf('empty individual')
        else
            fig = figure;
            drawtree(raster{perm(cx)}, dimension(perm(cx),:), offset(perm(cx),:), 0, 1, [0 0 0], 'fit');
            title(sprintf('individual: %03d', indicesClasificados(perm(cx))));
        end
    end
end
