function [individualsSpecies,protos,offsets]=clusteringMorfPlusFilo(basedir,poph,Generacion)

% poph = loadSimulation(basedir);
nPixelsInSoil = poph.nPixelsInSoil;
load(sprintf('%s/P%03d.mat',basedir,Generacion));
P=eval(sprintf('P%03d', Generacion));

numIndividuals = zeros(1,Generacion+1);
for i = 1 : Generacion+1
    numIndividuals(i) = sum(poph.generation == i-1); % numIndividuals(1) = numero de individuos en la generacion 0
end
nPixelsInSoil = poph.nPixelsInSoil;

ramasNegativas = nPixelsInSoil(sum(numIndividuals(1:Generacion))+1:sum(numIndividuals(1:Generacion+1)));

[Z matrixDist indicesClasificados]=HierarchicalClustering(P,1:numIndividuals(Generacion+1),ramasNegativas);
[eachgroup proto offsetproto sim numGroups] = numberGroup(P,Z,matrixDist,indicesClasificados);

individualsSpecies = eachgroup{numGroups};
protos = proto{numGroups};
offsets = offsetproto{numGroups};


leachgroup = length(eachgroup{numGroups});

nueva = leachgroup + 1;

for i = 1 : leachgroup
    if i == 12
        disp('hola')
    end
    if length(eachgroup{numGroups}{i}) > 1
        dist = [];
        [Z2 matrixDist2 indicesClasificados2] = HierarchicalClusteringGenome(P,poph,Generacion,numIndividuals,eachgroup{numGroups}{i});
        [eachgroup2 proto2 offsetproto2 sim2 numGroups2] = numberGroupGenome(P,Z2,matrixDist2,eachgroup{numGroups}{i});
        if numGroups2>1
            for k = 1 : numGroups2
                dist(k) = 1 - similitud2(protos{i},offsets{i},proto2{numGroups2}{k},offsetproto2{numGroups2}{k});
            end
            [mini pos] = min(dist);
            individualsSpecies{i} = eachgroup2{numGroups2}{pos};
            protos{i} = proto2{numGroups2}{pos};
            offsets{i} = offsetproto2{numGroups2}{pos};
            newSpecies = setdiff(1:numGroups2,pos);
            t = 1;
            for k = nueva:nueva+length(newSpecies)-1
                individualsSpecies{k} = eachgroup2{numGroups2}{newSpecies(t)};
                protos{k}=proto2{numGroups2}{newSpecies(t)};
                offsets{k}= offsetproto2{numGroups2}{newSpecies(t)};

                t = t + 1;
            end
            nueva = nueva+numGroups2-1;
        end
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Z matrixDist indicesClasificados] = HierarchicalClustering(P,individuos,ramasNega)


nonEmptyNeg = true; %no clasifico negativos

raster = P.raster(individuos,:);
dimensions = P.dimensions(individuos,:);
offset = P.offset(individuos,:);
ramasNegativas = ramasNega(individuos);
% las generaciones empiezan en 0, tal y como aparecen en poblacion.txt
numIndividuals = length(raster);
indicesClasificados = [];

if nonEmptyNeg

    negativos = find(ramasNegativas > 0);

    indicesNoClasificados = negativos;
    indicesClasificados =individuos(setdiff(1:numIndividuals,indicesNoClasificados));
    raster(indicesNoClasificados,:)=[];
    dimensions(indicesNoClasificados,:)=[];
    offset(indicesNoClasificados,:)=[];
    numIndividuals = length(raster);
end

dist = zeros(1,(numIndividuals*(numIndividuals-1))/2);
k = 1;
for i = 1 : numIndividuals-1
    i
    raster2 = raster{i};
    offset2 = offset(i,:);
    for j = i+1 : numIndividuals
        raster1 = raster{j};
        offset1 = offset(j,:);
        dist(k) = 1 - similitud2(raster1,offset1,raster2,offset2);
        k = k + 1;
    end
end

matrixDist = squareform(dist);
Z = linkage(dist);
% fig = figure;%('visible', 'off');
[H,T,perm]=dendrogram(Z,0);
set(H,'ButtonDownFcn', {@showInf,raster,dimensions,offset,H,perm,indicesClasificados}, 'Interruptible', 'off');
% set(gca,'XTick',[]);
% saveas(fig,[basedir,'hierarchicalClustering' num2str(Generacion, '%03g')],'png');
% close(fig);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [individualsSpecies proto offsetproto sim numGroups] = numberGroup(P,Z,matrixDist,indicesClasificados)


sZ = size(Z,1);
raster = P.raster;
offset = P.offset;
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

if sim == 0
    derivada = 0;
else
    derivada = abs(diff(sim));
end


%     h1 = figure;%('visible','off');
%
%     bar(sim);
%
%     xlabel('numero de grupos');
%     ylabel('dispersion');
%     saveas(h1,[basedir,'derivada' num2str(Generacion, '%03g')],'png');
%     close(h1);

%
%         h2 = figure;%('visible','off');
%         bar(derivada);
%         xlabel('numero de grupos');
%         ylabel('derivada de la dispersion normalizada');
%         saveas(h2,[basedir,'hierarchicalClustering' num2str(Generacion, '%03g')],'png');
%     close(h2);


%     h3 = figure;%('visible','off');
f = find(derivada <= 1);
if numel(f) == 0
    f = length(derivada)+1;
end

%     paintPrototypes(eachgroup{f(1)},proto{f(1)},offsetproto{f(1)},indicesClasificados)
%     saveas(h3,[basedir,'grupos' num2str(Generacion, '%03g')],'png');
%     close(h3);

numGroups = f(1);


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ances = ancestros(poph,indice,generacion)

numGenerations = max(poph.generation);

numIndividuals = zeros(1,numGenerations);
for i = 1 : numGenerations+1
    numIndividuals(i) = sum(poph.generation == i-1); % numIndividuals(1) = numero de individuos en la generacion 0
end

indicesPredecesores = poph.tree.iA;

predecesores = cell(1,numGenerations);
for i = 2 : numGenerations+1
    predecesores{i-1} = indicesPredecesores(sum(numIndividuals(1:i-1))+1:sum(numIndividuals(1:i)));
end

ances = ones(1,generacion);
for i = generacion:-1:1
    ances(i) = predecesores{i}(indice);
    indice = ances(i);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Z matrixDist indicesClasificados] = HierarchicalClusteringGenome(P,poph,generacion,numIndividuals,individuos)

changes = poph.tree.change;
numIndividual = length(individuos);
for i = 1 : numIndividual
    ances{i} = ancestros(poph,individuos(i),generacion);
end

indicesClasificados = [];

dist = zeros(1,(numIndividual*(numIndividual-1))/2);
k = 1;
for i = 1 : numIndividual-1
    i
    for j = i+1 : numIndividual
        dist(k) = meanChanges(ances{i},ances{j},numIndividuals,changes);
        k = k + 1;
    end
end

matrixDist = squareform(dist);
Z = linkage(dist);
% fig = figure;%('visible', 'off');
% [H,T,perm]=dendrogram(Z,0);
% set(H,'ButtonDownFcn', {@showInf,raster,dimensions,offset,H,perm,indicesClasificados}, 'Interruptible', 'off');
% set(gca,'XTick',[]);
% saveas(fig,[basedir,'hierarchicalClustering' num2str(Generacion, '%03g')],'png');
% close(fig);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = meanChanges(ances1,ances2,numIndividuals,changes)

t = find(ances1 == ances2 == 0);
if isempty(t)
    m = 0;
else
    antecesorComun = t(1)-1;

    s1 = 0;
    s2 = 0;

    for i = antecesorComun : length(ances1)
        if changes{sum(numIndividuals(1:i))+ances1(i)} ~= '='
            s1 = s1 + 1;
        end
        if changes{sum(numIndividuals(1:i))+ances2(i)} ~= '='
            s2 = s2 + 1;
        end
    end

    m = (s1 +s2)/2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [eachgroup proto offsetproto sim numGroups] = numberGroupGenome(P,Z,matrixDist,indicesClasificados)


sZ = size(Z,1);
raster = P.raster;
offset = P.offset;
c = 1;
simi = 1;
if ~isempty(indicesClasificados)
    raster= raster(indicesClasificados,:);
    offset = offset(indicesClasificados,:);
end
while c<=sZ+1 && simi >0

    % groups = cluster(Z,'cutoff',c);
    groups = cluster(Z,'maxclust',c);%,'criterion','distance');
    for i = 1 : c
        ind = find(groups == i);
        numI(i) = length(ind);
        [proto{c}{i} offsetproto{c}{i} meanDisrup{c}{i}]= prototype2(raster,offset,ind,matrixDist);
        eachgroup{c}{i} = indicesClasificados(ind);
    end


    sim(c) = sum([meanDisrup{c}{:}]./numI)/c;
    simi = sim(c);

    c = c + 1;
end


% if sim == 0
%     derivada = 0;
% else
%     derivada = abs(diff(sim));
% end


%     h1 = figure;%('visible','off');
%
%     bar(sim);
%
%     xlabel('numero de grupos');
%     ylabel('dispersion');
%     saveas(h1,[basedir,'derivada' num2str(Generacion, '%03g')],'png');
%     close(h1);

%
%         h2 = figure;%('visible','off');
%         bar(derivada);
%         xlabel('numero de grupos');
%         ylabel('derivada de la dispersion normalizada');
%         saveas(h2,[basedir,'hierarchicalClustering' num2str(Generacion, '%03g')],'png');
%     close(h2);


%     h3 = figure;%('visible','off');
% f = find(derivada <= 10);
% if numel(f) == 0
%     f = length(derivada)+1;
% end

%     paintPrototypes(eachgroup{f(1)},proto{f(1)},offsetproto{f(1)},indicesClasificados)
%     saveas(h3,[basedir,'grupos' num2str(Generacion, '%03g')],'png');
%     close(h3);

f = find(sim <= 10);

numGroups = f(1);