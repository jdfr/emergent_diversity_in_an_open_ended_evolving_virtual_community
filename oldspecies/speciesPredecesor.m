function [numSpecies individualsSpecies protos offsets]= speciesPredecesor(basedir,poph, umbral,numGenerations)

% method = 1, elbow[
% method = 2, maximizar funcion
% umbral = 1;
method = 1;
if isempty(poph)
  poph = loadSimulation(basedir);
end

% numGenerations = max(poph.generation);
numSpecies = cell(1,numGenerations);
individualsSpecies = cell(1,numGenerations);
protos = cell(1,numGenerations);
offsets = cell(1,numGenerations);

numIndividuals = zeros(1,numGenerations+1);
for i = 1 : numGenerations+1
    numIndividuals(i) = sum(poph.generation == i-1); % numIndividuals(1) = numero de individuos en la generacion 0
end

indicesPredecesores = poph.tree.iA;
nPixelsInSoil = poph.nPixelsInSoil;

numSpecies{1}(1) = 1;
individualsSpecies{1}{1}=1:numIndividuals(1); % todos los individuos de la generacion 0 pertenecen a la misma especie porque son el mismo individuo
load(sprintf('%s/P%03d.mat',basedir,0));
P=eval(sprintf('P%03d', 0));
protos{1}{1}=P.raster{1};
offsets{1}{1}= P.offset(1,:);

for i = 2 : numGenerations+1
    if i == 339
        disp('hola')
    end
    ramasNegativas = nPixelsInSoil(sum(numIndividuals(1:i-1))+1:sum(numIndividuals(1:i)));
  
    i
    sg = [];
    predecesores = indicesPredecesores(sum(numIndividuals(1:i-1))+1:sum(numIndividuals(1:i)));
    for j = 1 : length(numSpecies{i-1})
        individuosActuales = [];
        for k = 1 : length(individualsSpecies{i-1}{j})
            individuosActuales = [individuosActuales find(predecesores == individualsSpecies{i-1}{j}(k))'];
        end
        especie{j} = individuosActuales; % hay especie{j} individuos en la generacion actual cuyos predecesores pertenecian a la especie j en la generacion anterior
    end
    lVars = load(sprintf('%s/P%03d.mat',basedir,i-1));
    P=eval(sprintf('lVars.P%03d', i-1));
    % para cada especie clasificada en la generacion actual teniendo en
    % cuenta las especies de los predecesores, hacemos un dendrograma
    nueva = length(numSpecies{i-1});
    for j = 1 : length(numSpecies{i-1})
        if j == 369
            disp('hola')
        end
        Z = [];
        matrixDist = [];
        indicesClasificados = [];
        eachgroup = [];
        proto = [];
        offsetproto = [];
        sim = [];
        numGroups = [];
        if numel(especie{j})>1
            dist = [];
            [Z matrixDist indicesClasificados] = HierarchicalClustering(P,especie{j},ramasNegativas);
            if numel(Z) >0
                if method == 1
                    [eachgroup proto offsetproto sim numGroups] = numberGroup(P,Z,matrixDist,indicesClasificados,umbral);
                else
                    [eachgroup proto offsetproto C numGroups] = numberGroupCalinski(P,Z,matrixDist,indicesClasificados);
                end
                prototypes = proto{numGroups};
                offsetprototypes = offsetproto{numGroups};
                grupos = eachgroup{numGroups};
                for k = 1 : numGroups
                    dist(k) = 1 - similitud2(prototypes{k},offsetprototypes{k},protos{i-1}{j},offsets{i-1}{j});
                end
                [mini pos] = min(dist);

                %                 if mini < 1 % la especie continua
                numSpecies{i}(numSpecies{i-1}(j))=numSpecies{i-1}(j);
                individualsSpecies{i}{numSpecies{i-1}(j)} = grupos{pos};
                protos{i}{numSpecies{i-1}(j)} = prototypes{pos};
                offsets{i}{numSpecies{i-1}(j)} = offsetprototypes{pos};
                %                 end
                newSpecies = setdiff(1:numGroups,pos);
                t = 1;
                for k = nueva+1:nueva+numGroups-1
                    numSpecies{i}(k)=k;
                    individualsSpecies{i}{k} = grupos{newSpecies(t)};
                    protos{i}{k}=prototypes{newSpecies(t)};
                    offsets{i}{k}= offsetprototypes{newSpecies(t)};

                    t = t + 1;
                end
            nueva = nueva+numGroups-1;
            else
                if numel(indicesClasificados)>0
                    numSpecies{i}(numSpecies{i-1}(j)) = numSpecies{i-1}(j);
                    individualsSpecies{i}{numSpecies{i-1}(j)} = indicesClasificados;
                    protos{i}{numSpecies{i-1}(j)} = P.raster{indicesClasificados};
                    offsets{i}{numSpecies{i-1}(j)} = P.offset(indicesClasificados,:);
                else
                    numSpecies{i}(numSpecies{i-1}(j)) = 0;
                    individualsSpecies{i}{numSpecies{i-1}(j)} = indicesClasificados;
                    protos{i}{numSpecies{i-1}(j)} = indicesClasificados;
                    offsets{i}{numSpecies{i-1}(j)} = indicesClasificados;

                end
            end
        elseif numel(especie{j}) == 1
            numSpecies{i}(numSpecies{i-1}(j)) = numSpecies{i-1}(j);
            individualsSpecies{i}{numSpecies{i-1}(j)} = especie{j};
            protos{i}{numSpecies{i-1}(j)} = P.raster{especie{j}};
            offsets{i}{numSpecies{i-1}(j)} =P.offset(especie{j},:);
        end
        


    end
    s=0;for p = 1:length(individualsSpecies{i});s= s + length(individualsSpecies{i}{p});end
end

% h1 = figure;
% for i = 0 : numGenerations - 1
%     numSpecies{i+1}(find(numSpecies{i+1}==0))=[];
% numSpeciesperG(i) = length(find(numSpecies{i+1}>0))
%     plot(i,numSpecies{i+1},'b*');
%     hold on;
% end
% xlabel('generations')
% ylabel('species')

% h2 = figure;
% plot(numSpeciesperG);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Z matrixDist indicesClasificados]=HierarchicalClustering(P,especie,ramasNega)%
%,basedir,Generacion)

ramasNegativas = ramasNega(especie);
nonEmptyNeg = true; %no clasifico negativos

raster = P.raster(especie,:);
dimensions = P.dimensions(especie,:);
offset = P.offset(especie,:);
ramasNegativas = ramasNega(especie);
% las generaciones empiezan en 0, tal y como aparecen en poblacion.txt
numIndividuals = length(raster);
indicesClasificados = [];

if nonEmptyNeg

    negativos = find(ramasNegativas > 0);

    indicesNoClasificados = negativos;
    indicesClasificados = especie(setdiff(1:numIndividuals,indicesNoClasificados));
    raster(indicesNoClasificados,:)=[];
    dimensions(indicesNoClasificados,:)=[];
    offset(indicesNoClasificados,:)=[];
    numIndividuals = length(raster);
end
dist = zeros(1,(numIndividuals*(numIndividuals-1))/2);
k = 1;
for i = 1 : numIndividuals-1
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
if numel(dist) == 0
    Z = [];
else
    Z = linkage(dist);
end
% fig = figure('visible', 'off');
% [H,T,perm]=dendrogram(Z,0);
% set(H,'ButtonDownFcn', {@showInf,raster,dimensions,offset,H,perm,indicesClasificados}, 'Interruptible', 'off');
% set(gca,'XTick',[]);
% saveas(fig,[basedir,'hierarchicalClustering' num2str(Generacion, '%03g')],'png');
% close(fig);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eachgroup proto offsetproto sim numGroups] = numberGroup(P,Z,matrixDist,indicesClasificados,umbral)
%,basedir,Generacion)
    % Z = linkage(dist);
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
            eachgroup{c}{i} = indicesClasificados(eachgroup{c}{i});
        end


        sim(c) = sum([meanDisrup{i}{:}])/c;
        simi = sim(c);

        c = c + 1;
    end
    
    if sim == 0
        derivada = 0;
    else
        derivada = abs(diff(sim));
    end


%     h1 = figure('visible','off');
%     
%     bar(sim);
% 
%     xlabel('numero de grupos');
%     ylabel('dispersion');
% %     saveas(h1,[basedir,'derivada' num2str(Generacion, '%03g')],'png');
% %     close(h1);
% 
% 
%         h2 = figure('visible','off');
%         bar(derivada);
%         xlabel('numero de grupos');
%         ylabel('derivada de la dispersion normalizada');
% %         saveas(h2,[basedir,'hierarchicalClustering' num2str(Generacion, '%03g')],'png');
% %     close(h2);
% 
% 
%     h3 = figure('visible','off');
    f = find(derivada <= umbral);
    if numel(f) == 0
        f = length(derivada)+1;
    end

%     paintPrototypes(eachgroup{f(1)},proto{f(1)},offsetproto{f(1)},indicesClasificados)
%     saveas(h3,[basedir,'grupos' num2str(Generacion, '%03g')],'png');
%     close(h3);
  
    numGroups = f(1);

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



