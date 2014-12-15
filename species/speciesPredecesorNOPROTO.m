function ret= speciesPredecesorNOPROTO(basedir,poph,umbral,maxgen, mode, matrixDistsALL, indicesClasificadosALL)

matrixDistEXIST = exist('matrixDistsALL', 'var');

doTesting = false;

switch mode
  case 'edit'
    mode = 1;
  case 'miniedit'
    mode = 2;
  case 'morpho'
    mode = 3;
  otherwise
    error('mode %s not understood', any2str(mode));
end


if isempty(poph)
  if exist([basedir filesep 'poph.mat'], 'file')
    poph = load([basedir filesep 'poph.mat']);
    poph = poph.poph;
  else
    poph = loadSimulation(basedir);
  end
end

poph.tree.parents = poph.tree.parents(:);
if isfield(poph, 'minGenome')
  minGenome = poph.minGenome;
else
  if mode==2
    minGenome = cell(find(poph.generation==maxgen, 1, 'last'), 1);
    gn = poph.genome;
    for j=1:numel(minGenome)
      minGenome{j} = sanitizeString(gn{j});
    end
    clear gn;
  else
    minGenome = [];
  end
end

%WE RELY ON ThESE TWO FACTS TO EASE COMPUTATIONS:
if not(all(poph.idx==poph.tree.idxD(1:numel(poph.idx))))
  error('INDEXES INCOSISTENT FOR %s!!!!', basedir);
end
if not(all(poph.idx(poph.tree.parents(poph.tree.parents>0))==poph.tree.idxA(poph.tree.parents>0)))
  error('PARENT INDEXES INCOSISTENT FOR %s!!!!', basedir);
end

generations = 0:maxgen;

% numGenerations = max(poph.generation);
numSpecies = cell(1,numel(generations));
individualsSpecies = cell(1,numel(generations));
protos = cell(1,numel(generations));
if doTesting
  megaprotos = cell(1,numel(generations));
end

numIndividuals = zeros(1,numel(generations));
initg = zeros(1,numel(generations));
endg  = zeros(1,numel(generations));
for i = 1 : numel(generations)
  thisg = poph.generation == generations(i);
  initg(i) = find(thisg, 1);
  endg(i) = find(thisg, 1, 'last');
    numIndividuals(i) = sum(thisg);
end
clear thisg;

predecesores = poph.tree.parents;
indicesPredecesores = poph.tree.iA;
nPixelsInSoil = poph.nPixelsInSoil;

speciesInd = nan(size(poph.generation));

numSpecies{1}(1) = 1;
individualsSpecies{1}{1}=initg(1):endg(1);%1:numIndividuals(1); % todos los individuos de la generacion 0 pertenecen a la misma especie porque son el mismo individuo
speciesInd(initg(1):endg(1)) = 1; %THE SPECIES INDENTIFICATOR FOR EACH INDIVIDUAL IN THE SIMULATION
P = load(sprintf('%s/P%03d.mat',basedir,0));
P = P.(sprintf('P%03d', 0));
protos{1}(1)=1;
if doTesting
  megaprotos{1}{1}=P.raster{1};
end

Pold = P;
thisgidxOLD = 1;

%LAS ESPECIES VAN AÑADIENDOSE UNA A UNA AL REGISTRO DE ESPECIES. SE NUMERAN
%CONFORME APARECEN


for i = 2 : numel(generations)
    tg = generations(i);
    imd = i-1; %index in matridDists
    ramasNegativas = nPixelsInSoil(initg(i):endg(i));
    idxpredecesores = predecesores(initg(i):endg(i));
    thisgidx = initg(i):endg(i);
    especiesPredecesores = speciesInd(idxpredecesores); 
    especie = cell(size(individualsSpecies{i-1}));
    if doTesting
      especieZ = cell(size(individualsSpecies{i-1}));
      predecesoresANTIGUO = indicesPredecesores(initg(i):endg(i));%(sum(numIndividuals(1:i-1))+1:sum(numIndividuals(1:i)));
    end
    for j=1:numel(individualsSpecies{i-1})
      if not(isempty(individualsSpecies{i-1}{j}))
        especie{j} = find(especiesPredecesores==j); %individuals in  generation i whose ancestors in generation i-1 where in j-th species
      end
      if doTesting
        individuosActuales = [];
        zaq = individualsSpecies{i-1}{j}(:)';
        for k = 1 : length(zaq)
            individuosActuales = [individuosActuales find(predecesoresANTIGUO == zaq(k))']; %#ok<AGROW>
        end
        especieZ{j} = individuosActuales; % hay especie{j} individuos en la generacion actual cuyos predecesores pertenecian a la especie j en la generacion anterior
        if not(isequal(especie{j}(:), especieZ{j}(:)))
          error('jarllll!!!');
        end
      end
    end
    fprintf('GENERATION %d, %d species in last generation\n', tg, numel(numSpecies{i-1}));
    namep = sprintf('P%03d', tg);
    P = load([basedir filesep namep '.mat']);
    P =P.(namep);
    % para cada especie clasificada en la generacion actual teniendo en
    % cuenta las especies de los predecesores, hacemos un dendrograma
    numSpecies{i}= zeros(numel(numSpecies{i-1})+50,1);
    individualsSpecies{i}= cell(numel(individualsSpecies{i-1})+50,1);
    protos{i}= zeros(numel(protos{i-1})+50,1);
    if doTesting
      megaprotos{i}= cell(numel(megaprotos{i-1})+50,1);
    end
    nextSpecie = numel(numSpecies{i-1})+1;
    for j = 1 : numel(numSpecies{i-1})
        if numel(especie{j})>1
          %fprintf('GENERATION %d, SPECIES %d/%d\n', tg,j,numel(numSpecies{i-1}))
            if matrixDistEXIST
              if doTesting
                indicesClasificados2 = getIndicesClasificados(ramasNegativas);
                if not(isequal(indicesClasificados2, indicesClasificadosALL{imd}))
                  error('jarllll!!!');
                end
              end
              [Z matrixDist indicesClasificados] = HierarchicalClustering(especie{j},matrixDistsALL{imd}, indicesClasificadosALL{imd});
            else
              [Z matrixDist indicesClasificados] = HierarchicalClustering2(P,minGenome, thisgidx, especie{j},ramasNegativas, mode);          
            end
            if numel(Z) >0
                [grupos prototypes sim numGroups fullproto] = numberGroup(P,Z,matrixDist,indicesClasificados,umbral,doTesting);
                if numGroups>1
                  dist = zeros(numGroups,1);
                  pold = protos{i-1}(j);
                  if doTesting && not(isequal(Pold.raster{pold}, megaprotos{i-1}{j}))
                    error('jarrlll1');
                  end
                  switch mode
                    case 1
                      for k = 1 : numGroups
                        dist(k) = strdstMX(P.genome{prototypes(k)}, Pold.genome{pold});%strdist(genome{i}, genome{j});
                      end
                    case 2
                      for k = 1 : numGroups
                        dist(k) = strdstMX(minGenome{thisgidx(prototypes(k))}, minGenome{thisgidxOLD(pold)});%strdist(genome{i}, genome{j});
                      end
                    case 3
                      for k = 1 : numGroups
                        p = prototypes(k);
                        if doTesting && not(isequal(P.raster{p}, fullproto{k}))
                          error('jarrlll2');
                        end
                        dist(k) = 1 - similitud2(P.raster{p},P.offset(p,:),Pold.raster{pold}, Pold.offset(pold,:));
                      end
                  end
                  [mini pos] = min(dist);

                  %                 if mini < 1 % la especie continua
                  numSpecies{i}(j)=numSpecies{i-1}(j);
                  individualsSpecies{i}{j} = grupos{pos};
                  speciesInd(thisgidx(grupos{pos})) = j;
                  protos{i}(j) = prototypes(pos);
                  if doTesting
                    megaprotos{i}(j) = fullproto(pos);
                  end
                  %                 end
                  newSpecies = 1:numGroups; newSpecies(pos) = []; %newSpecies = setdiff(1:numGroups,pos);
                  for k = 1:numel(newSpecies) 
                    if numel(numSpecies{i})<nextSpecie
                      numSpecies{i}(end+50) = 0;
                      individualsSpecies{i}{end+50} = [];
                      protos{i}(end+50) = 0;
                      if doTesting
                        megaprotos{i}(end+50) = [];
                      end
                    end
                    numSpecies{i}(nextSpecie) = nextSpecie;
                    individualsSpecies{i}{nextSpecie} = grupos{newSpecies(k)};
                    speciesInd(thisgidx(grupos{newSpecies(k)})) = nextSpecie;
                    protos{i}(nextSpecie) = prototypes(newSpecies(k));
                    if doTesting
                      megaprotos{i}(nextSpecie) = fullproto(newSpecies(k));
                    end
                    nextSpecie = nextSpecie+1;
                  end
                else
                  if (numel(grupos)>1) || (numel(prototypes)>1)
                    error('requetejarlllll!!!!');
                  end
                  numSpecies{i}(j)=numSpecies{i-1}(j);
                  individualsSpecies{i}{j} = grupos{1};
                  speciesInd(thisgidx(grupos{1})) = j;
                  protos{i}(j) = prototypes(1);
                  if doTesting
                    megaprotos{i}(j) = fullproto(1);
                  end
                end
            else
                if numel(indicesClasificados)>0
                    numSpecies{i}(j) = numSpecies{i-1}(j);
                    individualsSpecies{i}{j} = indicesClasificados;
                    speciesInd(thisgidx(indicesClasificados)) = j;
                    protos{i}(j) = indicesClasificados(1);
                    if doTesting
                      megaprotos{i}{j} = P.raster{indicesClasificados};
                    end
                else
                    numSpecies{i}(j) = 0;
                    individualsSpecies{i}{j} = [];
                    protos{i}(j) = 0;
                    if doTesting
                      megaprotos{i}{j} = [];
                    end
                end
            end
        elseif numel(especie{j}) == 1
            numSpecies{i}(j) = numSpecies{i-1}(j);
            individualsSpecies{i}{j} = especie{j};
            speciesInd(thisgidx(especie{j})) = j;
            protos{i}(j) = especie{j};
            if doTesting
              megaprotos{i}{j} = megaprotos{i-1}{j};
            end
        end
        
    end
    Pold = P;
    thisgidxOLD = thisgidx;
    numSpecies{i}= numSpecies{i}(1:nextSpecie-1);
    individualsSpecies{i}= individualsSpecies{i}(1:nextSpecie-1);
    protos{i}= protos{i}(1:nextSpecie-1);
    if doTesting
      megaprotos{i}= megaprotos{i}(1:nextSpecie-1);
    end
    
    %s=0;for p = 1:numel(individualsSpecies{i});s= s + numel(individualsSpecies{i}{p});end
end

ret = struct('numSpecies', {numSpecies}, 'individualsSpecies', {individualsSpecies}, 'protos', {protos}, 'speciesInd', {speciesInd});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Z matrixDist indicesClasificados]=HierarchicalClustering(especie,matrixdGENERAL, indicesclaGENERAL)%
%,basedir,Generacion)

[esclasif idxinclasif] = ismember(especie, indicesclaGENERAL);

idxinclasif = idxinclasif(esclasif);

indicesClasificados = especie(esclasif);

matrixDist = matrixdGENERAL(idxinclasif, idxinclasif);

if isempty(matrixDist) || (numel(idxinclasif)==1)
    Z = [];
else
    Z = linkage(squareform(matrixDist));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Z matrixDist indicesClasificados]=HierarchicalClustering2(P,minGenome, thisg, especie,ramasNega, mode)%
%,basedir,Generacion)

nonEmptyNeg = true; %no clasifico negativos

switch mode
  case 1
    genomes = minGenome(thisg);
  case 2
    genomes = P.genome(especie,:);
  case 3
    raster = P.raster(especie,:);
    offset = P.offset(especie,:);
end
ramasNegativas = ramasNega(especie);
% las generaciones empiezan en 0, tal y como aparecen en poblacion.txt
numIndividuals = numel(especie);
indicesClasificados = [];

if nonEmptyNeg

    negativos = find(ramasNegativas > 0);

    indicesNoClasificados = negativos;
    indicesClasificados = especie(setdiff(1:numIndividuals,indicesNoClasificados));
    switch mode
      case {1, 2}
        genomes(indicesNoClasificados)=[];
      case 3
        raster(indicesNoClasificados,:)=[];
        offset(indicesNoClasificados,:)=[];
    end
    numIndividuals = numel(indicesClasificados);
end
dist = zeros(1,(numIndividuals*(numIndividuals-1))/2);

switch mode
  case {1, 2}
    k=1;
    for i=1:numIndividuals-1
      genome2 = genomes{i};
      for j=i+1:numIndividuals
        genome1 = genomes{j};
        dist(k) = strdstMX(genome1, genome2);
        k=k+1;
      end
    end
  case 3
    k = 1;
    %dist2 = dist;
    for i = 1 : numIndividuals-1
        raster2 = raster{i};
        offset2 = offset(i,:);
        for j = i+1 : numIndividuals
            raster1 = raster{j};
            offset1 = offset(j,:);
            dist(k) = 1 - similitud2(raster1,offset1,raster2,offset2);
            %dist2(k) = 1-similitud2OLD(raster1,offset1,raster2,offset2);
            k = k + 1;
        end
    end
end

matrixDist = squareform(dist);
%matrixDist2 = squareform(dist2);
if numel(dist) == 0
    Z = [];
else
    Z = linkage(dist);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sim = similitud2(raster1,offset1,raster2,offset2)%, show)
    pixel1 = [raster1{:}];
    pixel2 = [raster2{:}];
    v = isempty(pixel1) + isempty(pixel2);
if v==2
    sim = 1;
elseif v==1
    sim = 0;
else
    
    if size(pixel1,1) == 1 
        pixel1 = [];

        pixel1(:,1) = raster1{1};
        pixel1(:,2) = raster1{2};
    end
    if size(pixel2,1) == 1 
        pixel2 = [];

        pixel2(:,1) = raster2{1};
        pixel2(:,2) = raster2{2};
    end
    
    pixel1(:,1) = pixel1(:,1) - offset1(1);
    pixel1(:,2) = pixel1(:,2) - offset1(2);

    pixel2(:,1) = pixel2(:,1) - offset2(1);
    pixel2(:,2) = pixel2(:,2) - offset2(2);
    
    sim = similitudFast(pixel1, pixel2);
    %sim = size(intersect(pixel1,pixel2,'rows'),1)/size(union(pixel1,pixel2,'rows'),1);
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
function [eachgroupO protoO sim numGroups fullprotoO offsetprotoO] = numberGroup(P,Z,matrixDist,indicesClasificados,umbral,doTesting)
%,basedir,Generacion)
    if doTesting
      raster = P.raster;
      offset = P.offset;
      if ~isempty(indicesClasificados)
          raster= raster(indicesClasificados,:);
          offset = offset(indicesClasificados,:);
      end
    end
    % Z = linkage(dist);
    sZ = size(Z,1);
    c = 1;
    sim = nan;
    firstone = true;
    continuar = true;
    goon = c<=sZ+1 && continuar;
    while goon

        % groups = cluster(Z,'cutoff',c);
        groups = cluster(Z,'maxclust',c,'criterion','distance');
        proto = zeros(c,1);
        meanDisrup = zeros(c,1); 
        eachgroup = cell(c,1);
        if doTesting
          meanDisrup2 = zeros(c,1);
          fullproto = cell(c,1);
          offsetproto = cell(c,1);
        end
        for i = 1 : c
            eachgroup{i}= find(groups == i);
            [proto(i) meanDisrup(i)]= prototype2(eachgroup{i},matrixDist);
            if doTesting
              [fullproto{i} offsetproto{i} meanDisrup2(i)]= prototype2FULL(raster,offset,eachgroup{i},matrixDist);
              if meanDisrup(i)~=meanDisrup2(i)
                error('requete!!!');
              end
            end
            eachgroup{i} = indicesClasificados(eachgroup{i});
        end

        oldsim = sim;
        sim = sum(meanDisrup)/c;
        derivada = abs(oldsim-sim);
        continuar = (sim>0) && (firstone || (derivada>umbral));
        c = c + 1;
        goon = c<=sZ+1 && continuar;
        if firstone || goon
          oldeachgroup = eachgroup;
          oldproto = proto;
          if doTesting
            oldfullproto = fullproto;
            oldoffsetproto = offsetproto;
          end
        end
        firstone = false;
    end
    
    numGroups = numel(oldeachgroup);
    protoO = indicesClasificados(oldproto);
    eachgroupO = oldeachgroup;
    if doTesting
      offsetprotoO = oldoffsetproto;
      fullprotoO = oldfullproto;
    else
      offsetprotoO = [];
      fullprotoO = [];
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [minidisruptive mini]= prototype2(eachgroup,matrixDist)

% prototype = less disruptive individual of the group
leachgroup=length(eachgroup);

disruptive = zeros(leachgroup,1);
for i = 1 : leachgroup
    disruptive(i) = sum(matrixDist(eachgroup(i),eachgroup));
end

[mini pos]= min(disruptive);
minidisruptive = eachgroup(pos);

% meanDisrup = sum(sum(matrixDist(eachgroup,eachgroup)))/leachgroup^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [proto offsetproto mini]= prototype2FULL(raster,offset,eachgroup,matrixDist)

% prototype = less disruptive individual of the group
leachgroup=length(eachgroup);

for i = 1 : leachgroup
    disruptive(i) = sum(matrixDist(eachgroup(i),eachgroup));
end

[mini pos]= min(disruptive);
minidisruptive = eachgroup(pos);
proto = raster{minidisruptive};
offsetproto = offset(minidisruptive,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dist = getDistance(P, p, Pold, pold, mode)
switch mode
  case 1
    dist = strdstMX(P.genome{p}, Pold.genome{pold});%strdist(genome{i}, genome{j});
  case 2
    dist = strdstMX(P.minGenome{p}, Pold.minGenome{pold});%strdist(genome{i}, genome{j});
  case 3
    dist = 1 - similitud2(P.raster{p},P.offset(p,:),Pold.raster{pold}, Pold.offset(pold,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sim = similitud2OLD(raster1,offset1,raster2,offset2)%, show)

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

function indicesClasificados = getIndicesClasificados(ramasNegativas)
numIndividual = numel(ramasNegativas);
        negativos = find(ramasNegativas > 0);

        indicesNoClasificados = negativos;
        indicesClasificados = setdiff(1:numIndividual,indicesNoClasificados);

