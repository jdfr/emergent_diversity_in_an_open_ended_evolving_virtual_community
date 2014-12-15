function HierarchicalClusteringEDITDISTANCEMINI(basedir,poph,GenIn,GenFin,parallel,prefix)

if not(exist('parallel', 'var'))
  parallel = false;
end
if not(exist('prefix', 'var'))
  prefix = '';
end

% poph = loadSimulation(basedir);

if parallel
  jobArgs={'Tag', 'clusteringEDITDISTANCE', 'PathDependencies',  {pwd}};
  gens = array2cell((GenIn:GenFin)');
  pophs = cell(size(gens));
  basedirs = repmat({basedir}, size(gens));
  prefixes = repmat({prefix}, size(gens));
  dummy = my_dfeval(@doAGeneration, pophs, basedirs, gens, prefixes, jobArgs{:});
else
  for a = GenIn : GenFin
    doAGeneration(poph, basedir, a, prefix);
  end
end

close all;

function dummy = doAGeneration(poph, basedir, a, prefix)
dummy = a;
    tosave1 = [prefix 'Z' num2str(a, '%03g')];
    tosave2 = [prefix 'matrixDist' num2str(a, '%03g')];
    tosave3 = [prefix 'indicesClasificados' num2str(a, '%03g')];
    if exist([basedir filesep tosave1], 'file') && exist([basedir filesep tosave2], 'file') && exist([basedir filesep tosave3], 'file')
      return
    end
%     numIndividuals = zeros(1,a+1);
%     for i = 1 : a+1
%         numIndividuals(i) = sum(poph.generation == i-1); % numIndividuals(1) = numero de individuos en la generacion 0
%     end
% 
if isempty(poph)
  if exist([basedir filesep 'poph.mat'], 'file')
    poph = load([basedir filesep 'poph.mat']);
    poph = poph.poph;
  else
    poph = loadSimulation(basedir, false, false, false);
    poph.minGenome = cellfun(@sanitizeString, poph.genome, 'uniformoutput', false);
  end
end
nPixelsInSoil = poph.nPixelsInSoil;
    ramasNegativas = nPixelsInSoil(poph.generation==a);%sum(numIndividuals(1:a))+1:sum(numIndividuals(1:a+1)));

    nonEmptyNeg = true; %no clasifico vacios y negativos

    genome = poph.minGenome(poph.generation==a);
    % las generaciones empiezan en 0, tal y como aparecen en poblacion.txt
    numIndividual = numel(genome);
    indicesClasificados = [];

    if nonEmptyNeg

        negativos = find(ramasNegativas > 0);

        indicesNoClasificados = negativos;
        indicesClasificados = setdiff(1:numIndividual,indicesNoClasificados);
        genome(indicesNoClasificados,:)=[];
        numIndividual = length(genome);
    end

    dist = zeros(1,(numIndividual*(numIndividual-1))/2);
    k = 1;
    if numIndividual > 1
        for i = 1 : numIndividual-1
            i
            for j = i+1 : numIndividual
                dist(k) = strdist(genome{i}, genome{j});
                k = k + 1;
            end
        end

        matrixDist = squareform(dist);
        Z = linkage(dist);
        % fig = figure;%('visible', 'off');
        %[H,T,perm]=dendrogram(Z,0);
        %set(H,'ButtonDownFcn', {@showInf,raster,dimensions,offset,H,perm,indicesClasificados}, 'Interruptible', 'off');
        % set(gca,'XTick',[]);
    else
        Z = [];
        matrixDist = [];
    end
    % saveas(fig,[basedir,'hierarchicalClustering' num2str(Generacion, '%03g')],'png');
    % close(fig);

    eval([tosave1 ' = Z;']);
    eval([tosave2 ' = matrixDist;']);
    eval([tosave3 ' = indicesClasificados;']);

    save([basedir filesep tosave1],tosave1);
    save([basedir filesep tosave2],tosave2);
    save([basedir filesep tosave3],tosave3);

    clear('Z', 'matrixDist','indicesClasificados', tosave1,tosave2,tosave3);

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
return
if plotTree
    if isempty(indicesClasificados)
        if isempty(raster{perm(cx)}{1})
            sprintf('empty individual')
        else
            fig = figure;%('visible', 'off');
            drawtree(raster{perm(cx)}, dimension(perm(cx),:), offset(perm(cx),:), 0, 1, [0 0 0], 'fit');
            title(sprintf('individual: %03d', perm(cx)));
        end
    else
        if isempty(raster{perm(cx)}{1})
            sprintf('empty individual')
        else
            fig = figure;%('visible', 'off');
            drawtree(raster{perm(cx)}, dimension(perm(cx),:), offset(perm(cx),:), 0, 1, [0 0 0], 'fit');
            title(sprintf('individual: %03d', indicesClasificados(perm(cx))));
        end
    end
end