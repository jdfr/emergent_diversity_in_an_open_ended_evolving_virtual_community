function HierarchicalClusteringNUMMUTS(basedir,poph,GenIn,GenFin, parallel)

if not(exist('parallel', 'var'))
  parallel = false;
end

if isempty(poph)
  if exist([basedir filesep 'poph.mat'], 'file');
    poph = load([basedir filesep 'poph.mat']);
    poph = poph.poph;
  else
    poph = loadSimulation(basedir);
  end
end
% poph = loadSimulation(basedir);

if parallel
  poph =  makeNumMutations(poph);
  numMuts = poph.tree.numMuts;
  mutDists = poph.tree.mutDists;
  save([basedir filesep 'mutdists.mat'], 'numMuts', 'mutDists');
  jobArgs={'Tag', 'clusteringNUMMUTS', 'PathDependencies',  {pwd}};
  gens = array2cell((GenIn:GenFin)');
  pophs = cell(size(gens));
  basedirs = repmat({basedir}, size(gens));
  dummy = my_dfeval(@doAGeneration, pophs, basedirs, gens, jobArgs{:});
else
  poph =  makeNumMutations(poph);
  for a = GenIn : GenFin
    doAGeneration(poph, basedir, a);
  end
end

function dummy = doAGeneration(poph, basedir, a)
dummy = a;
    tosave1 = ['Z' num2str(a, '%03g')];
    tosave2 = ['matrixDist' num2str(a, '%03g')];
    tosave3 = ['indicesClasificados' num2str(a, '%03g')];
    if exist([basedir filesep tosave1], 'file') && exist([basedir filesep tosave2], 'file') && exist([basedir filesep tosave3], 'file')
      return
    end
if isempty(poph)
  poph = loadSimulation(basedir);
  zs = load([basedir filesep 'mutdists.mat']);
  poph.tree.numMuts = zs.numMuts;
  poph.tree.mutDists = zs.mutDists;
  clear zs;
end

% poph = loadSimulation(basedir);
nPixelsInSoil = poph.nPixelsInSoil;

%     numIndividuals = zeros(1,a+1);
%     for i = 1 : a+1
%         numIndividuals(i) = sum(poph.generation == i-1); % numIndividuals(1) = numero de individuos en la generacion 0
%     end
% 
    ramasNegativas = nPixelsInSoil(poph.generation==a);%sum(numIndividuals(1:a))+1:sum(numIndividuals(1:a+1)));

    nonEmptyNeg = true; %no clasifico vacios y negativos

    genome = poph.genome(poph.generation==a);
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
                dist(k) = mutDist(poph, a, indicesClasificados(i), indicesClasificados(j));
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

close all;

function poph = makeNumMutations(poph)

numMuts = zeros(size(poph.tree.change));
change = poph.tree.change;
for k=1:numel(numMuts)
  c = upper(change{k});
  numMuts(k) = sum((c>='A')&(c<='Z')&(c~='G'));
end
poph.tree.numMuts = numMuts;
% cumNumMuts = zeros(size(poph.tree.change));
% for k=1:numel(cumNumMuts)
%   if poph.tree.idxA(k)~=0
%     cumNumMuts(k) = cumNumMuts(find(poph.tree.idxD==poph.tree.idxA(k)))+numMuts(k);
%   end
% end
% poph.tree.cumNumMuts = cumNumMuts;
mutDists = cell(max(poph.generation)+1,1);
gen=0;
mutDists{1} = nan(sum(poph.generation==gen));
for z=1:size(mutDists{1},1)
  mutDists{1}(z,z) = 0;
end
for k=2:numel(mutDists)
  gen = gen+1;
  mm = mutDists{k-1};
  thisGen = find(poph.tree.gD==gen);
  m = zeros(numel(thisGen));
  z = size(m,1);
  inds=poph.tree.iD(thisGen);
  if not(all(inds==(1:z)'))
    error('not consecutive!!\n%s', mat2str(inds));
  end
  for a=1:z
    for b=1:z
      m(a,b) = mm(poph.tree.iA(thisGen(a)), poph.tree.iA(thisGen(b)))+numMuts(thisGen(a))+numMuts(thisGen(b));
    end
    m(a,a)=0;
  end
  mutDists{k}=m;
end
poph.tree.mutDists=mutDists;

function d=mutDist(poph, gen, idx1, idx2)
d=poph.tree.mutDists{gen+1}(idx1, idx2);

function d=mutDist2(poph, gen, idx1, idx2)
n1 = 0;
n2 = 0;
i1 = find((poph.tree.gD==gen)&(poph.tree.iD==idx1));
i2 = find((poph.tree.gD==gen)&(poph.tree.iD==idx2));
while i1~=i2
  n1 = n1+poph.tree.numMuts(i1);
  n2 = n2+poph.tree.numMuts(i2);
  i1 = find(poph.tree.idxD==poph.tree.idxA(i1));
  i2 = find(poph.tree.idxD==poph.tree.idxA(i2));
end
d=n1+n2;

  