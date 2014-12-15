function [comGs comLZgenome comLZphenotype comLZfigure comLZgenomeLongitud comLZphenotypeLongitud comLZfigureLongitud genomeI]= evolutionComplexity(basedir,Generacion,individuo)

poph = loadSimulation(basedir);
load(sprintf('%s/P%03d.mat',basedir,Generacion));
P=eval(sprintf('P%03d', Generacion));

numIndividuals = zeros(1,Generacion);
for i = 1 : Generacion+1
    numIndividuals(i) = sum(poph.generation == i-1); % numIndividuals(1) = numero de individuos en la generacion 0
end

ances = ancestros(poph,individuo,Generacion,numIndividuals);
indices = [ances individuo];
for i = 1 : Generacion + 1
    i
    [comGs(i) comLZgenome(i) comLZphenotype(i) comLZfigure(i) comLZgenomeLongitud(i) comLZphenotypeLongitud(i) comLZfigureLongitud(i) genomeI{i}] = complexity(poph,basedir,indices(i),i-1);
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function ances = ancestros(poph,indice,generacion,numIndividuals)

indicesPredecesores = poph.tree.iA;

predecesores = cell(1,generacion);
for i = 2 : generacion+1  
    predecesores{i-1} = indicesPredecesores(sum(numIndividuals(1:i-1))+1:sum(numIndividuals(1:i)));
end

ances = ones(1,generacion);
for i = generacion:-1:1
    ances(i) = predecesores{i}(indice);
    indice = ances(i);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [comGs comLZgenome comLZphenotype comLZfigure comLZgenomeLongitud comLZphenotypeLongitud comLZfigureLongitud genomeI] = complexity(poph,basedir,individuo,generacion)

lVars = load(sprintf('%s/P%03d.mat',basedir,generacion));
P=eval(sprintf('lVars.P%03d', generacion));

raster = P.raster{individuo};
dimension = P.dimensions(individuo,:);
offset = P.offset(individuo,:);

numPixeles = size(raster{1},1);

altura = poph.height;
anchura = poph.width;
numRamas = poph.nbranches;
genome = poph.genome;

numGenerations = max(poph.generation);
numIndividuals = zeros(1,numGenerations);
for i = 1 : numGenerations+1
    numIndividuals(i) = sum(poph.generation == i-1); % numIndividuals(1) = numero de individuos en la generacion 0
end

if generacion == 0
    alturaI = altura(individuo);
    anchuraI = anchura(individuo);
    numRamasI = numRamas(individuo);
    genomeI = genome(individuo);
else
    alturaI = altura(sum(numIndividuals(1:generacion))+individuo);
    anchuraI = anchura(sum(numIndividuals(1:generacion))+individuo);
    numRamasI = numRamas(sum(numIndividuals(1:generacion))+individuo);
    genomeI = genome(sum(numIndividuals(1:generacion))+individuo);
end

if numRamasI == 0
    comGs = 0;
    comLZgenomeLongitud = 0;
    comLZgenome = 0;
    comLZphenotypeLongitud = 0;
    comLZphenotype = 0;
    comLZfigureLongitud = 0;
    comLZfigure = 0;
else

    % complexity depending on the number of Gs
    comGs = numRamasI;

    % Lempel-Ziv complexity of the genome
    sgenome = length(uint8(genomeI{1}));
    LZgenome = norm2lzw(uint8(genomeI{1}));
    sLZgenome = length(LZgenome);
    comLZgenomeLongitud = sLZgenome;
    comLZgenome = (sLZgenome*100)/sgenome;


    % Lempel-Ziv complexity of the phenotype (three iterations of the genome)
    phenotype = ls2(genomeI{1},'canonical');
    sphenotype = length(uint8(phenotype));
    LZphenotype = norm2lzw(uint8(phenotype));
    sLZphenotype = length(LZphenotype);
    comLZphenotypeLongitud = sLZphenotype;
    comLZphenotype = (sLZphenotype*100)/sphenotype;

    % Lempel-Ziv complexity of the phenotype (as a figure)
    [colors array] = drawtree(raster,dimension,offset,0,1,[0 0 0]);
    sfigure = length(array(:));
    LZfigure = norm2lzw(uint8(array(:)));
    sLZfigure = length(LZfigure);
    comLZfigureLongitud = sLZfigure;
    comLZfigure = (sLZfigure*100)/sfigure;
end

% title({['COMPLEXITYGs: ',mat2str(comGs)],['COMPLEXITYLZgenome: ',mat2str(comLZgenome)],['COMPLEXITYLZphenotype: ',mat2str(comLZphenotype)],['COMPLEXITYLZfigure: ',mat2str(comLZfigure)],['COMPLEXITYLZgenomeLongitud: ',mat2str(comLZgenomeLongitud)],['COMPLEXITYLZphenotypeLongitud: ',mat2str(comLZphenotypeLongitud)],['COMPLEXITYLZfigureLongitud: ',mat2str(comLZfigureLongitud)],['DIMENSION: ',mat2str(dimension)]},'FontSize',6);
% drawnow;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function axiom = ls2(regla,T)

if(~corchetesBalanceados(regla))
    error('ls2: regla con corchetes no balanceados: %s', regla);
end

ProductionSystem = ['G=>',regla];

if any(upper(ProductionSystem)~=ProductionSystem)
    error('This code only can handle production systems whose letters are uppercase');
end

rules = regexp(lower(ProductionSystem), '([^ ])=>([^ ]*)', 'tokens'); %1xnrules cell array of 1x2 cell arrays
n_Rules = numel(rules);

rules = vertcat(rules{:});
pres  = upper(horzcat(rules{:,1})); %precedents  are uppercase
randomExpression = isstruct(T) && (isfield(T, 'randomize')) && logical(T.randomize);
posts = cellfun(@(x)sanitizeString(x,randomExpression), rules(:,2), 'uniformoutput', false);                 %consequents are lowercase
if ~all(isletter(pres)) || any(pres=='�') || any(pres=='�') % '�' and '�' are deemed letters by matlab, but they have not uppercase counterparts
  error('This code only can handle rules whose precedent are letters!');
end

% starting string (axiom)
axiom = 'G';

n_Repeats = 3;

for i=1:n_Repeats
  for j=1:n_Rules
    axiom=strrep(axiom, pres(j), posts{j});
  end
  axiom=upper(axiom);
end