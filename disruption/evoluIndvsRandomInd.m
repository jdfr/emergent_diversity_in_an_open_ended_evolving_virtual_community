function evoluIndvsRandomInd(basedirs, symbols,finaldir)
% Estudio de la robustez: disrupción de los operadores de mutación en
% individuos aleatorios vs evolucionados

opFuncs = {@op_alteracionG, ...
           @op_eliminacionG, ...
           @op_insercionG, ...
           @op_dupaleatoriaG, ...
           @op_duplAleatoriaEdit, ...
           @op_duplSecuenciaEdit};

plotXTickLabels = 'A|D|I|DR|DRE|DTE';

numOpFuncs = size(opFuncs, 2);

total = size(basedirs,2);

for i = 1 : total;
    poph{i} = loadSimulation(basedirs{i});
    numGenerations = max(poph{i}.generation);
    genom = poph{i}.genome(poph{i}.generation == numGenerations);
    s = size(genom,1);
    genEvoR = randperm(s);
    for j = 1 : 20
        genome{i}{j} = genom{genEvoR(j)};
    end
    sgenome = size(genome{i},2);
    numGs = zeros(1, sgenome);
    numBra = zeros(1, sgenome);
    numSillyBra = zeros(1, sgenome);
    numAnid = zeros(1, sgenome);
    longGenome = zeros(1, sgenome);
    for j = 1 : sgenome
        numGs(j) = sum(genome{i}{j} == 'G');
        numBra(j) = sum(genome{i}{j} == '[');
        numSillyBra(j) = sum(genome{i}{j}(1:end-1) == '[' & genome{i}{j}(2:end) == ']');
        numAnid(j) = max(cumsum((genome{i}{j} == '[') + (-1 * (genome{i}{j} == ']'))));
        longGenome(j) = length(genome{i}{j});
    end
    
    meanGs(i) = mean(numGs);
    meanBra(i) = mean(numBra);
    meanSillyBra(i) = mean(numSillyBra);
    meanAnid(i) = mean(numAnid);
    meanLongGenome(i) = mean(longGenome);

    if(meanAnid > 1)
        genomeRandom{i} = generateGenomeRTipo2(sgenome,round(meanGs(i)),round(meanLongGenome(i)),round(meanBra(i)), round(meanAnid(i)));
    else
        genomeRandom{i} = generateGenomeRTipo3(sgenome,round(meanGs(i)),round(meanLongGenome(i)),round(meanBra(i)));
    end

    numGsR = zeros(1, sgenome);
    numBraR = zeros(1, sgenome);
    numSillyBraR = zeros(1, sgenome);
    numAnidR = zeros(1, sgenome);
    longGenomeR = zeros(1, sgenome);
    
    for j = 1 : sgenome
        numGsR(j) = sum(genomeRandom{i}{j} == 'G');
        numBraR(j) = sum(genomeRandom{i}{j} == '[');
        numSillyBraR(j) = sum(genome{i}{j}(1:end-1) == '[' & genome{i}{j}(2:end) == ']');
        numAnidR(j) = max(cumsum((genomeRandom{i}{j} == '[') + (-1 * (genomeRandom{i}{j} == ']'))));
        longGenomeR(j) = length(genomeRandom{i}{j});
    end

    meanGsR(i) = mean(numGsR);
    meanBraR(i) = mean(numBraR);
    meanSillyBraR(i) = mean(numSillyBra);
    meanAnidR(i) = mean(numAnidR);
    meanLongGenomeR(i) = mean(longGenomeR);
    
end

genome = [genome{:}];
genomeRandom = [genomeRandom{:}];

times = round(2000/(sgenome*total));

mode = struct('shadowing','transparentBranches','angle',22);

fLog = fopen(sprintf('%sRobustez_%s.log',finaldir, datestr(now, 30)), 'w');
fprintf(fLog, 'Evolved genomes: meanGs=%f meanBra=%f meanSillyBra=%f meanAnid=%f meanLongGenome=%f\r\n', mean(meanGs), mean(meanBra), mean(meanSillyBra), mean(meanAnid), mean(meanLongGenome));
fprintf(fLog, 'Random Genomes: meanGs=%f meanBra=%f meanSillyBra=%f meanAnid=%f meanLongGenome=%f\r\n\r\n', mean(meanGsR), mean(meanBraR), mean(meanSillyBraR), mean(meanAnidR), mean(meanLongGenomeR));

disrup = zeros(numOpFuncs, sgenome*total, times);
% evolution genome
for i = 1 : sgenome*total
    fprintf('Evolution genome %d/%d\n', i, sgenome*total);
    
    cad = genome{i};
    fprintf(fLog, 'Evolution genome %d/%d: %s\r\n', i, sgenome*total, cad);
    [raster, offset]=ls2(cad, mode);
    %raster = [];
    %offset = 0;
    for j = 1 : times
        for(k = 1:numOpFuncs)
            disrup(k,i,j) = calcDisrup(opFuncs{k}, cad, symbols, mode, raster, offset);
        end
    end
end
meanDisrup = mean(mean(disrup, 3), 2);
fprintf(fLog, 'meanDisrup = %s \r\n', mat2str(meanDisrup));

disrupR = zeros(numOpFuncs, sgenome*total, times);
% random genome
for i = 1 : sgenome*total
    fprintf('Random genome %d/%d\n', i, sgenome*total);
    cad = genomeRandom{i};
    fprintf(fLog, 'Random genome %d/%d: %s\r\n', i, sgenome*total, cad);
    [raster, offset]=ls2(cad, mode);
    for j = 1 : times
        for(k = 1:numOpFuncs)
            disrupR(k,i,j) = calcDisrup(opFuncs{k}, cad, symbols, mode, raster, offset);
        end
    end
end
meanDisrupR = mean(mean(disrupR, 3), 2);
fprintf(fLog, 'meanDisrupR = %s \r\n', mat2str(meanDisrupR));

fprintf(fLog, 'Done!\r\n');
fclose(fLog);

Y = [meanDisrup, meanDisrupR];

%Y = [meanDisrupGenomeA meanDisrupGenomeRA; meanDisrupGenomeE meanDisrupGenomeRE; meanDisrupGenomeI meanDisrupGenomeRI; meanDisrupGenomeDA meanDisrupGenomeRDA; meanDisrupGenomeDN meanDisrupGenomeRDN; meanDisrupGenomeDAE meanDisrupGenomeRDAE; meanDisrupGenomeDNE meanDisrupGenomeRDNE; meanDisrupGenomeDSE meanDisrupGenomeRDSE];
fig = figure;
bar(Y,'group');
xlabel('operators');
ylabel('mean disruption');
set(gca,'XTickLabel', plotXTickLabels);
legend('Evolved Individuals','Random Individuals');
saveas(fig,sprintf('%sRobustez_%s',finaldir, datestr(now, 30)),'epsc');
close(fig);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [disrup] = calcDisrup(funcOp, cad, symbols, mode, raster, offset)
cadOp = funcOp(cad, symbols);
if(~corchetesBalanceados(cadOp))
    error('regla mutada (%s) con corchetes no balanceados:\n original: %s\nmutante: %s', func2str(funcOp), cad, cadOp);
end
%disrup = 0;
%return;
[raster2, offset2] = ls2(cadOp, mode);
disrup = (1 - similitud2(raster, offset, raster2, offset2));
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cad t] = op_dupnivelG(cad1,symbols)

%selecciona la posición de inicio del segmento a duplicar
posCorAb = find(cad1=='[');
if(~isempty(posCorAb))
    pos = posCorAb(randint(1,1,[1,length(posCorAb)]));
    seg = extraer_segmento(cad1,pos);


    %busca en cad las posiciones que tienen el mismo nivel de
    %anidamiento que el segmento seleccionado y selecciona una de ellas
    %aleatoriamente
    nivelesAnid = cumsum((cad1 == '[') + (-1 * (cad1 == ']')));
    if(pos>1)
        nivelAnidPos = nivelesAnid(pos-1);
    else
        nivelAnidPos = 0;
    end
    posIgualNivel = find(nivelesAnid == nivelAnidPos);
    p = randint(1,1,[1,length(posIgualNivel)]);
    n = posIgualNivel(p) + 1;
    t = p/length(cad1);
    if t > 1
        disp('hola')
    end
    cad = [cad1(1:n-1) seg cad1(n:end)];
else
    cad = cad1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cad t] = op_duplSecuenciaEdit(cad1,symbols)

%selecciona aleatoriamente una de las posiciones susceptibles a
%duplicar
posCorAb = find(cad1=='[');
if(~isempty(posCorAb))
    pos = posCorAb(randint(1,1,[1,length(posCorAb)]));
    seg = extraer_segmento(cad1, pos);

    %inserta el segmento selecccionado en secuencia con el original
    n = pos + length(seg);
    cad = [cad1(1:n-1) seg cad1(n:end)];
    t = (n-1)/length(cad1);

    if length(seg) > 2
        r = randint(1,1,[1,3]);

        if r == 1
            %alteracion
            seg2 = symbols(randint(1,1,[1,length(symbols)]));
            pos = find((seg~='[') & (seg~=']'));
            if(~isempty(pos))
                m = randint(1,1,[1 length(pos)]); % seleccion de una posicion dentro del segmento duplicado
                cad(n+pos(m)-1) = seg2;
            end
        elseif r == 2
            %insercion
            seg2 = symbols(randint(1,1,[1,length(symbols)]));
            if seg2=='['
                seg2 = '[]';
            end
            m = randint(1,1,[n+1 n+length(seg)-2]); % seleccion de una posicion dentro del segmento duplicado
            cad = [cad(1:m-1) seg2 cad(m:end)];
        else
            %eliminacion
            pos = find((seg~='[') & (seg~=']'));
            if(~isempty(pos))
                m = randint(1,1,[1 length(pos)]); % seleccion de una posicion dentro del segmento duplicado
                cad(n+pos(m)-1) = [];
            end
        end
    end
else
    cad = cad1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [cad t] = op_duplNivelEdit(cad1,symbols)

%selecciona la posición de inicio del segmento a duplicar
posCorAb = find(cad1=='[');
if(~isempty(posCorAb))
    pos = posCorAb(randint(1,1,[1,length(posCorAb)]));
    seg = extraer_segmento(cad1,pos);
    nivelesAnid = cumsum((cad1 == '[') + (-1 * (cad1 == ']')));
    if(pos>1)
        nivelAnidPos = nivelesAnid(pos-1);
    else
        nivelAnidPos = 0;
    end
    posIgualNivel = find(nivelesAnid == nivelAnidPos);
    p = randint(1,1,[1,length(posIgualNivel)]);
    n = posIgualNivel(p) + 1;
    t = p/length(cad1);
    cad = [cad1(1:n-1) seg cad1(n:end)];

    if length(seg) > 2
        r = randint(1,1,[1,3]);

        if r == 1
            %alteracion
            seg2 = symbols(randint(1,1,[1,length(symbols)]));
            pos = find((seg~='[') & (seg~=']'));
            if(~isempty(pos))
                m = randint(1,1,[1 length(pos)]); % seleccion de una posicion dentro del segmento duplicado
                cad(n+pos(m)-1) = seg2;
            end
        elseif r == 2
            %insercion
            seg2 = symbols(randint(1,1,[1,length(symbols)]));
            if seg2=='['
                seg2 = '[]';
            end
            m = randint(1,1,[n+1 n+length(seg)-2]); % seleccion de una posicion dentro del segmento duplicado
            cad = [cad(1:m-1) seg2 cad(m:end)];
        else
            %eliminacion
            pos = find((seg~='[') & (seg~=']'));
            if(~isempty(pos))
                m = randint(1,1,[1 length(pos)]); % seleccion de una posicion dentro del segmento duplicado
                cad(n+pos(m)-1) = [];
            end
        end
    end
else
    cad = cad1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cad t] = op_duplAleatoriaEdit(cad1,symbols)

posCorAb = find(cad1=='[');
if(~isempty(posCorAb))
    pos = posCorAb(randint(1,1,[1,length(posCorAb)]));
    seg = extraer_segmento(cad1,pos);

    n = randint(1,1,[1,length(cad1)]);
    t = n/length(cad1);
    cad = [cad1(1:n-1) seg cad1(n:end)];

    if length(seg) > 2
        r = randint(1,1,[1,3]);

        if r == 1
            %alteracion
            seg2 = symbols(randint(1,1,[1,length(symbols)]));
            pos = find((seg~='[') & (seg~=']'));
            if(~isempty(pos))
                m = randint(1,1,[1 length(pos)]); % seleccion de una posicion dentro del segmento duplicado
                cad(n+pos(m)-1) = seg2;
            end
        elseif r == 2
            %insercion
            seg2 = symbols(randint(1,1,[1,length(symbols)]));
            if seg2=='['
                seg2 = '[]';
            end
            m = randint(1,1,[n+1 n+length(seg)-2]); % seleccion de una posicion dentro del segmento duplicado
            cad = [cad(1:m-1) seg2 cad(m:end)];
        else
            %eliminacion
            pos = find((seg~='[') & (seg~=']'));
            if(~isempty(pos))
                m = randint(1,1,[1 length(pos)]); % seleccion de una posicion dentro del segmento duplicado
                cad(n+pos(m)-1) = [];
            end
        end
    end
else
    cad = cad1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cad t] = op_dupaleatoriaG(cad1,symbols)

%selecciona aleatoriamente un segmento a duplicar

posCorAb = find(cad1=='[');
if(~isempty(posCorAb))
    pos = posCorAb(randint(1,1,[1,length(posCorAb)]));
    seg = extraer_segmento(cad1,pos);

    %inserta el segmento selecccionado en una posición de cad aleatoria
    n = randint(1,1,[1,length(cad1)]);
    t = n/length(cad1);
    cad = [cad1(1:n-1) seg cad1(n:end)];
else
    cad = cad1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cad t] = op_alteracionG(cad1,symbols)

% symbols = ['G'];%['G','+','-'];
seg = symbols(randint(1,1,[1,length(symbols)]));

% selecciona aleatoriamente uno de las posiciones de la cadena cad
% sin incluir corchetes
pos = find(cad1 < '[');

n = pos(randint(1,1,[1 length(pos)])); % seleccion de una posicion

cad = cad1;
cad(n) = seg;
t = n/length(cad1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cad t] = op_insercionG(cad1,symbols)

%selecciona una posición de cad aleatoria donde se realizará la
%inserción del símbolo también seleccionado aleatoriamente

% symbols                              = ['+','-']; %['G','+','-','['];

seg = symbols(randint(1,1,[1,length(symbols)]));

if seg=='['
    seg = '[]';
end
n = randint(1,1,[1,length(cad1)]);
cad = [cad1(1:n-1) seg cad1(n:end)];
t = n/length(cad1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cad t]= op_eliminacionG(cad1,symbols)
% symbols                              = ['G']; %['G','+','-','['];
perm                                 = randperm(numel(symbols));
operAplicada                         = false;

for z=1:numel(symbols)
    seg = symbols(perm(z));
    if(seg=='[')
        % implica eliminación de un corchete vacío []
        k                              = strfind(cad1, '[]');
        if ~isempty(k)
            n                            = k(randint(1,1,[1, length(k)]));
            cad                          = cad1;
            cad(n:(n+1))                 = [];
            operAplicada                 = true;
        end
    else
        k = find(cad1 == seg);
        if(~isempty(k))
            n                              = k(randint(1,1,[1, length(k)]));
            cad                            = cad1;
            cad(n)                         = [];
            operAplicada                   = true;
        end
    end
    if operAplicada
        break;
    end
end
t = n/length(cad1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sim = similitud2(raster1,offset1,raster2,offset2)%, show)

if isempty([raster1{:}]) && isempty([raster2{:}])
    sim = 1;
elseif (isempty([raster1{:}]) && ~isempty([raster2{:}])) || (~isempty([raster1{:}]) && isempty([raster2{:}]))
    sim = 0;
else
    pixel1 = [raster1{:}];
    pixel2 = [raster2{:}];
    
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
    
    sim = size(intersect(pixel1,pixel2,'rows'),1)/size(union(pixel1,pixel2,'rows'),1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function genomeRandom=generateGenomeRTipo3(sgenome,meanGs,meanLongGenome,meanBra)

genomeRandom = cell(1, sgenome);

meanUGs = meanGs/meanLongGenome;
meanUBra = meanBra/meanLongGenome;

for i = 1 : sgenome
    sizeGenome = randint(1,1,[meanLongGenome-10,meanLongGenome+10]);
    genome = repmat('x', 1, sizeGenome);
    
    numG = round(meanUGs * sizeGenome);
    %numG = randint(1,1,[meanGs-10, min(meanGs+10, sizeGenome-2)]); % -2 in order to have at least a pair of []
    permutaciones = randperm(sizeGenome);
    posG = permutaciones(1:numG);
    
    numBra = round(2 * meanUBra * sizeGenome);
    %numBra = randint(1,1,[2 sizeGenome-numG]);
    posBra = permutaciones(numG+1:min(sizeGenome,numG+numBra));
    if(mod(length(posBra),2)~=0)
        posBra(end) = [];
        numBra = numBra - 1;
    end
    posBra = sort(posBra);

    posSignos = permutaciones(numG+numBra+1:sizeGenome);

    genome(posG) = 'G';
    genome(posBra(1:2:end-1)) = '['; % genoma Tipo 3
    genome(posBra(2:2:end)) = ']';
    symbols = ['+','-'];
    genome(posSignos) = symbols(randint(1,length(posSignos),[1,length(symbols)]));


    noPlus = isempty(find(genome == '+', 1));
    noMinus = isempty(find(genome == '-', 1));
    if noPlus && noMinus
        genome(posBra(1)) = '+';
        genome(posBra(2)) = 'G';
    end
    genomeRandom{i} = genome;
    
    %fprintf('%s\n', genome);
    if(~corchetesBalanceados(genome))
        error('regla con corchetes no balanceados: %s', genome);
    end
    if(isempty(find(genome=='[', 1)))
        error('genoma sin corchetes: %s', genome);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function genomeRandom=generateGenomeRTipo2(sgenome,meanGs,meanLongGenome,meanBra, meanAnid)

genomeRandom = cell(1, sgenome);

meanUGs = meanGs/meanLongGenome;
meanUBra = meanBra/meanLongGenome;

for i = 1 : sgenome
    sizeGenome = randint(1,1,[meanLongGenome-10,meanLongGenome+10]);
    genome = repmat('x', 1, sizeGenome);
    
    numG = round(meanUGs * sizeGenome);
    %numG = randint(1,1,[meanGs-10, min(meanGs+10, sizeGenome-2)]); % -2 in order to have at least a pair of []
    permutaciones = randperm(sizeGenome);
    posG = permutaciones(1:numG);
    
    numBra = round(2 * meanUBra * sizeGenome);
    %numBra = randint(1,1,[2 sizeGenome-numG]);
    posBra = permutaciones(numG+1:min(sizeGenome,numG+numBra));
    if(mod(length(posBra),2)~=0)
        posBra(end) = [];
        numBra = numBra - 1;
    end
    posBra = sort(posBra);
    %posBra = reshape(sort(reshape(posBra, 2, [])), 1, []);
    genome(posBra(numBra/2)) = '[';
    genome(posBra(numBra/2 + 1)) = ']';
    rp = randperm(numBra/2 - 1);
    leftAnidIndex = rp(1:meanAnid-1);
    rp = randperm(numBra/2 - 1) + numBra/2 + 1;
    rightAnidIndex = rp(1:meanAnid-1);
	genome(posBra(leftAnidIndex)) = '[';	
    genome(posBra(rightAnidIndex)) = ']';
    posBra([numBra/2, numBra/2 + 1, leftAnidIndex, rightAnidIndex]) = [];
    for(j=1:2:length(posBra))
        genome(posBra(j)) = '[';
        genome(posBra(j+1)) = ']';
    end
    
    posSignos = permutaciones(numG+numBra+1:sizeGenome);

    genome(posG) = 'G';
    %genome(posBra(1:2:end-1)) = '['; % genoma Tipo 3
    %genome(posBra(2:2:end)) = ']';
    symbols = ['+','-'];
    genome(posSignos) = symbols(randint(1,length(posSignos),[1,length(symbols)]));


    noPlus = isempty(find(genome == '+', 1));
    noMinus = isempty(find(genome == '-', 1));
    if noPlus && noMinus
        genome(posBra(1)) = '+';
        genome(posBra(2)) = 'G';
    end
    genomeRandom{i} = genome;
    %fprintf('%s\n', genome);
    if(~corchetesBalanceados(genome))
        error('regla con corchetes no balanceados: %s', genome);
    end
    if(isempty(find(genome=='[', 1)))
        error('genoma sin corchetes: %s', genome);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function seg = extraer_segmento(cad,pos)
% Extrae toda una rama del árbol empezando en pos (cad(pos) = '[')

subCad = cad(pos:end);
posCor = (subCad == '[') + (-1 * (subCad == ']'));

finNivel = find(cumsum(posCor)==0, 1);
seg = subCad(1:finNivel);
