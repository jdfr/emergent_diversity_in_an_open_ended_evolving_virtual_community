function [meanDisrup plotXTickLabels disrup] = calculateDisruption(genome, times, seeds, dst, doparallel)

opFuncs = {@op_alteracionG, ...
           @op_eliminacionG, ...
           @op_insercionG, ...
           @op_dupaleatoriaG, ...
           @op_dupnivelG, ...
           @op_duplAleatoriaEdit, ...
           @op_duplNivelEdit, ...
           @op_duplSecuenciaEdit};

plotXTickLabels = 'A|D|I|R|L|R+E|L+E|T+E';

sgenome = numel(genome);

numOpFuncs = numel(opFuncs);

fLog = 1;

% evolution genome
if not(isempty(doparallel))
  jobArgs = doparallel.jobArgs;

  if isempty(doparallel.numJobs)
    nj = numel(genome);
  else
    nj = doparallel.numJobs;
  end
  [ranges rangesizes] = calculateRanges(numel(genome), nj);
  genomes = mat2cell(genome, rangesizes, 1);
  seedss = mat2cell(seeds, rangesizes, 1);
  timess = repmat({times}, numel(rangesizes), 1);
  opFuncss = repmat({opFuncs}, numel(rangesizes), 1);
  dsts = repmat({dst}, numel(rangesizes), 1);
  fprintf('GOING PARALLEL...\n');
  disrupss = my_dfeval(@doitseveralTimes, seedss, genomes, timess, opFuncss, dsts, jobArgs{:});
  fprintf('RECEIVING PARALLEL...\n');
  %save('disrupss.mat', 'disrupss');
  disrup = zeros(numOpFuncs, sgenome, times);
  for k=1:numel(disrupss)
    %fprintf('MIRA:  %s: %s, %s\n', mat2str(size(disrup)), mat2str(size(disrupss{k})), mat2str(ranges(k,:)));
    disrup(:,ranges(k,1):ranges(k,2),:) = disrupss{k};
  end
else
  disrup = doitseveralTimes(seeds, genome, times, opFuncs, dst);
end
meanDisrup = mean(mean(disrup, 3), 2);
fprintf(fLog, 'meanDisrup = %s \r\n', mat2str(meanDisrup));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function disrups = doitseveralTimes(seeds, cads, times, opFuncs, dst)
numOpFuncs = numel(opFuncs);
disrups =  zeros(numOpFuncs, numel(cads), times);
sgenome = numel(cads);
for k=1:numel(cads)
  fprintf('Evolution genome %d/%d\n', k, sgenome);
  disrups(:,k,:) = doit(seeds(k), times, opFuncs, cads{k}, dst);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function disrup = doit(seed, times, opFuncs, cad, dst)
mode = struct('shadowing', 'canonical', 'angle', 22, 'justPaint', true);
symbols = 'G+-';
numOpFuncs = numel(opFuncs);
    rand('twister', seed);
    disrup = zeros(numOpFuncs, 1, times);
    [raster, offset, counts]=ls2(cad, mode); %counts  = [nramas nleafs nPixelsNeg];
    minicad = sanitizeString(cad);

%     bc = fprintf(' ');
    for j = 1 : times
            fprintf('%d,', mod(j, 10));
        for k = 1:numOpFuncs
%             fprintf(repmat('\b', 1, bc));
%             bc = fprintf('j=%d/%d, k=%d/%d', j, times, k, numOpFuncs);
%             fprintf('j=%d/%d, k=%d/%d\n', j, times, k, numOpFuncs);
            disrup(k,1,j) = calcDisrup(opFuncs{k}, cad, minicad, symbols, mode, raster, offset, counts, dst);
        end
    end
    fprintf('\n');
%             fprintf(repmat('\b', 1, bc));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [disrup] = calcDisrup(funcOp, cad, minicad, symbols, mode, raster, offset, counts, dst)
cadOp = funcOp(cad, symbols);
minicadop = sanitizeString(cadOp);
if isempty(minicadop)||not(any(upper(minicadop)=='G'))
  disrup=1;
  return;
end
if(~corchetesBalanceados(cadOp))
    error('regla mutada (%s) con corchetes no balanceados:\n original: %s\nmutante: %s', func2str(funcOp), cad, cadOp);
end
if strcmp(minicad, minicadop)
  disrup = 0;
else
  %disrup = 0;
  %return;
  switch dst
    case 'MP'
      [raster2, offset2] = ls2(minicadop, mode);
      disrup = (1 - similitud2(raster, offset, raster2, offset2));
    case 'LOGBN'
      [raster2, offset2, counts2] = ls2(minicadop, mode); %counts  = [nramas nleafs nPixelsNeg];
      disrup = abs(log10(counts(1))-log10(counts2(1))); %abs(diff(log10(nramas),1,2))
  end
end        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cad t] = op_dupnivelG(cad1,symbols) %#ok<INUSD>

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

function [cad t] = op_dupaleatoriaG(cad1,symbols) %#ok<INUSD>

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function seg = extraer_segmento(cad,pos)
% Extrae toda una rama del árbol empezando en pos (cad(pos) = '[')

subCad = cad(pos:end);
posCor = (subCad == '[') + (-1 * (subCad == ']'));

finNivel = find(cumsum(posCor)==0, 1);
seg = subCad(1:finNivel);

