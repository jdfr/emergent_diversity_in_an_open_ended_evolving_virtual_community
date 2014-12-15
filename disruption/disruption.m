function disrup = disruption(numoperator,cad1,cad2)

mode = 'canonical';
[raster1, offset1, counts1, dimension1, isLeaf1, maxiy1]=ls2(cad1, mode);

if numoperator == 1
    for i = 1 : 100
        cad{i} = op_eliminacionG(cad1);
    end
elseif numoperator == 2
    for i = 1 : 100
        cad{i} = op_insercionG(cad1);
    end
elseif numoperator == 3
    for i = 1 : 100
        cad{i} = op_alteracionG(cad1);
    end
elseif numoperator == 4
    for i = 1 : 100
        cad{i} =op_dupaleatoriaG(cad1);
    end
elseif numoperator == 5
    for i = 1 : 100
        cad{i} =op_dupnivelG(cad1);
    end
elseif numoperator == 6
    for i = 1 : 100
        cad{i} =op_dupsecuenciaG(cad1);
    end
elseif numoperator == 7
    for i = 1 : 100
        cad{i} =op_transferenciaG(cad1,cad2);
    end
end


disrup = 0;
for i = 1 : 100
    [raster2, offset2, counts2, dimension2, isLeaf2, maxiy2]=ls2(cad{i}, mode);
    disrup = disrup + (1 - similitud2(raster1,offset1,raster2,offset2));
end

disrup = disrup/100;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    sim = size(intersect(pixel1,pixel2,'rows'),1)/size(union(pixel1,pixel2,'rows'),1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cad = op_eliminacionG(cad1)
symbols                              = ['G','+','-','['];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cad = op_insercionG(cad1)

%selecciona una posición de cad aleatoria donde se realizará la
%inserción del símbolo también seleccionado aleatoriamente

symbols                              = ['G','+','-','['];

seg = symbols(randint(1,1,[1,length(symbols)]));

if seg=='['
    seg = '[]';
end
n = randint(1,1,[1,length(cad1)]);
cad = [cad1(1:n-1) seg cad1(n:end)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cad = op_alteracionG(cad1)

symbols = ['G','+','-'];
seg = symbols(randint(1,1,[1,length(symbols)]));

% selecciona aleatoriamente uno de las posiciones de la cadena cad
% sin incluir corchetes
pos = find(cad1 < '[');

n = pos(randint(1,1,[1 length(pos)])); % seleccion de una posicion

cad = cad1;
cad(n) = seg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cad = op_dupaleatoriaG(cad1)

%selecciona aleatoriamente un segmento a duplicar

posCorAb = find(cad1=='[');
pos = posCorAb(randint(1,1,[1,length(posCorAb)]));
seg = extraer_segmento(cad1,pos);

%inserta el segmento selecccionado en una posición de cad aleatoria
n = randint(1,1,[1,length(cad1)]);
cad = [cad1(1:n-1) seg cad1(n:end)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cad = op_dupnivelG(cad1)

%selecciona la posición de inicio del segmento a duplicar
posCorAb = find(cad1=='[');
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
cad = [cad1(1:n-1) seg cad1(n:end)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cad = op_dupsecuenciaG(cad1)

%selecciona aleatoriamente una de las posiciones susceptibles a
%duplicar
posCorAb = find(cad1=='[');
pos = posCorAb(randint(1,1,[1,length(posCorAb)]));
seg = extraer_segmento(cad1, pos);

%inserta el segmento selecccionado en secuencia con el original
n = pos + length(seg);
cad = [cad1(1:n-1) seg cad1(n:end)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cad = op_transferenciaG(cad1,cad2)


%--------------------------------------------------------------------------------------------------------------
% Operador de transferencia  21/01/2009
% Como generalización del operador de Cruzamiento, la Transferencia
% consistirá en seleccionar dos individuos y hacer un descendiente
% igual a uno de ellos, con cualquier segmento (sección del genoma entre corchetes,
% no importa lo que contenga) del otro copiado en el genoma del primero.
% Inicialmente no vamos a considerar variantes, el lugar de inserción será aleatorio.
% p.e. G-G+[GG+[+G]] y GG+[-G[G]] son seleccionados y se hace uno nuevo G-G+[G[-G[G]]G+[+G]],
% con el primero insertando [-G[G]].
%
% Simularás este nuevo operador sin combinar con duplicación (poner a 0 sus probabilidades en main.m).
% Vamos a ver algunos resultados y decidimos por dónde tirar.

% extraemos un seg aleatorio de cad2
posCorAb = find(cad2=='[');
pos = posCorAb(randint(1,1,[1,length(posCorAb)]));

seg = extraer_segmento(cad2,pos);


%inserta el segmento selecccionado de cad2 en una posición de cad1 aleatoria
n = randint(1,1,[1,length(cad1)]);
cad = [cad1(1:n-1) seg cad1(n:end)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function seg = extraer_segmento(cad,pos)
% Extrae toda una rama del árbol empezando en pos (cad(pos) = '[')

subCad = cad(pos:end);
posCor = (subCad == '[') + (-1 * (subCad == ']'));

finNivel = find(cumsum(posCor)==0, 1);
seg = subCad(1:finNivel);
