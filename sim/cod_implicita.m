function [cad,seg,n,operatorChar,disruption] = cod_implicita(oper,cad1,cad2,maxLongs,params)
%La función recibe como parámetros los siguientes elementos:
%   oper: String que determina que operador genético se ejecutará.
%   cad: String que contiene el individuo sobre el cual se aplicarán los
%   operadores genéticos.
%   seg: String que representa el segmento a insertar en cad según el
%   operador genético que se aplique.
%   pos: Int. Representa la posición a partir de la cual se aplicará el
%   operador genético.

%Esta función contiene los operadores genéticos descritos a continuación:
%-- Alteracion (A): se cambia un único símbolo (G,F,+ o -) por otro.  Este
%operador no afecta a los corchetes.

%-- Duplicación aleatoria(DA): igual que la DN, pero se puede insertar en
%cualquier sitio.

%-- Duplicación por nivel(DN): igual que la DS, pero se puede insertar
%en cualquier sitio dentro del mismo nivel de anidamiento, es decir,
%que debe tener por delante los mismos corchetes abiertos.
%Por ejemplo, ...[...[...[G+F]]]...]... => ...[...[...[G+F]]]...]...[[[...[G+F]]...]]

%-- Duplicación secuencia (DS): consiste en copiar un segmento entre
%corchetes, tal y como está y a continuación. Por ejemplo, ...][G+F][... => ...][G+F][G+F][...
%Suponemos que las cadenas que recibe como parámetro esta función son
%cadenas correctamente balanceadas

%-- Eliminación (E): se elimina un único símbolo (G,F,+ o -), o bien un par
%de corchetes vacíos ([])

%-- Inserción (I): se inserta un sólo símbolo aleatoriamente (G,F,+ o -), o
%bien un par de corchetes vacíos ([]).

%Estos operadores están sujetos a una probabilidad de ejecución relativaa una
%probabilidad de mutación (pM) y otra relativa a cada operador (pDS, etc).

%--------------------------------------------------------------------------

operAplicada = false;
operatorChars = ' ARLCDIT'; % == {ELITISM (NO CHAR USED)} ALTERATION, RANDOM, LEVELED, CONTIGUOUS, DELETION; INSERTION, TRANSFER
operatorChar = operatorChars(oper);
cad = '';
doNotAllowPCD = (~params.allowPrefixChangeDir) && (~firstIsChangeDir(cad1));
switch oper
  case 2 %op_alteracion
    %selecciona aleatoriamente un símbolo del conjunto de símbolos
    %posibles
    symbols = ['G','+','-'];

    % selecciona aleatoriamente uno de las posiciones de la cadena cad
    % sin incluir corchetes
    pos = find((cad1~='[')&(cad1~=']'));
    if(~isempty(pos))
        continuar = true; contador = 0;
        while continuar
          operAplicada = false;

          n = pos(randint(1,1,[1 length(pos)])); % seleccion de una posicion

          cad = cad1;
          %selecciona un símbolo distinto del que ya hay
          seg = symbols(symbols~=cad(n));
          seg = seg(ceil(rand*numel(seg)));
          cad(n) = seg;

          contador = contador+1;
          continuar = ((~isempty(maxLongs.maxGs)) && (sum(cad=='G')>maxLongs.maxGs));
          if contador>100; break; end
          if continuar || all(cad~='G'); continue; end
          if doNotAllowPCD && firstIsChangeDir(cad); continue; end

          operAplicada = true;
        end
    end
        
  case 3    %op_dupaleatoria
    %selecciona aleatoriamente un segmento a duplicar
    posCorAb = find(cad1=='[');
    if(~isempty(posCorAb))
        continuar = true; contador = 0;
        while continuar
          operAplicada = false;
          pos = posCorAb(randint(1,1,[1,length(posCorAb)]));
          seg = extraer_segmento(cad1,pos);

          [continuar contador] = condicion(cad1, seg, maxLongs, contador);
          if contador>100; break; end
          if continuar; continue; end

          %inserta el segmento selecccionado en una posición de cad aleatoria
          n = randint(1,1,[1,length(cad1)]);
          cad = [cad1(1:n-1) seg cad1(n:end)];
          operAplicada = true;
        end
    end

  case 4    %op_dupnivel
    %selecciona la posición de inicio del segmento a duplicar
    posCorAb = find(cad1=='[');
    if(~isempty(posCorAb))
        continuar = true; contador = 0;
        while continuar
          operAplicada = false;
          pos = posCorAb(randint(1,1,[1,length(posCorAb)]));
          seg = extraer_segmento(cad1,pos);

          [continuar contador] = condicion(cad1, seg, maxLongs, contador);
          if contador>100; break; end
          if continuar; continue; end

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
          operAplicada = true;
        end
    end

  case 5    %op_dupsecuencia
    %selecciona aleatoriamente una de las posiciones susceptibles a
    %duplicar
    posCorAb = find(cad1=='[');
    if(~isempty(posCorAb))
        continuar = true; contador = 0;
        while continuar
          operAplicada = false;
          pos = posCorAb(randint(1,1,[1,length(posCorAb)]));
          seg = extraer_segmento(cad1, pos);

          [continuar contador] = condicion(cad1, seg, maxLongs, contador);
          if contador>100; break; end
          if continuar; continue; end

          %inserta el segmento selecccionado en secuencia con el original
          n = pos + length(seg);
          cad = [cad1(1:n-1) seg cad1(n:end)];
          operAplicada = true;
        end
    end

  case 6    %op_eliminacion
    if (numel(cad1)>1) && ((numel(cad1)>2) || (cad1(1)~='['))
      % selecciona aleatoriamente el símbolo a eliminar de la cadena cad firstIsChangeDir(cad)
      [cad,seg,n,operAplicada] = op_eliminacion(cad1, params.op_eliminacion_extended, params.selectSymbolForMutation, doNotAllowPCD);
    end

  case 7    %op_insercion
    if iscell(params.mode_op_insercion)
      [cad,seg,n, operAplicada] = op_insercion(params, params.mode_op_insercion, cad1,maxLongs,doNotAllowPCD);
    else
      switch params.mode_op_insercion(1)
        case 'n' %'normal'
          [cad,seg,n, operAplicada] = op_insercion(params, {'G','+','-','[]'}, cad1,maxLongs,doNotAllowPCD);
        case 'o' %'one_level'
          if rand>0.25
            %classic insertion
            [cad,seg,n, operAplicada] = op_insercion(params, {'G','+','-'}, cad1,maxLongs,doNotAllowPCD);      
          else
            %insert empty brackets
            seg = '[]';
            if ~condicion(cad1, seg, maxLongs, 0) %if we are allowed to do it:
              %find suitable positions to put the empty brackets
              positions = find([0 cumsum((cad1=='[')-(cad1==']'))]==0);
              n         = positions(ceil(rand*numel(positions)));
              cad       = [cad1(1:(n-1)) seg cad1(n:end)];
            end
            operAplicada = true;
          end
        otherwise
          error('Bad parameter mode_op_insercion!!!');
      end
    end
      
  case 8 %op_transferencia
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
    if(~isempty(posCorAb))
        continuar = true; contador = 0;
        while continuar
          operAplicada = false;
          pos = posCorAb(randint(1,1,[1,length(posCorAb)]));

          seg = extraer_segmento(cad2,pos);

          [continuar contador] = condicion(cad1, seg, maxLongs, contador);
          if contador>100; break; end
          if continuar; continue; end

          %inserta el segmento selecccionado de cad2 en una posición de cad1 aleatoria
          if(~isempty(cad1))
              n = randint(1,1,[1,length(cad1)]);
              cad = [cad1(1:n-1) seg cad1(n:end)];
          else
              n = 1;
              cad = seg;
          end
          operAplicada = true;
        end
    end

  otherwise
    error('Error: cod_implicita: operador erróneo (%d)', oper);
end

if(~operAplicada)
  % Si no se consigue aplicar el operador seleccionado, se aplica el
  % operador inserción, que siempre es posible aplicar.
  if oper~=7 %make sure that we don't get trapped into an infinite recursion in limit cases
    [cad,seg,n,operatorChar] = cod_implicita(7,cad1, '', maxLongs, params);
  else
    cad=cad1;
    seg='=';
    n=0;
  end
end

if nargout>=5
  disruption =1 - similitud(cad1,cad);
end


%--------------------------------------------------------------------------
function seg = extraer_segmento(cad,pos)
% Extrae toda una rama del árbol empezando en pos (cad(pos) = '[')

subCad = cad(pos:end);
posCor = (subCad == '[') + (-1 * (subCad == ']'));

finNivel = find(cumsum(posCor)==0, 1);
seg = subCad(1:finNivel);


%--------------------------------------------------------------------------
function [continuar contador] = condicion(cad1, seg, maxLongs, contador)
contador = contador+1;
continuar = ...%(contador<100) && ...
  ( ((~isempty(maxLongs.maxInsert)) && (numel(seg)>maxLongs.maxInsert)) || ...
    ((~isempty(maxLongs.maxGs)) && ((sum(seg=='G')+sum(cad1=='G'))>maxLongs.maxGs)) || ...
    ((~isempty(maxLongs.maxString)) && ((numel(seg)+numel(cad1))>maxLongs.maxString)) ...
  );
%continuar=continuar;


%--------------------------------------------------------------------------
function [cad,seg,n,operAplicada] = op_insercion(params, symbols,cad1,maxLongs,doNotAllowPCD)
continuar = true; contador = 0;
while continuar
  operAplicada = false;
  cad = cad1;
  n   = 0;
  seg = symbols{randint(1,1,[1,length(symbols)])};

  [continuar contador] = condicion(cad1, seg, maxLongs, contador);
  if contador>100; break; end
  if continuar; continue; end

  if params.insertOnlyNested
    %insert only inside brackets
    positions = find(cumsum((cad1=='[')-(cad1==']'))>0);
    if ~isempty(positions)
      n   = positions(ceil(rand*numel(positions)));
      cad = [cad1(1:n) seg cad1(n+1:end)];
      n   = n-1;
      operAplicada = true;
    end
  else
    %selecciona una posición de cad aleatoria donde se realizará la
    %inserción del símbolo también seleccionado aleatoriamente
    if(~isempty(cad1))
        n = randint(1,1,[1,length(cad1)]);
        cad = [cad1(1:n-1) seg cad1(n:end)];
    else
        n = 1;
        cad = seg;
    end
    operAplicada = true;
  end
  operAplicada = operAplicada && ((~doNotAllowPCD) || (~firstIsChangeDir(cad)));
end

%--------------------------------------------------------------------------
function [cad,seg,n,operAplicada]    = op_eliminacion(cad1, deleteNesting, selectSymbolForMutation, doNotAllowPCD)
operAplicada                         = false;
switch selectSymbolForMutation
  case 'bySymbol'
    symbols                          = ['G','+','-','['];
    bySymbol                         = true;
  case 'byPos'
    symbolPos                        = find(cad1~=']');
    symbols                          = cad1(symbolPos);
    bySymbol                         = false;
  otherwise
    error('Unknown value <%s> for option selectSymbolForMutation!!!!!', any2str(params.selectSymbolForMutation));
end
perm                             = randperm(numel(symbols));
for z=1:numel(symbols)
  seg = symbols(perm(z));
  if ~bySymbol
    n                            = symbolPos(perm(z));
  end
  if seg=='['   % implica eliminación de un nivel de anidamiento
    if deleteNesting
      if bySymbol
        k                        = find(cad1 == '[');
        n                        = k(randint(1,1,[1, length(k)]));
      else
        k                        = 1;
      end
      ndel                       = find(cumsum((cad1(n:end)=='[')-(cad1(n:end)==']'))==0, 1);
    else        % implica eliminación de un corchete vacío []
      ndel                       = 2;
      if bySymbol
        k                        = strfind(cad1, '[]');
        n                        = k(randint(1,1,[1, length(k)]));
      else
        k                        = ndel(cad1(n+1)==']');
      end
    end
  else
    ndel                         = 1;
    if bySymbol
      k                          = find(cad1 == seg);
      n                          = k(randint(1,1,[1, length(k)]));
    else
      k                          = 1;
    end
  end
  if ~isempty(k)
    cad                          = cad1;
    cad(n:(n+ndel-1))            = [];
    operAplicada                 = any(cad=='G') && ((~doNotAllowPCD) || (~firstIsChangeDir(cad)));
  end
  if operAplicada
    break;
  end
end

if ~operAplicada
  n                                  = 0;
  cad                                = '';
end

%--------------------------------------------------------------------------
function firstIsChgDir = firstIsChangeDir(cad)

cad0depth = cad(cumsum((cad=='[')-(cad==']'))==0);

pG  = find(cad0depth=='G', 1, 'first');
pCD = find((cad0depth=='+') | (cad0depth=='-'), 1, 'first');

clear cad0depth;

if isempty(pG)
  firstIsChgDir = ~isempty(pCD);
else
  firstIsChgDir = (~isempty(pCD)) && (pG>pCD);
end

