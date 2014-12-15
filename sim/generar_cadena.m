function cad = generar_cadena(num_simb, num_anidam)

%inicialización de variables
symbols = 'G+-[]'; %'GF+-[]';
ls = length(symbols);
cad =repmat(' ', 1, num_simb);
idx_end = 0;
cont_ab = 0;

num_simb = ceil(rand*num_simb);
num_anidam = ceil(rand*num_anidam);

nsimbNC = 0;

%si quedan posiciones suficientes para balancear los parentesis entramos a
%generar un nuevo símbolo aleatorio. En otro caso, completamos el resto
%de cadena con cierre de corchetes para balancear la cadena

while (nsimbNC < num_simb)

   %p = randint(1,1,[1,ls]);          %genera un numero aleatorio entre 1 y 6
   p = floor(rand*ls)+1; %far faster
   e = symbols(p);                  %devuelve el elemento del vector symbols que ocupa la posicion definida por p

  switch e
       case '['
           if (cont_ab < num_anidam)
               cont_ab = cont_ab + 1;
               idx_end = idx_end + 1;
               cad(idx_end) = e;
           end

       case ']'
           % Caso: Si existen corchetes abiertos en la cadena, el último elemento de ésta no es un abre-corchete ni signos '+' o '-'
           if ((idx_end>0) && (cont_ab > 0)) && all('[+-'~=cad(idx_end))
              idx_end = idx_end + 1;
              cad(idx_end) = e;
              cont_ab = cont_ab - 1;
           end

       % Caso ('+' & '-'): Si se concatenan varios '+' y '-' seguidos se anulan. Por
       % lo que sólo concatenaré el símbolo si el último elemento de la
       % cadena no es su inverso.
       case '+'
          % Si estamos en posiciones cercanas al valor max, debemos
          % comprobar, para concatenar el signo '+', que podemos además
          % concatenar un símbolo 'F' o 'G'
          if ((idx_end>0) && ((nsimbNC + 1) < num_simb)) && (cad(idx_end)~='-') %fixed: operator + takes precedence over <
             idx_end = idx_end + 1;
             cad(idx_end) = e;
             nsimbNC = nsimbNC + 1;
          end
       case '-'
          if ((idx_end>0) && ((nsimbNC + 1) < num_simb)) && (cad(idx_end)~='+') %fixed: operator + takes precedence over <
              idx_end = idx_end + 1;
              cad(idx_end) = e;
              nsimbNC = nsimbNC + 1;
          end

       otherwise
         idx_end = idx_end + 1;
         cad(idx_end) = e;
         nsimbNC = nsimbNC + 1;
  end %end switch

end %end while

%if longitud(cad) < num_simb
  cad = [cad(1:idx_end), repmat(']', 1, cont_ab)];
%end

%cad = ['[',cad,']'];

%cad = ['G->',cad];         eliminado temporalmente



% %-------------------------------------------------------------------------
% function num = longitud(cad,idx_end)
% %función que calcula el nº de símbolos sin contar los corchetes
% num = sum((cad(1:idx_end)~='[')&(cad(1:idx_end)~=']'));
