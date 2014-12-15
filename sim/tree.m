function [ramas,nramas,alto,ancho,leaf,maxiy]=tree(axiom, mode)
% axiom      cadena a representar
% pintar     representar graficamente si pintar=true
% tiempo     tiempo (s) para visualizar la grafica
%
% alto       alto del bounding box
% ancho      ancho del bounding box
% ramas      numero de ramas

chanceFun = @()true;
useTransparency  = false;
switch mode
    case 'canonical'
    case 'probabilistic'
        probability      = 0.8;
        chanceFun        = @()rand<=probability;
    case 'transparentBranches'
        useTransparency  = true;
    otherwise
        error('Mode %s not understood!!!', anytostr(mode));
end

minlong = 0.00001;    % alto o ancho minimo del bounding box, para evitar divisiones por cero

% angle: '+' = rotate anticlockwise
%        '-' = rotate clockwise
alpha = -5; % degrees (angle)

% lengths of the lines F and G
% length_F = 1;
length_G = 10;

% initialization
xT = 0;
yT = 0;
aT = 0;
da = alpha/180*pi; % convert deg to rad
stkPtr = 1;
stack = zeros(50,3);

colorG = double('b'); % pinta las ramas en azul
% colorF = double('r');


%nodos = [xT yT];
ramas = [];
g = 1;
for i = 1:length(axiom)
    cmdT = axiom(i);
    switch cmdT
        % It is possible to add multiple cases here
        % in order to expand the program.
        case 'G'
            if g == 1 % the first G
                newxT = xT + length_G*cos(aT);
                newyT = yT + length_G*sin(aT);
                %nodos = cuenta_filas(nodos,[newxT newyT], 2);
                %ramas = cuenta_filas(ramas,[yT newyT xT newxT colorG], 4);
                ramas = sortedInsertNoRepeat(ramas,[yT newyT xT newxT double('b')], 4); % origen del grafo (xT,yT)
                xT = newxT;
                yT = newyT;
                g = 2;
            else
                newxT = xT + length_G*cos(aT);
                newyT = yT + length_G*sin(aT);
                %nodos = cuenta_filas(nodos,[newxT newyT], 2);
                %ramas = cuenta_filas(ramas,[yT newyT xT newxT colorG], 4);
                ramas = sortedInsertNoRepeat(ramas,[yT newyT xT newxT colorG], 4);
                xT = newxT;
                yT = newyT;
            end

        case '+' % rotate anticlockwise
            aT = aT + da;
        case '-' % rotate clockwise
            aT = aT - da;
        case '[' % save current position
            stack(stkPtr,:) = [xT yT aT];
            stkPtr = stkPtr +1 ;
        case ']' % return to former position (last save)
            stkPtr = stkPtr -1 ;
            xT = stack(stkPtr,1);
            yT = stack(stkPtr,2);
            aT = stack(stkPtr,3);
        case 'F'
            newxT = xT + length_F*cos(aT);
            newyT = yT + length_F*sin(aT);
            %nodos = cuenta_filas(nodos,[newxT newyT], 2);
            %ramas = cuenta_filas(ramas,[yT newyT xT newxT colorF], 4);
            ramas = sortedInsertNoRepeat(ramas,[yT newyT xT newxT colorF], 4);
            xT = newxT;
            yT = newyT;
        otherwise
            disp('error');
            return
    end

end

if ~isempty(ramas)
    r = ramas(:,1:4);
    r = round(r);
    minimax = min(min(r(:,1:2)));
    x = 1 - minimax;
    r(:,1:2) = r(:,1:2)+x;
    ramas = [r ramas(:,5)]; % me situa el arbol en x >= 1 (discretizado)
end

if nargout<=1
    return
end

nramas = size(ramas,1);      % numero de ramas diferentes

if nargout<=2
    return
end

if nramas==0
    alto  = 0;
    ancho = 0;
else
    %cálculo de bounding box
    mx    = min(min(ramas(:,3:4)));   % x minimo
    Mx    = max(max(ramas(:,3:4)));   % x maximo
    my    = min(min(ramas(:,1:2)));   % y minimo
    My    = max(max(ramas(:,1:2)));   % y maximo
    alto  = Mx-mx;
    if alto<=minlong
        alto = 0;
    end
    ancho = My-my;
    if ancho<=minlong
        ancho = 0;
    end
end

if nargout<=4
    return
end

if isempty(ramas)
    maxiy = [];
    leaf = [];
else

    leaf = leafTree(ramas);

    if nargout<=5
        return
    end

    if useTransparency %USE TRANSPARENT BRANCHES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        unileaf = unique(leaf(:,1));
        lunileaf = length(unileaf);
        leafwithourrep = leaf;
        if lunileaf < length(leafwithourrep)
            for i = 1 : lunileaf
                rep(i) = sum(leafwithourrep(:,1) == unileaf(i));
            end
            frep = find(rep > 1);
            for j = 1 : length(frep)
                fleaf = find(leafwithourrep(:,1) == unileaf(frep(j)));
                [maxy pos] = max(leafwithourrep(fleaf,2));
                leafwithourrep(setdiff(fleaf,fleaf(pos)),:) = [];
            end
        end
        maxxl = max(leafwithourrep(:,1));
        maxiy = zeros(1,maxxl);
        maxiy(leafwithourrep(:,1)) = leafwithourrep(:,2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else %USE CANONICAL ALGORITHM, MIGHT BE SLIGHTLY MODIFIED FOR ALLOWING CHANCE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sr = size(r,1);
        maxx = max(max(r(:,1:2)));
        maxiy = zeros(1,maxx);
        for j = 1 : sr
            x0 = r(j,1);
            x1 = r(j,2);
            y0 = r(j,3);
            y1 = r(j,4);
            if x0 == x1
                y = max(y0,y1);
                if (maxiy(x0)< y) && chanceFun()
                    maxiy(x0)=y;
                end
            elseif x0<x1
                for k = x0:x1
                    y = y0 + (((k-x0)*(y1-y0))/(x1-x0));
                    if (maxiy(k) < y) && chanceFun()
                        maxiy(k) = y;
                    end
                end
            else
                for k = x1:x0
                    y = y0 + (((k-x0)*(y1-y0))/(x1-x0));
                    if (maxiy(k) < y) && chanceFun()
                        maxiy(k) = y;
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% función para contar las ramas y los nodos. las que se superponen solo se cuentan una vez
function filas = cuenta_filas(filas,nueva, firstcols) %#ok<DEFNU>

%explanation:
%  bsxfun(@eq, filas, nueva) == [filas(:,1)==nueva(1), filas(:,2)==nueva(2), ..., filas(:,firstcols)==nueva(firstcols)]
%  all(x, 2) == x(:,1) & x(:,2) & ... & x(:,firstcols)
%  ~any(z) == all elements in z are 0
%so, the expression ~any(all(bsxfun(@eq, filas, nueva), 2)) is true iff
%there is now row in 'filas' matching the single row 'nueva' in the first
%'firstcols' columns
if isempty(filas) || ~any(all(bsxfun(@eq, filas(:,1:firstcols), nueva(:,1:firstcols)), 2))
    filas = [filas;nueva];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function does the same as cuenta_filas, but makes sure to operate on
%an ordered list of rows, in order to (hopefully) be more efficient. An
%even more efficient version can be made if using binary search, but it
%would be very cumersome to implement taking into account special cases
function filas = sortedInsertNoRepeat(filas, nueva, firstcols)

if isempty(filas)
    filas = nueva;
else
    idx1 = 1;
    idx2 = size(filas,1);
    %test values for each column. After all loops, idx will point to the
    %first row equal to or greater than 'nueva'
    for k=1:firstcols
        eqs = nueva(k)==filas(idx1:idx2,k);
        newidx1 = find(eqs, 1, 'first');
        if isempty(newidx1) %there is no row with this k-th component: ordered insert
            pos = find(nueva(k)<filas(idx1:idx2,k), 1, 'first')+idx1-1;
            if isempty(pos) %k-th component is greater than all already inserted: to the end of the segment
                filas = [filas(1:idx2,:); nueva; filas((idx2+1):end,:)];
            else %insert k-th component ordered
                filas = [filas(1:pos-1,:); nueva; filas(pos:end,:)];
            end
            return;
        end
        newidx2 = find(eqs, 1, 'last');
        oldidx1 = idx1;
        idx1    = newidx1+oldidx1-1; %try it again in the range with this k-th component
        idx2    = newidx2+oldidx1-1;
    end
    if idx1~=idx2
        error('At this stage, the window must have been narrowed down to a single element of ''filas'' to be identical to ''nueva''!!!');
    end
    if ~all(nueva(1:firstcols)==filas(idx1,1:firstcols))
        error('The detected element of ''filas'' is expected to be the same as ''nueva''!!!');
    end
    %do nothing, to silently ignore the repeated element
end