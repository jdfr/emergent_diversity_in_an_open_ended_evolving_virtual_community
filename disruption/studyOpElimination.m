function [pos ejey] = studyOpElimination(symbols)
mode = 'canonical';
for i = 1 : 10000
    i
    cad1 = [];
    cad = [];
    raster1 = [];
    raster2 = [];
    offset1 = [];
    offset2 = [];
    cad1 = generateGenome;
    [cad t(i)]= op_eliminacionG(cad1,symbols);
    [raster1, offset1]=ls2(cad1, mode);
    [raster2, offset2]=ls2(cad, mode);
    
    disrup(i) = (1 - similitud2(raster1,offset1,raster2,offset2));
end

a = unique(t);

for i = 1 : length(a)
    f{i} = find(t == a(i));
    pos(i) = a(i);
    ejey(i) = sum(disrup(f{i}))/length(f{i});
end


    
bar(pos,ejey);
xlabel('N/L, donde N es la posicion donde se inserta el segmento y L la longitud del genoma')
ylabel('disrupcion media')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cad t]= op_eliminacionG(cad1,symbols)
% symbols                              = ['G']; %['G','+','-','['];
perm                                 = randperm(numel(symbols));
operAplicada                         = false;

for z=1:numel(symbols)
    seg = symbols(perm(z));
    if(seg=='[')
        % implica eliminaci�n de un corchete vac�o []
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

function genome=generateGenome

sizeGenome = randint(1,1,[4 40]);

numBra = randint(1,1,[1 round(sizeGenome/2)]);
if mod(numBra,2)~=0
    numBra = numBra + 1;
end
permutaciones = randperm(sizeGenome);
posBra = sort(permutaciones(1:numBra));

symbols                              = ['G','+','-'];
genome(setdiff(1:sizeGenome,posBra)) = symbols(randint(1,length(setdiff(1:sizeGenome,posBra)),[1,length(symbols)]));
genome(posBra(1:2:end-1)) = '[';
genome(posBra(2:2:end)) = ']';

k1 = find(genome == '+');
k2 = find(genome == '-');
k3 = find(genome == 'G');
if (isempty(k1) && isempty(k2)) ||isempty(k3)
    genome(posBra(1)) = '+';
    genome(posBra(2)) = 'G';
end

