function sparseM = landscape(R,ran)

% R es un cell que contine la matriz del grafo pixelado de todos los individuos de la
% poblacion
if sum(size(R)) == 0
    s =1
end
se = size(R{1},2);
sR = size(R,2);
t = floor((se*sR)/2);
sparseM = sparse(zeros(size(R{1},1),t));

for i = 1 : sR
    if ~isempty(R{i})% puede ser vacio en l-systems del tipo [--][-]
        
        m = R{i}*i;
        pos = randint(1,1,[1,t-se+1]);
        for k = pos : se+pos-1
            for j = 1 : size(R{i},1)
                if sparseM(j,k) == 0
                    sparseM(j,k) = m(j,k-pos+1);
                elseif sparseM(j,k) ~= 0 &&  m(j,k-pos+1)~= 0
                    sparseM(j,k) = -1;
                end
            end
        end
    end
end

