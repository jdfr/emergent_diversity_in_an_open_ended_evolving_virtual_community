function especie=agrupaciones(basedir,numSpecies1,individualsSpecies1,protos1,offsets1,numSpecies2,individualsSpecies2,protos2,offsets2,generation)

novacios1 = find(numSpecies1{generation}~=0);
novacios2 = find(numSpecies2{generation}~=0);

individualsSpecies1 = individualsSpecies1{generation}(novacios1);
proto1 = protos1{generation}(novacios1);
offset1 = offsets1{generation}(novacios1);

individualsSpecies2 = individualsSpecies2{generation}(novacios2);
proto2 = protos2{generation}(novacios2);
offset2 = offsets2{generation}(novacios2);


for i = 1 : length(novacios1)
    i
    for j = 1 : length(individualsSpecies1{i}) 
        j
        ind = individualsSpecies1{i}(j);
        for k = 1 : length(novacios2)
            if ~isempty(find(individualsSpecies2{k}==ind));
                especie{i}(j)=k;
            end
        end
    end
end

for i = 1 : length(novacios2)
    spe{i} = [];
    for j = 1 : length(novacios1)
            if length(find(especie{j} == i))>0
                spe{i} =[spe{i} j];
            end
    end
end


for i = 1 : length(novacios2)
    h = figure;
    lagru = length(spe{i});
    m = ceil(sqrt(lagru));
    if(m * (m-1) >= lagru)
        if lagru == 0
            n = 1;
        else
            n = m-1;
        end
    else
        n = m;
    end
    m = m + 1;

    subplot(m,n,1:n);
    dimension = [max(proto2{i}{1}) max(proto2{i}{2})];
    if isempty(dimension)
        dimension = zeros(1,2);
    end
    %     subAxis = [0 dimension(1,1) 0 dimension(1,2)];
    drawtree(proto2{i},dimension,offset2{i},0,1,[0 0 0], 'fit');
    title({mat2str(individualsSpecies2{i}),['DIMENSION: ',mat2str(dimension)]},'FontSize',6);
    for j = 1 : lagru
        subplot(m,n,n+j);
        dimension = [max(proto1{spe{i}(j)}{1}) max(proto1{spe{i}(j)}{2})];
        if isempty(dimension)
            dimension = zeros(1,2);
        end

        %         subAxis = [0 dimension(1,1) 0 dimension(1,2)];
        drawtree(proto1{spe{i}(j)},dimension,offset1{spe{i}(j)},0,1,[0 0 0], 'fit');
        title({['ESPECIE: ',mat2str(spe{i}(j))],['DIMENSION: ',mat2str(dimension)]},'FontSize',6);
    end
    drawnow;
saveas(h,[basedir,'especie' num2str(i, '%03g')],'png');
end




        
        
        
