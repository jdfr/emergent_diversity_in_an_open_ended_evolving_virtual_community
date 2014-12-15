function paintPrototypes(eachgroup,proto,offsetproto,indicesClasificados,basedir,Generacion)

grandote = true;
separados = false;

leachgroup = length(eachgroup);
drawing.FontSize = 6;

m = ceil(sqrt(leachgroup));
if(m * (m-1) >= leachgroup)
    n = m-1;
else
    n = m;
end

dimension = zeros(leachgroup, 2);
notEmptyTrees = cellfun(@(x) (~isempty(x{1})), proto);

dimension(notEmptyTrees, :) = cell2mat(cellfun(@(x) ([max(x{1}) max(x{2})]), proto, 'UniformOutput', false)');

maxDim = max(dimension);

        h = figure('visible','off');

for i = 1 : leachgroup
 
    if(grandote)
        subAxis = [0 dimension(i,2) 0 dimension(i,1)];
    else
        subAxis = [[0 maxDim(2)]-((maxDim(2)-dimension(i,2))/2) 0 maxDim(1)];
    end
    if separados
        drawtree(proto{i},dimension(i,:),offsetproto{i},0,1,[0 0 0], subAxis);
        saveas(h,[basedir filesep 'grupos' num2str(Generacion, '%03g') num2str(i, '%03g')],'png');
        close(h);
    else
        subplot(m,n,i);%,'align');
        drawtree(proto{i},dimension(i,:),offsetproto{i},0,1,[0 0 0], subAxis);
    if ~isempty(indicesClasificados)
        title({mat2str(indicesClasificados(eachgroup{i})),['DIMENSION: ',mat2str(dimension(i,:))]},'FontSize',drawing.FontSize);
    else
        title({mat2str(eachgroup{i}),['DIMENSION: ',mat2str(dimension(i,:))]},'FontSize',drawing.FontSize);
    end
    end

end
if ~separados
    saveas(h,[basedir filesep 'grupos' num2str(Generacion, '%03g')],'png');
    close(h);
end
drawnow;
