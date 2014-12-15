function [lmingenome comLZfigureLongitud]=longGenomavsCom(basedir,Generation)

load(sprintf('%s/comLZfigureLongitud%03d.mat',basedir,Generation));
comLZfigureLongitud = eval(sprintf('comLZfigureLongitud%03d', Generation));

speciationResumen = load([basedir filesep 'speciationResumen.mat']);

s = size(speciationResumen.zdata.groups.individualsByProto{Generation},2);

for i = 1 : s
    tam = [];
    genomes = speciationResumen.zdata.groups.individualsByProto{Generation}{i}.genomes;
    raster{i}=speciationResumen.zdata.groups.individualsByProto{Generation}{i}.raster;
    dimensions(i,:)=speciationResumen.zdata.groups.individualsByProto{Generation}{i}.dimensions;
    offset{i}=speciationResumen.zdata.groups.individualsByProto{Generation}{i}.offset;

    for j = 1 : size(genomes,1)
        tam(j) = length(genomes{j});
    end
    lmingenome(i) = min(tam);
end

dim = [max(dimensions(:,1)),max(dimensions(:,2))];

for i = 1 : s
    h = figure('visible','off');
    drawtree(raster{i}{1},dim,offset{i},0,1,[0 0 0]);
    drawnow
    saveas(h,[basedir filesep 'proto' num2str(Generation, '%03g') num2str(i, '%03g')],'png');
    close(h)
end
tosave1 = ['lmingenome' num2str(Generation, '%03g')];
eval([tosave1 ' = lmingenome;']);
save([basedir filesep tosave1],tosave1);

plot(lmingenome,comLZfigureLongitud,'k.','MarkerSize',20);

xlabel('genome length');
ylabel('complexity');