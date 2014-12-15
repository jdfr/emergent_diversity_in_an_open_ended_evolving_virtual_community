

function res = complexityProtoNEWCOMPLEXITY(basedir,GenI,GenF, level)

generations = (GenI : GenF)';
comLZfigureLongitudALL = cell(size(generations));

for z = 1:numel(generations)
  i=generations(z);
  fprintf('generation %d -> %d -> %d\n',  GenI, i, GenF);
    load(sprintf('%s/numGroups%03d.mat',basedir,i));
    numGroups = eval(sprintf('numGroups%03d', i));
    clear(sprintf('numGroups%03d', i));

%     load(sprintf('%s/offsetproto%03d.mat',basedir,i));
%     offsetproto = eval(sprintf('offsetproto%03d', i));
  
    fprintf('loading proto ...\n');
    load(sprintf('%s/proto%03d.mat',basedir,i));
    proto = eval(sprintf('proto%03d', i));
    clear(sprintf('proto%03d', i));
    proto = proto{numGroups};
    comLZfigureLongitud = zeros(numGroups,1);
    fprintf('calculating complexities ...\n');
    for j = 1 : numGroups
        dimension = [max(proto{j}{1}) max(proto{j}{2})];
        comLZfigureLongitud(j) = getCompression(proto{j}, dimension, level);
    end
    comLZfigureLongitudALL{z} = comLZfigureLongitud;
    clear('comLZfigureLongitud');
end

res = struct('gs', {generations}, 'complexities', {comLZfigureLongitudALL});
