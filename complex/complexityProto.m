function complexityProto(basedir,GenI,GenF)

for i = GenI : GenF
    comLZfigureLongitud = [];
    load(sprintf('%s/numGroups%03d.mat',basedir,i));
    numGroups = eval(sprintf('numGroups%03d', i));

    load(sprintf('%s/offsetproto%03d.mat',basedir,i));
    offsetproto = eval(sprintf('offsetproto%03d', i));

    load(sprintf('%s/proto%03d.mat',basedir,i));
    proto = eval(sprintf('proto%03d', i));

    for j = 1 : numGroups
        dimension = [max(proto{numGroups}{j}{1}) max(proto{numGroups}{j}{2})];
        array = booleanTree(proto{numGroups}{j},dimension);
        LZfigure = norm2lzw(uint8(array(:)));
        sLZfigure = length(LZfigure);
        comLZfigureLongitud(j) = sLZfigure;
    end
    tosave1 = ['comLZfigureLongitud' num2str(i, '%03g')];
    eval([tosave1 ' = comLZfigureLongitud;']);
    save([basedir filesep tosave1],tosave1);
    clear('comLZfigureLongitud', tosave1);
end

function array = booleanTree(raster,dimensions)

array = false(dimensions);

if isempty(raster)
  return
end

ys = raster{1}; if size(ys,2)>1; ys = ys(:); end;
xs = raster{2}; if size(xs,2)>1; xs = xs(:); end;

array(sub2ind(dimensions, ys, xs)) = true;