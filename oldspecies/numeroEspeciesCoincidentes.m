function numeroEspeciesCoincidentes(basedir,G)

% G = number of generations

s = zeros(1,G-1);


for i = 1 : G
    load(sprintf('%s/proto%03d.mat',basedir,i));
    load(sprintf('%s/offsetproto%03d.mat',basedir,i));
    load(sprintf('%s/numGroups%03d.mat',basedir,i));
    numGrupos(i)=eval(sprintf('numGroups%03d', i));

end


for i = 1 : G-1
    c = 0;
    proto1t = eval(sprintf('proto%03d', i));
    proto2t = eval(sprintf('proto%03d', i+1));
    offsetproto1t = eval(sprintf('offsetproto%03d', i));
    offsetproto2t = eval(sprintf('offsetproto%03d', i+1));
    proto1 = proto1t{numGrupos(i)};
    proto2 = proto2t{numGrupos(i+1)};
    offsetproto1 = offsetproto1t{numGrupos(i)};
    offsetproto2 = offsetproto2t{numGrupos(i+1)};

    for j = 1 : numGrupos(i)

        for k = 1 : numGrupos(i+1)
            dist = 1 - similitud2(proto1{j},offsetproto1{j},proto2{k},offsetproto2{k});
            if dist <= 0.2
                c = c + 1;
            end
        end
    end
    s(i) = c;
end
h1 = figure;
plot(numGrupos,'k');
hold on

plot(s);
xlabel('generations')
ylabel('groups')

label('number of groups','number of coincident groups between the generation n and n+1')







function sim = similitud2(raster1,offset1,raster2,offset2)%, show)

if isempty([raster1{:}]) && isempty([raster2{:}])
    sim = 1;
elseif (isempty([raster1{:}]) && ~isempty([raster2{:}])) || (~isempty([raster1{:}]) && isempty([raster2{:}]))
    sim = 0;
else
    pixel1 = [raster1{:}];
    pixel2 = [raster2{:}];

    pixel1(:,1) = pixel1(:,1) - offset1(1);
    pixel1(:,2) = pixel1(:,2) - offset1(2);

    pixel2(:,1) = pixel2(:,1) - offset2(1);
    pixel2(:,2) = pixel2(:,2) - offset2(2);

    %     sraster1 = size(pixel1,1);
    %     sraster2 = size(pixel2,1);

    sim = size(intersect(pixel1,pixel2,'rows'),1)/size(union(pixel1,pixel2,'rows'),1);
end

