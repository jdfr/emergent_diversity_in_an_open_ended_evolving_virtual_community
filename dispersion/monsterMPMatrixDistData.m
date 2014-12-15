function res = monsterMPMatrixDistData(dirs)
%for morphological (jaccard, lobrnaches) distances across several simulations

mxs = 5000;

%[siz1 siz2 off1 off2 nramas nleafs npixelsneg]
%[01   02   03   04   05     06     07        ]
datos = zeros(mxs, 7, 'int32'); 
indivs = cell(mxs, 1);
ni = mxs;

popindexes = cell(size(dirs));
gens       = cell(size(dirs));

n = 0;

%add initial individual
n=n+1;
P000 = load([dirs{1} filesep 'P000.mat']);
P000 = P000.P000;
datos(n,:) = [P000.dimensions, P000.offset, P000.counts];
indivs{n} = P000.raster{1};
clear P000;

%find unique individuals across all simulations and generations
for k=1:numel(dirs)
  d = dir([dirs{k} filesep 'P*.mat']);
  d = cellfun(@(x)str2double(x(2:4)), {d.name}');
  ming = min(d);
  maxg = max(d);
  popindexes{k} = cell(maxg-ming+1,1);
  gens{k} = (ming:maxg)';
  z=1;
  for g=ming:maxg
    fprintf('GEN %03g, %s\n', g, dirs{k});
    pn = sprintf('P%03g', g);
    P = load([dirs{k} filesep pn '.mat']);
    P = P.(pn);
    Pdatos = int32([P.dimensions, P.offset, P.counts]);
    indexes = zeros(numel(P.raster),1);
    for m=1:numel(P.raster)
      if Pdatos(m,end)==0 %valid tree
        maybe = find( (datos(1:n,1)==Pdatos(m,1)) & (datos(1:n,2)==Pdatos(m,2)) & (datos(1:n,3)==Pdatos(m,3)) & (datos(1:n,4)==Pdatos(m,4)) & (datos(1:n,5)==Pdatos(m,5)) );
        found = false;
        Rm = P.raster{m};
        for q = 1:numel(maybe)
          found = isequal(indivs{maybe(q)}, Rm);
          if found
            break
          end
        end
        if found
          indexes(m) = maybe(q);
        else
          n=n+1;
          if n>ni
            datos(end+mxs,1)=0;
            indivs{end+mxs} = [];
            ni = numel(indivs);
          end
          datos(n,:) = Pdatos(m,:);
          indivs{n}  = Rm;
          indexes(m) = n;
        end
      end
    end
    popindexes{k}{z} = indexes;
    z=z+1;
  end
end

indivs = indivs(1:n);
datos  = datos(1:n,:);

res = struct('indivs', {indivs}, 'datos', {datos}, 'popindexes', {popindexes}, 'gens', gens);  


