function res = monsterGenomeMatrixDistData(pophs, type)
%for reduced genome comparison across several simulations

switch type
  case 'minGenome'
  case 'genome'
  otherwise
    error('type %s not understood!!!!', type);
end

mxs = 10000;

%[NG N[ N+ N- N ]
%[01 02 03 04 05]
datos = zeros(mxs, 5, 'int32'); 
indivs = cell(mxs, 1);
ni = mxs;

popindexes = cell(size(pophs));
gens       = cell(size(pophs));

n = 0;

%find unique individuals across all simulations and generations
for k=1:numel(pophs)
  poph = pophs{k};
  ming = min(poph.generation);
  maxg = max(poph.generation);
  popindexes{k} = cell(maxg-ming+1,1);
  gens{k} = (ming:maxg)';
  z=1;
  for g=ming:maxg
    fprintf('GEN %03g, %s\n', g, poph.ACIParams.nomdir);
    thisg = find(poph.generation==g);
    nsoil = poph.nPixelsInSoil(thisg);
    genomas = poph.(type)(thisg);
    Pdatos = genomeChars(genomas);
    indexes = zeros(numel(thisg),1);
    for m=1:numel(indexes)
      if nsoil(m)==0 %valid tree
        maybe = find( (datos(1:n,1)==Pdatos(m,1)) & (datos(1:n,2)==Pdatos(m,2)) & (datos(1:n,3)==Pdatos(m,3)) & (datos(1:n,4)==Pdatos(m,4)) );
        found = find(strcmp(genomas{m}, indivs(maybe)));
        if not(isempty(found))
          indexes(m) = maybe(found);
        else
          n=n+1;
          if n>ni
            datos(end+mxs,1)=0;
            indivs{end+mxs} = [];
            ni = numel(indivs);
          end
          datos(n,:) = Pdatos(m,:);
          indivs{n}  = genomas{m};
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

res = struct('indivs', {indivs}, 'datos', {datos}, 'popindexes', {popindexes}, 'gens', {gens});  


