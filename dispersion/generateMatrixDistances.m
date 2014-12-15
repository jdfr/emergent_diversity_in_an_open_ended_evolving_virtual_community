function matrixDistss = generateMatrixDistances(dirs, res, dists, prefix)
%from structures representing a giant distance matrix between all
%individuals in one or several simulations, construct the lists of
%distance matrices for each simulation. The arguments for this function
%come from monsterMPMatrixDistData and
%monsterMPMatrixDist/monsterGenomeMatrixDist

popindexes = res.popindexes;
gens       = res.gens;

matrixDistss = cell(size(dirs));
for k=1:numel(dirs)
  fprintf('Going for %s...\n', dirs{k});
  matrixDist = cell(numel(gens{k})-1,1);
  for z=1:numel(matrixDist);
    p = popindexes{k}{z+1};
    p = p(p>0);
    matrixDist{z} = dists(p, p);
  end
  matrixDistss{k} = matrixDist;
  save([dirs{k} filesep prefix 'matrixDist.mat'], 'matrixDist');
end
