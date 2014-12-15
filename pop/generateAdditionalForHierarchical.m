function lastgen = generateAdditionalForHierarchical(basedir)

if exist([basedir filesep 'poph.mat'], 'file')
  poph = load([basedir filesep 'poph.mat']);
  poph = poph.poph;
  changed = false;
else
  poph = loadSimulation(basedir);
%  save([basedir filesep 'poph.mat'], 'poph');
  changed = true;
end

if not(isfield(poph, 'minGenome'))
  fprintf('sanitizing...\n');
  poph.minGenome = cellfun(@sanitizeString, poph.genome, 'uniformoutput', false);
  changed = true;
  fprintf('sanitized...\n');
end

lastgen = max(poph.generation);

if changed
  fprintf('saving poph...\n');
  save([basedir filesep 'poph.mat'], 'poph');
  fprintf('poph saved...\n');
end

if not(exist([basedir filesep 'mutdists.mat'], 'file'))
  fprintf('calculating numuts...\n');
  [numMuts mutDists] = makeNumMutations(poph);
  mutDists = cellfun(@(x)uint16(x), mutDists, 'uniformoutput', false);
  fprintf('saving numuts...\n');
  save([basedir filesep 'mutdists.mat'], 'numMuts', 'mutDists');
  fprintf('numuts saved...\n');
end


function d=mutDist(poph, gen, idx1, idx2)
d=poph.tree.mutDists{gen+1}(idx1, idx2);

function d=mutDist2(poph, gen, idx1, idx2)
n1 = 0;
n2 = 0;
i1 = find((poph.tree.gD==gen)&(poph.tree.iD==idx1));
i2 = find((poph.tree.gD==gen)&(poph.tree.iD==idx2));
while i1~=i2
  n1 = n1+poph.tree.numMuts(i1);
  n2 = n2+poph.tree.numMuts(i2);
  i1 = find(poph.tree.idxD==poph.tree.idxA(i1));
  i2 = find(poph.tree.idxD==poph.tree.idxA(i2));
end
d=n1+n2;

