function [numMuts mutDists] = makeNumMutations(poph)
%get data on mutation distance between individuals
numMuts = zeros(size(poph.tree.change));
change = poph.tree.change;
for k=1:numel(numMuts)
  c = upper(change{k});
  numMuts(k) = sum((c>='A')&(c<='Z')&(c~='G'));
end
mutDists = cell(max(poph.generation)+1,1);
gen=0;
mutDists{1} = nan(sum(poph.generation==gen));
for z=1:size(mutDists{1},1)
  mutDists{1}(z,z) = 0;
end
for k=2:numel(mutDists)
  gen = gen+1;
  mm = mutDists{k-1};
  thisGen = find(poph.tree.gD==gen);
  m = zeros(numel(thisGen));
  z = size(m,1);
  inds=poph.tree.iD(thisGen);
  if not(all(inds==(1:z)'))
    error('not consecutive!!\n%s', mat2str(inds));
  end
  for a=1:z
    for b=1:z
      m(a,b) = mm(poph.tree.iA(thisGen(a)), poph.tree.iA(thisGen(b)))+numMuts(thisGen(a))+numMuts(thisGen(b));
    end
    m(a,a)=0;
  end
  mutDists{k}=m;
end
