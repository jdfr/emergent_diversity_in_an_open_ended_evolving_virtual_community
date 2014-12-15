function g=getCommonAncestorGen(poph, gen, idx1, idx2)
n1 = 0;
n2 = 0;
i1 = find((poph.tree.gD==gen)&(poph.tree.iD==idx1));
i2 = find((poph.tree.gD==gen)&(poph.tree.iD==idx2));
while i1~=i2
  i1 = find(poph.tree.idxD==poph.tree.idxA(i1));
  i2 = find(poph.tree.idxD==poph.tree.idxA(i2));
end
g = poph.tree.gD(i1);
