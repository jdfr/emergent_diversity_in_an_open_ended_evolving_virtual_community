function d = permuteSQM(i, o)
%given a square matrix and a new list of rows/cols, it makes the new matrix
%from the old one
d = zeros(numel(i));
for k=1:numel(i)
  d(:,k) = o(i, i(k));
end
