function out = packLogical(matrix)

[a b] = size(matrix);

out = zeros(a,b/8,'uint8');

za = 1;
zb = 8;
for k=1:size(out,2)
  chunk = matrix(:,za:zb);
  chunk = bwpack(chunk')';
  out(:,k) = uint8(chunk);
  za = za+8;
  zb=zb+8;
end  