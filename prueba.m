function prueba

bolas = 1:6;

n = 100000;
res1 = zeros(n,1);
res2 = zeros(n,1);

for k=1:n
  res2(k) = ceil(rand*numel(bolas));
  res1(k) = bolas(1);
  for z=2:numel(bolas)
    prob = (z-1)./z;
    rnd  = rand>prob;
    if rnd
      res1(k) = bolas(z);
    end
  end
end

figure;

subplot(2,1,1); hist(res1, 0:6); grid on;
subplot(2,1,2); hist(res2, 0:6); grid on;
