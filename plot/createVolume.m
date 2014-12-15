function createVolume(res, space2, k)
%a hack to show how a bidimensional measure of a population through time
%(3rd dimension). Not very helpful

popindexes = res.popindexes{k};
gens       = res.gens{k};

n=50;

xgrid = linspace(min(space2(:,1)), max(space2(:,1)), n);
ygrid = linspace(min(space2(:,2)), max(space2(:,2)), n);


xygrid = {xgrid, ygrid};

volume = zeros(n, n, numel(popindexes));

for z=1:numel(popindexes);
  inds = popindexes{z};
  inds = inds(inds>0);
  volume(:,:,z) = hist3(space2(inds,:), xygrid);
end

[x y z] = meshgrid(xgrid,ygrid,gens);
p = patch(isosurface(x,y,z,volume,1));
x=x;