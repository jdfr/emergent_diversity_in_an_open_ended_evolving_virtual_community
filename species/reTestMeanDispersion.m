function [mmins allmmins allallmmins allallallmmins] = reTestMeanDispersion(matrixDist, indicesClasificados, clusterss)
%this was coded to demonstrate myself that I didn't make any fucking 
%mistake when I implemented Rene's proposal of the clustering algorithm
%(namely, using the minimal mean disperson instead of the minimal sum of
%dispersions). I WAS FUCKING RIGHT

idx = 1:max(indicesClasificados);
idx(indicesClasificados) = 1:numel(indicesClasificados);

mmins = zeros(size(clusterss));
allmmins = cell(size(clusterss));
allallmmins = cell(size(clusterss));
allallallmmins = cell(size(clusterss));
for k=1:numel(clusterss)
  clusters = clusterss{k};
  minds = zeros(size(clusters));
  allallmmins{k} = cell(size(clusters));
  allallallmmins{k} = cell(size(clusters));
  for z=1:numel(clusters)
    cluster = idx(clusters{z});
    dispersions = zeros(size(cluster));
    allallallmmins{k}{z} = cell(size(cluster));
    for m=1:numel(cluster)
      dispersions(m) = mean(matrixDist(cluster(m), cluster));
      allallallmmins{k}{z}{m} = matrixDist(cluster(m), cluster);
    end
    allallmmins{k}{z} = dispersions;
    mind = min(dispersions);
    minds(z) = mind;
  end
  mmin = mean(minds);
  mmins(k) = mmin;
  allmmins{k} = minds;
end
