function saveTree(axiom, fname, mode, varargin)

if ~exist('mode', 'var')
  mode = 'canonical';
end
[raster, offset, counts, dimension, isLeaf, maxiy, leafs]=ls2(axiom, mode);

rastermaster = zeros(dimension, 'uint8');

rastermaster(sub2ind(size(rastermaster), raster{1}, raster{2})) = 1;
rastermaster(sub2ind(size(rastermaster), leafs{1},  leafs{2}))  = 2;

imwrite(rastermaster(end:-1:1, :), [1 1 1; 0 0 0; 1 0 0; 1 0 0], fname);
