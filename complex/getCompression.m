function c = getCompression(raster, dimensions, level)


if isempty(raster) || isempty(raster{1})
  c = 0;
  return
end

ys = raster{1};
xs = raster{2};
bitmap = false(dimensions);
bitmap(ys+(dimensions(1)*(xs-1))) = true;
%bitmap(sub2ind(dimensions, ys, xs)) = true;

if not(all(size(bitmap)==dimensions))
  error('The individual has dimensions %s, but the bitmap has dimensions %s!!!!', mat2str(dimensions), mat2str(size(bitmap)));
end

bitmap = bitmap(end:-1:1,:);

c = compressedImgSize(bitmap, level);
