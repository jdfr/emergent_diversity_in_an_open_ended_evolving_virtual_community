function showGenomes(axioms, mode, newfig)
%for several genomes, show their corresponding phenotypes
angle=22.5;
tumbar=false;
%1 3 8 26 27

if not(exist('newfig', 'var'))
  newfig = false;
end

for z=1:numel(axioms)
  [raster, offset, counts]=ls2(axioms{z}, mode);
  dimension = [max(raster{1}) max(raster{2})];
  array = drawraster(raster, dimension);
  array=uint8(array);
  %array(dimension(1)-offset(1)+1, offset(2))=2; %STEM
  if tumbar
    array=array';
  end
  if newfig
    figure;
  end
  imshow(array, [0 0 0; 1 1 1;1 0 0]);
  title(sprintf('%s, idx: %d/%d', mat2str(dimension), z, numel(axioms)));
  pause;
end

function array = drawraster(raster, dimensions)
ys = raster{1}; if size(ys,2)>1; ys = ys(:); end;
xs = raster{2}; if size(xs,2)>1; xs = xs(:); end;
%array = accumarray([ys,xs],uint16(indexInColors+2), uint16(dimensions));
array = true(dimensions);
array(sub2ind(dimensions, ys, xs)) = false;
array = array(end:-1:1,:);
