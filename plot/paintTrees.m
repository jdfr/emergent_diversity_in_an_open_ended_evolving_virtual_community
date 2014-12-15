function allarray = paintTrees(rasters, dimensions, doshow)
%show several trees together

arrays = cell(size(rasters));

m1 = 0;
m2 = 0;

for k=1:numel(arrays)
  raster = rasters{k};
  if iscell(raster)
    dimension = dimensions(k,:);
    ys = raster{1}; if size(ys,2)>1; ys = ys(:); end;
    xs = raster{2}; if size(xs,2)>1; xs = xs(:); end;
    array = true(dimension);
    array(sub2ind(dimension, ys, xs)) = false;
    array = array(end:-1:1,:);
  else
    array = raster;
  end
  arrays{k} = array;
  m1 = max(m1, size(array,1));
  m2 = max(m2, size(array,2));
end

for k=1:numel(arrays)
  s = size(arrays{k});
  if s(1)<=m1
    arrays{k} = [true(m1-s(1)+1, s(2)); arrays{k}];
    %arrays{k}(end+1:end+(m1-s(1)+1),:) = true;
    %arrays{k}(end+(m1-s(1)+1),1) = false;
  end
%   if s(2)<=m2
%     %arrays{k}(1,end+(m2-s(2)+1)) = false;
%     arrays{k}(:,end+1:end+(m2-s(2)+1)) = true;
%   end
end

allarray = horzcat(arrays{:});

if doshow
  showSuperLandscapeIMG(allarray, 800);
%figure;
%imshow(allarray);
end
