function treeToMatch = generateOddRaster1(heightAll, heightElbow, heightHorzCanal, widthAll, widthVertLeft, widthVertRight, fillMode)

if ~exist('fillMode', 'var')
  fillMode = 'full';
end

raster = false(heightAll, widthAll);

raster(heightElbow+1:heightElbow+heightHorzCanal,:)                      = true;
raster(heightElbow+heightHorzCanal+1:end,        1:widthVertLeft)        = true;
raster(1:heightElbow,                            end-widthVertRight:end) = true;

switch fillMode
  case 'full'
  case 'chess'
    for k=1:heightAll
      raster(k,mod(k,2)+1:2:widthAll) = false;
    end
  otherwise
    error('fillMode <%s> not understood!!!!', fillMode);
end

[ys xs]     = find(raster);
offset      = [1, widthAll-round(widthVertRight/2)];
% raster2=uint8(raster); raster2(offset(1), offset(2)) = 2; imshow(raster2(end:-1:1,:), [0 0 0; 1 1 1; 1 0 0; 0 1 0], 'InitialMagnification', 'fit');
% return;
dimension   = [heightAll, widthAll];
maxiy       = repmat(-1, 1, widthAll);
isLeaf      = false(size(maxiy));
raster      = {ys xs};
counts      = [0 0 0];
treeToMatch = {raster, offset, counts, dimension, isLeaf, maxiy};