%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawTwoTreesTogether(axiom1, axiom2, T)

if nargin<3
  T = 'transparentBranches';
end

[raster1, offset1, counts1, dimension1, isLeaf1, maxiy1]=ls2(axiom1, T); %#ok<ASGLU,NASGU>
[raster2, offset2, counts2, dimension2, isLeaf2, maxiy2]=ls2(axiom2, T); %#ok<NASGU>

% if max(dimension1)>max(dimension2)
  colors = [2 3];
  rastermaster = zeros(dimension1, 'uint8');
  rastermaster(sub2ind(size(rastermaster), raster1{1}, raster1{2})) = colors(1);
  offsetMaster = offset1;
  rasterTree = raster2;
  offsetTree = offset2;
% else
%   colors = [3 2];
%   rastermaster = zeros(dimension2, 'uint8');
%   rastermaster(sub2ind(size(rastermaster), raster2{1}, raster2{2})) = colors(1);
%   offsetMaster = offset2;
%   rasterTree = raster1;
%   offsetTree = offset1;
% end

ys = rasterTree{1}-offsetTree(1)+offsetMaster(1);
xs = rasterTree{2}-offsetTree(2)+offsetMaster(2);

minY = min(ys);
minX = min(xs);

if minY<=0
  rastermaster = [zeros(1-minY, size(rastermaster,2), 'uint8'); rastermaster];
  ys = ys+1-minY;
end
if minX<=0
  rastermaster = [zeros(size(rastermaster,1), 1-minX, 'uint8'), rastermaster];
  xs = xs+1-minX;
end

maxY = max(ys);
maxX = max(xs);

if maxY>size(rastermaster,1)
  rastermaster = [rastermaster; zeros(maxY-size(rastermaster,1), size(rastermaster,2), 'uint8')];
end
if maxX>size(rastermaster,2)
  rastermaster = [rastermaster, zeros(size(rastermaster,1), maxX-size(rastermaster,2), 'uint8')];
end

rastermaster(sub2ind(size(rastermaster), ys, xs)) = colors(2);
%n=5;rastermaster=rastermaster(:,reshape(repmat(1:size(rastermaster,2), n, 1),[], 1));
imshow(rastermaster, [1 1 1; 0 0 0; 1 0 0; 0 0 1], 'InitialMagnification', 'fit');
set(gca, 'YDir', 'normal');
