function sim=compara(axiom1, axiom2, T)

sim = similitudFastGenoma(axiom1,axiom2,T)
if nargin<3
  T = 'transparentBranches';
end

[raster1, offset1, counts1, dimension1, isLeaf1, maxiy1]=ls2(axiom1, T); %#ok<ASGLU,NASGU>
[raster2, offset2, counts2, dimension2, isLeaf2, maxiy2]=ls2(axiom2, T); %#ok<NASGU>

% if max(dimen)>max(dimension2)
colors = [2 3 4];
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

tree = rastermaster;
tree(find(rastermaster>0))=0;


%%%
rastermaster2 = zeros(dimension2, 'uint8');

  rastermaster2(sub2ind(size(rastermaster2), raster2{1}, raster2{2})) = colors(1);
  offsetMaster2 = offset2;
  rasterTree2 = raster1;
  offsetTree2 = offset1;

ys2 = rasterTree2{1}-offsetTree2(1)+offsetMaster2(1);
xs2 = rasterTree2{2}-offsetTree2(2)+offsetMaster2(2);

minY = min(ys2);
minX = min(xs2);

if minY<=0
  rastermaster2 = [zeros(1-minY, size(rastermaste2,2), 'uint8'); rastermaster2];
  ys2 = ys2+1-minY;
end
if minX<=0
  rastermaster2 = [zeros(size(rastermaster2,1), 1-minX, 'uint8'), rastermaster2];
  xs2 = xs2+1-minX;
end

maxY = max(ys2);
maxX = max(xs2);

if maxY>size(rastermaster2,1)
  rastermaster2 = [rastermaster2; zeros(maxY-size(rastermaster2,1), size(rastermaster2,2), 'uint8')];
end
if maxX>size(rastermaster2,2)
  rastermaster2 = [rastermaster2, zeros(size(rastermaster2,1), maxX-size(rastermaster2,2), 'uint8')];
end

ind1 = sub2ind(size(rastermaster2), ys2, xs2);
tree1(ind1)=colors(1);

%%%%


rastermaster(sub2ind(size(rastermaster), ys, xs)) = colors(2);


ind2 = sub2ind(size(rastermaster), ys, xs);
int = intersect(ind1,ind2);
rastermaster(int) = colors(3);

tree1 = tree;
tree2 = tree;
tree3 = tree;

tree1(ind2) = colors(2);
tree2(ind1) = colors(1);
tree3(int) = colors(3);

h1 = figure;
imshow(tree2, [1 1 1; 0 0 0;0.85 0.85 0.85; 0.6 0.6 0.6 ;0 0 0], 'InitialMagnification', 'fit');
set(gca, 'YDir', 'normal');

h2 = figure;
imshow(tree1, [1 1 1; 0 0 0;0.85 0.85 0.85; 0.6 0.6 0.6 ;0 0 0], 'InitialMagnification', 'fit');
set(gca, 'YDir', 'normal');

h3= figure;
imshow(tree3, [1 1 1; 0 0 0;0.85 0.85 0.85; 0.6 0.6 0.6 ;0 0 0], 'InitialMagnification', 'fit');
set(gca, 'YDir', 'normal');

h4 = figure;
imshow(rastermaster, [1 1 1; 0 0 0;0.85 0.85 0.85; 0.6 0.6 0.6 ;0 0 0], 'InitialMagnification', 'fit');
set(gca, 'YDir', 'normal');

h5 = figure;
subplot(2,2,1)
imshow(tree2, [1 1 1; 0 0 0;0.85 0.85 0.85; 0.6 0.6 0.6 ;0 0 0], 'InitialMagnification', 'fit');
set(gca, 'YDir', 'normal');
title(['tree1: ',axiom1]);
subplot(2,2,2)
imshow(tree1, [1 1 1; 0 0 0;0.85 0.85 0.85; 0.6 0.6 0.6 ;0 0 0], 'InitialMagnification', 'fit');
set(gca, 'YDir', 'normal');
title(['tree2: ',axiom2]);
subplot(2,2,3)
imshow(tree3, [1 1 1; 0 0 0;0.85 0.85 0.85; 0.6 0.6 0.6 ;0 0 0], 'InitialMagnification', 'fit');
set(gca, 'YDir', 'normal');
title('mixed');

subplot(2,2,4)
imshow(rastermaster, [1 1 1; 0 0 0;0.85 0.85 0.85; 0.6 0.6 0.6 ;0 0 0], 'InitialMagnification', 'fit');
set(gca, 'YDir', 'normal');

