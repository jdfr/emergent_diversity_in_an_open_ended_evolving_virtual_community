function scriptShowTree
%a hack to show the development of a tree, or several trees, or several
%trees and their leaves...

angle=19;22.5;
angle = 22;
nreps=3;
tumbar=false;
axioms = {...
  '++G', ...
  ...'[+GG]G', ...
  ...'G[--GG][G]', ...
  ...'G[G][--G][++G]', ...
  ...'G[G][--G]++', ...
  ...
  ...'G-[++G]-G+[-G]+G',...
  ...%'G+G+';...%'G-[++G]-G+[-G]+G'
...%   '[[+G]G[G]+G]', ...
...%   '[+G]G[G]+G', ...
  ...%'[+G]G'; ...
...%    'GG[+G[G]-G[-G-G]++G]G';...
...%    'GG[+[G]-G[-G-G]++G]G';...
...%    'G[G+[G][-G]]';...
...%    'G[G+[G][--G]]';...
...%    'G[G+[G+GG][-[-GG]+GG]]G';...
  ...%'GGGGG[GG[-G-G]G][GG][G+[G+GG][G[GG]+GG][G+GG+G]GG][GGG+G][G+G][GG]G';...
  ...%'[GG-][[+][-][-][-]][+][-][+]G[-]';...
  ...%'GGGG[-G[-GG][-G[-G[G][+GG]-GG][-G]+G][G][+GG]-GG][-G][+G]';...
  ...%'G[++G][--G][G]'; ... %centro
  ...%'G[--GG][G]'; ... %izquierda
  ...%'G[G][--G]++'; ... %derecha
  };

% exps1 = expandString(axioms{1}, true);
% [raster1, offset1]=treeRaster(exps1{end}, struct('shadowing', 'canonical', 'angle', 22, 'justPaint', true));
% dimension1 = [max(raster1{1}) max(raster1{2})];
% array1 = drawraster(raster1, dimension1);
% exps2 = expandString(axioms{2}, true);
% [raster2, offset1]=treeRaster(exps2{end}, struct('shadowing', 'canonical', 'angle', 22, 'justPaint', true));
% dimension2 = [max(raster2{1}) max(raster2{2})];
% array2 = drawraster(raster2, dimension2);
% 
% if size(array1,1)>size(array2,1)
%   array2 = [true(size(array1,1)-size(array2,1),size(array2,2)); array2];
% end
% if size(array2,1)>size(array1,1)
%   array1 = [true(size(array2,1)-size(array1,1),size(array1,2)); array1];
% end
% 
% zz=20;%round(((size(array1,2)+size(array2,2))/10));
% array = [array1, true(size(array1,1), zz), array2];
% 
% z=20;
% array = [true(size(array,1),z), array, true(size(array,1),z)];
% array = [true(z, size(array,2)); array; true(z, size(array,2))];
% 
% %imwrite
% 
% h = figure;
% 
% fs=20;
% if size(array1,2)>size(array2,2)
%   array2 = [array2 true(size(array2,1), size(array1,2)-size(array2,2))];
% end
% if size(array2,2)>size(array1,2)
%   array1 = [array1 true(size(array2,1), size(array2,2)-size(array1,2))];
% end
% subplot(1,2,1);imshow(array1); xlabel(axioms{1}, 'fontsize', fs);
% subplot(1,2,2);imshow(array2); xlabel(axioms{2}, 'fontsize', fs);
% 
% img = imcapture(h, 'all', 600);
% imwrite(img, 'exampledisruption.png');
% close(h);
mxdim = [0 0];
for z=1:numel(axioms)
  axiom = axioms{z};
%   fprintf('AXIOM:         %s\n', axiom);
  exps = [{'G'}; expandString(axiom, true, nreps)];
%   for k=1:numel(exps);
%     fprintf('   EXPANDED %d: %s\n', k, exps{k});
%   end
for m=1:numel(exps)
%   [raster, offset, counts, dimension, isLeaf, maxiy, leafs]=treeRaster(exps{end}, struct('shadowing', 'canonical', 'angle', 19));
%   fprintf('NRAMAS: %d, NCHARS: %s, NGS: %s\n', counts(1), mat2str(cellfun('prodofsize', exps)), mat2str(cellfun(@(x)sum(x=='G'), exps)));
%   [raster, offset]=treeRaster(exps{m}, struct('shadowing', 'canonical', 'angle', 19, 'justPaint', true));
  [raster, offset counts, dimension, isLeaf, maxiy, leafs]=treeRaster(exps{m}, struct('shadowing', 'transparentBranches', 'angle',angle, 'justPaint', false, 'sameLevelRule', false, 'transparentBranches', true));
  dimension = [max(raster{1}) max(raster{2})];

  array = drawraster(raster, dimension, leafs);
  array=uint8(array);
  %array(dimension(1)-offset(1)+1, offset(2))=2;
  h=figure('name', sprintf('N+,NG,NR,NL: %s,%s,%s,%s', mat2str(sum(exps{m}=='+')), mat2str(sum(exps{m}=='G')), mat2str(counts(1)), mat2str(counts(2)) ));
  if tumbar
    array=array';
  end
  imshow(array, [0 0 0; 1 1 1;1 0 0]);
%   arr = '{\downarrow}';
%   vals = [exps(:)'; exps(:)'];
%   vals(2,:) = {arr};
%   title(sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s', 'G', arr, vals{:}), 'fontsize', 14);
  title(exps{m});
  fprintf('%s\n', exps{m});
%   img = imcapture(h, 'all', 900);
%   imwrite(img, [axiom '.png']);
%   close(h);
arrays{m} = array;
mxdim = max(mxdim, double(dimension));
end
% arrays{z} = array;
%mxdim = max(mxdim, double(dimension));
  
  
%   tiempo = 0;
%   [colors array]= drawtree(raster, dimension,offset,tiempo);
  %title(mat2str(dimension));
end


% close all hidden;
% sheet = ones(mxdim);
% arrays = cellfunc(@(x)putOnSheet(x, sheet), arrays);
% arrays = cellfunc(@(x)magnify(x, 10), arrays);
% for k=1:numel(arrays);
%    imwrite(logical(arrays{k}), sprintf('minitree_%d.png', k));
% end

%find . -type f -name \*.m -exec grep -i "381" '{}' \; -print

function array = drawraster(raster, dimensions, leafs)
ys = raster{1}; if size(ys,2)>1; ys = ys(:); end;
xs = raster{2}; if size(xs,2)>1; xs = xs(:); end;
%array = accumarray([ys,xs],uint16(indexInColors+2), uint16(dimensions));
array = ones(dimensions, 'uint8');
rsinds = sub2ind(dimensions, ys, xs);
array(rsinds) = false;
if exist('leafs', 'var')
  lsinds = sub2ind(dimensions, leafs{1}, leafs{2});
  array(lsinds) = 2;
end
array = array(end:-1:1,:);

function frame = putOnSheet(img, sheet)
  sizesheet = size(sheet);
  pxthis            = sizesheet(1)-size(img,1)+1;
  pythis            = floor((sizesheet(2)-size(img,2))/2+1);
  frame             = copyOn(sheet, [pxthis pythis], img);

function frame = putOnSheetIzq(img, sheet)
  sizesheet = size(sheet);
  pxthis            = sizesheet(1)-size(img,1)+1;
  pythis            = sizesheet(2)-size(img,2)+1;
  frame             = copyOn(sheet, [pxthis pythis], img);

function newarray = magnify(array, n)
newarray = zeros(size(array)*n, class(array));
qq = [1 n];
for q=1:size(array,1)
  ww = [1 n];
  for w=1:size(array,2)
    newarray(qq(1):qq(2),ww(1):ww(2)) = array(q,w);
    ww = ww + n;
  end
  qq = qq + n;
end
