function [colors array]= drawtree(raster,dimensions,offset,tiempo,indexInColors,colors, subAxis, leafs)

% dimensions = [1240 669]
% colors = [0 0 0]
if isempty(raster)
  return
end

if ~exist('colors', 'var')
    colors = rand(1,3)*0.8;
    indexInColors = 1;
end

ys = raster{1}; if size(ys,2)>1; ys = ys(:); end;
xs = raster{2}; if size(xs,2)>1; xs = xs(:); end;
%array = accumarray([ys,xs],uint16(indexInColors+2), uint16(dimensions));
array = zeros(dimensions, 'uint16');
array(sub2ind(dimensions, ys, xs)) = uint16(indexInColors+1);
yOffset = offset(1);
if yOffset>1;
  array(yOffset-1,:) = 2;
end
if exist('leafs', 'var')
  array(sub2ind(size(array), leafs{1}, leafs{2})) = indexInColors+1+2;
end
array = array(end:-1:1,:);
if size(array,2)==1
  array(end,2)=0;
end
if (~exist('subAxis', 'var')) || isempty(subAxis) || (islogical(subAxis) && subAxis)
  imshow(array, [1 1 1; 0 0 0; colors],'InitialMagnification', 'fit');
elseif ischar(subAxis)
  imshow(array, [1 1 1; 0 0 0; colors], 'InitialMagnification', subAxis);
elseif islogical(subAxis) && (~subAxis)
  imshow(array, [1 1 1; 0 0 0; colors]);
else
  imshow(array, [1 1 1; 0 0 0; colors]);
  axis(subAxis);
end
%set(gca, 'YDir', 'normal');

if (tiempo==-1)     % esperar tecla
    pause
elseif tiempo>0
    pause(tiempo);
end
