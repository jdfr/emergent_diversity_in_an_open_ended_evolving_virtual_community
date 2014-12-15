function [colors array]= drawtreeleafed(raster, leaves, dimensions,offset,tiempo,indexInColors,colors, subAxis)


if isempty(raster)
  return
end

if ~exist('colors', 'var')
    colors = rand(1,3)*0.8;
    indexInColors = 1;
end

ys = raster{1}; if size(ys,2)>1; ys = ys(:); end;
xs = raster{2}; if size(xs,2)>1; xs = xs(:); end;
array = accumarray([ys,xs],uint16(indexInColors+2), dimensions);
array(sub2ind(dimensions, leaves{1}, leaves{2})) = array(sub2ind(dimensions, leaves{1}, leaves{2})) + 1;
yOffset = offset(1);
if yOffset>1;
  array(yOffset-1,:) = 2;
end
array = array(end:-1:1,:);
if size(array,2)==1
  array(end,2)=0;
end
if (~exist('subAxis', 'var')) || isempty(subAxis) || ~any(subAxis)
  imshow(array, [1 1 1; 0 0 0; colors],'InitialMagnification', 'fit');
elseif ischar(subAxis)
  imshow(array, [1 1 1; 0 0 0; colors], 'InitialMagnification', subAxis);
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
