function stackEnvironmentsTrueScale(imgnames, imgwrite, sep)
%take several images, and stack them, showing them true scale
%combs = {init1 end1 height1 col1 width1; init2 end2 height2 col2 width2;...}

[imgs maps]= cellfunc(@imread, imgnames);
for k=1:numel(maps)
  if ~isempty(maps{k})
    imgs{k}          = myind2rgb(imgs{k}, maps{k});
  end
end

widths = cellfun(@(x)size(x,2), imgs);

maxw = max(widths);

sep = ones([sep, maxw,3],  class(imgs{1}));
if not(isfloat(sep))
  sep(1:end) = inf;
end

for k=1:numel(imgs)
  imgs{k} = putOnSheet(imgs{k}, makeBackground(imgs{k}, maxw));
end
for k=1:numel(imgs)-1
  imgs{k} = [imgs{k}; sep];
end

newimg = vertcat(imgs{:});

imwrite(newimg, imgwrite, 'png');

function frame = putOnSheet(img, sheet)
  sizesheet = size(sheet);
  pxthis            = sizesheet(1)-size(img,1)+1;
  pythis            = floor((sizesheet(2)-size(img,2))/2+1);
  frame             = copyOn(sheet, [pxthis pythis], img);

function back = makeBackground(img, maxw)
 if isfloat(img)
   back = ones(size(img,1), maxw, 3);
 else
   back = cast(inf, class(img));
   back = back(ones(size(img,1),1),  ones(maxw,1), [1 1 1]);
 end

function thisim            = replaceSky(thisim, imageSkyColor)
  u255 = uint8(255);
  for k=1:size(thisim,1)
    l = thisim(k,:,:);
    areSky = find(all(l == u255, 3));
    nl = numel(l)/3;
    for z=1:3
      l(areSky) = imageSkyColor(z);
      areSky = areSky+nl;
    end
    thisim(k,:,:) = l;
  end

function o = myind2rgb(a,cm)
  if ~isa(a, 'double')
      a = a+1;    % Switch to one based indexing
  end

  if isa(cm, 'double')
      cm = uint8(cm*255);
  elseif isa(cm, 'uint8')
  else
    error('unexpected!!!!!!');
  end

  zz = double(max(max(a)));
  zb = size(cm,1);
  if zz>zb
    cm = [cm; cm(10:(9+zz-zb),:)];
  end

  % % Make sure A is in the range from 1 to size(cm,1)
  % a = max(1,min(a,size(cm,1)));


  r = zeros(size(a), 'uint8'); r(:) = cm(a,1);
  g = zeros(size(a), 'uint8'); g(:) = cm(a,2);
  b = zeros(size(a), 'uint8'); b(:) = cm(a,3);

  o = cat(3, r, g, b);
