function showDetailedImage(imgname, imgwrite, sep, width, combs, doReplaceSky, showUpper)
%combs = {init1 end1 height1 col1 width1; init2 end2 height2 col2 width2;...}

[img map]= imread(imgname);

[img map]      = imread(imgname);
if ~isempty(map)
  img          = myind2rgb(img, map);
end

inits   = cell2mat(combs(:,1));
ends    = cell2mat(combs(:,2));
heights = cell2mat(combs(:,3));
cols    = combs(:,4);
widths  = cell2mat(combs(:,5));

widths(isnan(widths))=width;
heights(isnan(heights))=1;%size(img,1); %the heights are calculated from above (that is to say, from the first row)

maxw = max([width;widths(:)]);

if doReplaceSky
  for k=1:numel(inits)
    a = inits(k);
    b = ends(k);
    img(:,a:b,:) = replaceSky(img(:,a:b,:), cols{k});
    img(:,a,:) = 0;
    img(:,b,:) = 0;
  end
end

imgwidth = size(img,2);

factor = imgwidth/width;

resizeParams      = {'Antialiasing', true, 'Method', 'lanczos2'};

if showUpper
  imgw = imresize(img, 'OutputSize', [nan width], resizeParams{:});

  if size(imgw,2)~=maxw
    imgw = putOnSheet(imgw, makeBackground(imgw, maxw));
  end
end

imgs    = cell(size(combs, 1),1);

for k=1:numel(imgs)
  a = inits(k);
  b = ends(k);
  h = heights(k);
  w = widths(k);
  subimg = img(h:end,a:b,:);
  imgs{k} = imresize(subimg, 'OutputSize', [nan w], resizeParams{:});
  if size(imgs{k},2)~=maxw
    imgs{k} = putOnSheet(imgs{k}, makeBackground(imgs{k}, maxw));
  end
end

if showUpper
  sep = zeros([sep, maxw,3],  class(imgw));

  newimg = vertcat(imgw, sep, imgs{:});
else
  newimg = vertcat(imgs{:});
end

imwrite(newimg, imgwrite, 'png');

function putOnSheet(img, sheet)
  sizesheet = size(sheet);
  pxthis            = sizesheet(1)-size(img,1)+1;
  pythis            = floor((sizesheet(2)-size(img,2))/2+1);
  frame             = copyOn(sheet, [pxthis pythis], img);

function back = makeBackground(img, maxw)
 if isfloat(img)
   back = ones(size(img,1), maxw);
 else
   back = cast(inf, class(img));
   back = back(ones(size(img,1),1),  ones(maxw,1));
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
