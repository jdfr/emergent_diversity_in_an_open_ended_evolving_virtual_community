function showDetailedImageByLines(imgname, imgwrite, scale, finalwidth, heights)
%variantion over showDetailedImage where, instead of adjusting the scale to
%fit all lines to the same length, the scale is established beforehand, and
%the last line might have a different length, as Rene wished
%combs = {init1 end1 height1 col1 width1; init2 end2 height2 col2 width2;...}

[img map]      = imread(imgname);
if ~isempty(map)
  img          = myind2rgb(img, map);
end

origlen = size(img,2);
finallen = origlen*scale;
origwidth = round(finalwidth/scale);

heights(isnan(heights))=1;%size(img,1); %the heights are calculated from above (that is to say, from the first row)
resizeParams      = {'Antialiasing', true, 'Method', 'lanczos2'};
pointer = 1;
k=1;
imgs = {};
while pointer<=origlen
  sect = pointer:min(origlen, pointer+origwidth-1);
  h = heights(min(k, numel(heights)));
  subimg = img(h:end,sect,:);
  imgs{k} = imresize(subimg, 'OutputSize', [nan numel(sect)*scale], resizeParams{:}); %#ok<AGROW>
  if size(imgs{k},2)<finalwidth
    imgs{k} = [imgs{k}, makeBackground(imgs{k}, finalwidth-size(imgs{k},2))]; %#ok<AGROW>
  elseif size(imgs{k},2)==(finalwidth+1)
    imgs{k} = imgs{k}(:,1:end-1,:);
  else
    error('jarlll111111!!!!!!!!');
  end
  if size(imgs{k},2)~=finalwidth
    error('jarlll!!!!!!!!');
  end
  k=k+1;
  pointer = sect(end)+1;
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
