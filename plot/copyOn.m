function newi = copyOn(background, pos, image, alphaColor)

newi = background;

if exist('alphaColor', 'var')
  if numel(size(image))>2
    [alphax alphay] = find((image(:,:,1)~=alphaColor(1)) & ...
                           (image(:,:,2)~=alphaColor(2)) & ...
                           (image(:,:,3)~=alphaColor(3))  );
  else
    [alphax alphay] = find((image(:,:,1)~=alphaColor(1)) & ...
                           (image(:,:,2)~=alphaColor(2)) );
  end

  for k=1:numel(size(image))
    indexesb = sub2ind(size(background), alphax+pos(1)-1, alphay+pos(2)-1, repmat(k, size(alphax)));
    indexesi = sub2ind(size(image),      alphax,        alphay,        repmat(k, size(alphax)));

    newi(indexesb) = image(indexesi);
  end
else
  if pos(1)+size(image,1)-1 > size(background, 1)
    error('background has too few rows!!!');
  end
  if pos(2)+size(image,2)-1 > size(background, 2)
    error('background has too few columns!!!');
  end
  newi(pos(1):pos(1)+size(image,1)-1, pos(2):pos(2)+size(image,2)-1, :) = image;
end

%figure; imshow(newi);
