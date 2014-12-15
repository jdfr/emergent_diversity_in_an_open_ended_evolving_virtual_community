function sim = similitud(varargin)

if ischar(varargin{1})
  [raster1, offset1] = ls2(varargin{1}, 'canonical');
  [raster2, offset2] = ls2(varargin{2}, 'canonical');
else
  raster1 = varargin{1};
  offset1 = varargin{2};
  raster2 = varargin{3};
  offset2 = varargin{4};
end

pixel1 = [raster1{:}];
pixel2 = [raster2{:}];

pixel1(:,1) = pixel1(:,1) - offset1(1);
pixel1(:,2) = pixel1(:,2) - offset1(2);

pixel2(:,1) = pixel2(:,1) - offset2(1);
pixel2(:,2) = pixel2(:,2) - offset2(2);

if isempty(pixel1) && isempty(pixel2)
    sim = 1;
elseif (~isempty(pixel1) && isempty(pixel2)) || (isempty(pixel1) && ~isempty(pixel2))
    sim = 0;
else
    sim = similitudFast(pixel1, pixel2);
end

  
