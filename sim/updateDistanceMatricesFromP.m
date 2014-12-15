function [dists indicesClasificados] = updateDistanceMatricesFromP(P, dists, mutated, changed, update)

if isempty(update)
  dists = cell(size(update));
  return
end

if isempty(dists)
  indicesClasificados = [];
  return
end

gnm = P.genome;
mng = P.mingnm;
raster = P.raster;
offset = P.offset;
nPixelsInSoil = P.counts(:,3);
negativos = nPixelsInSoil > 0;

%indicesNoClasificados = find(negativos);
indicesClasificados = find(not(negativos));

ipht = find(strcmp('pht', update));
ignm = find(strcmp('gnm', update));
imng = find(strcmp('mng', update));
dopht = not(isempty(ipht));
dognm = not(isempty(ignm));
domng = not(isempty(imng));
%case12 = dognm && not(domng);

%do them in this order
idxs = [imng ipht];

%first, mutated individuals have their genome distances re-evaluated
if dognm
  for k=1:numel(mutated)
    if mod(k,100)==0
      fprintf('matrixdist mutated %d/%d\n', k,numel(mutated));
    end
    c = mutated(k);
    if negativos(c) %no need to calculate anything for plants with underground branches
      continue
    end
    for j=1:numel(indicesClasificados) %no need to calculate anything for plants with underground branches
      i = indicesClasificados(j);
      v = strdstMX(gnm{c}, gnm{i});
      dists{ignm}(c,i) = v;
      dists{ignm}(i,c) = v;
    end
  end
end

%after that, we calculate red.genome distances and phenotypic distances for
%individuals that might have changed
for k=1:numel(changed)
  if mod(k,100)==0
    fprintf('matrixdist changed %d/%d\n', k,numel(changed));
  end
  c = changed(k);
  if negativos(c) %no need to calculate anything for plants with underground branches
    continue
  end
  if dopht
    raster2 = raster{c};
    offset2 = offset(c,:);
  end
  for j=1:numel(indicesClasificados) %no need to calculate anything for plants with underground branches
    i = indicesClasificados(j);
    for mm = 1:numel(idxs)
      m = idxs(mm);
      switch update{m}(1)
        case 'm' %mng
          v = strdstMX(mng{c}, mng{i});
        case 'p' %pht
          if (domng && (dists{imng}(i,c)==0)) %|| (case12 && (dists{ignm}(i,c)==0))
            v = 0;
          else
            raster1 = raster{i};
            offset1 = offset(i,:);
            v = 1 - similitud2(raster1,offset1,raster2,offset2);
          end
      end
      dists{m}(c,i) = v;
      dists{m}(i,c) = v;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sim = similitud2(raster1,offset1,raster2,offset2)%, show)
    pixel1 = [raster1{:}];
    pixel2 = [raster2{:}];
    v = isempty(pixel1) + isempty(pixel2);
if v==2
    sim = 1;
elseif v==1
    sim = 0;
else
    
    if size(pixel1,1) == 1 
        pixel1 = [];

        pixel1(:,1) = raster1{1};
        pixel1(:,2) = raster1{2};
    end
    if size(pixel2,1) == 1 
        pixel2 = [];

        pixel2(:,1) = raster2{1};
        pixel2(:,2) = raster2{2};
    end
    
    pixel1(:,1) = pixel1(:,1) - offset1(1);
    pixel1(:,2) = pixel1(:,2) - offset1(2);

    pixel2(:,1) = pixel2(:,1) - offset2(1);
    pixel2(:,2) = pixel2(:,2) - offset2(2);
    
    sim = similitudFast(pixel1, pixel2);
    %sim = size(intersect(pixel1,pixel2,'rows'),1)/size(union(pixel1,pixel2,'rows'),1);
end

% if isempty([raster1{:}]) && isempty([raster2{:}])
%     sim = 1;
% elseif (isempty([raster1{:}]) && ~isempty([raster2{:}])) || (~isempty([raster1{:}]) && isempty([raster2{:}]))
%     sim = 0;
% else
%     pixel1 = [raster1{:}];
%     pixel2 = [raster2{:}];
% 
%     pixel1(:,1) = pixel1(:,1) - offset1(1);
%     pixel1(:,2) = pixel1(:,2) - offset1(2);
% 
%     pixel2(:,1) = pixel2(:,1) - offset2(1);
%     pixel2(:,2) = pixel2(:,2) - offset2(2);
% 
%     %     sraster1 = size(pixel1,1);
%     %     sraster2 = size(pixel2,1);
% 
%     sim = size(intersect(pixel1,pixel2,'rows'),1)/size(union(pixel1,pixel2,'rows'),1);
% end


