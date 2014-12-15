function [dists indicesClasificados] = generateDistanceMatricesFromP(P, generate)

if isempty(generate)
  dists = cell(size(generate));
  return
end

gnm = P.genome;
mng = P.mingnm;
raster = P.raster;
offset = P.offset;
nPixelsInSoil = P.counts(:,3);
ramasNegativas = nPixelsInSoil;
numIndividual = numel(gnm);

negativos = find(ramasNegativas > 0);

indicesNoClasificados = negativos;
indicesClasificados = setdiff(1:numIndividual,indicesNoClasificados);
numIndividual = length(indicesClasificados);

gnm(indicesNoClasificados,:)=[];
mng(indicesNoClasificados,:)=[];
raster(indicesNoClasificados,:)=[];
offset(indicesNoClasificados,:)=[];

dists = repmat( {zeros(1,(numIndividual*(numIndividual-1))/2)}, size(generate) );
ipht = find(strcmp('pht', generate));
ignm = find(strcmp('gnm', generate));
imng = find(strcmp('mng', generate));
dopht = not(isempty(ipht));
dognm = not(isempty(ignm));
domng = not(isempty(imng));
case12 = dognm && not(domng);

%do them in this order
idxs = [ignm imng ipht];

k = 1;
if numIndividual > 1
  for i = 1 : numIndividual-1
    if mod(i,100)==0
      fprintf('matrixdist %d/%d\n', i,numIndividual);
    end
    if dopht
      raster2 = raster{i};
      offset2 = offset(i,:);
    end
    for j = i+1 : numIndividual
      for mm = 1:numel(idxs)
        m = idxs(mm);
        switch generate{m}(1)
          case 'g' %gnm
            dists{m}(k) = strdstMX(gnm{i}, gnm{j});
          case 'm' %mng
            if dognm && (dists{ignm}(k)==0)
              dists{m}(k) = 0;
            else
              dists{m}(k) = strdstMX(mng{i}, mng{j});
            end
          case 'p' %pht
            if (domng && (dists{imng}(k)==0)) || (case12 && (dists{ignm}(k)==0))
              dists{m}(k) = 0;
            else
              raster1 = raster{j};
              offset1 = offset(j,:);
              dists{m}(k) = 1 - similitud2(raster1,offset1,raster2,offset2);
            end
        end
      end
      k = k + 1;
    end
  end
  dists = cellfunc(@squareform, dists);
elseif numIndividual==1
  for m = 1:numel(generate)
    dists{m} = 0;
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


