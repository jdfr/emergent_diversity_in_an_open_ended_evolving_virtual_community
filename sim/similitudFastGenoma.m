function sim = similitudFastGenoma(genoma1,genoma2,mode)


[raster1, offset1, counts1, dimension1, isLeaf1, maxiy1, leaves1]=ls2(genoma1, mode);
[raster2, offset2, counts2, dimension2, isLeaf2, maxiy2, leaves2]=ls2(genoma2, mode);

pixel1 = [raster1{:}];
pixel2 = [raster2{:}];

pixel1(:,1) = pixel1(:,1) - offset1(1);
pixel1(:,2) = pixel1(:,2) - offset1(2);

pixel2(:,1) = pixel2(:,1) - offset2(1);
pixel2(:,2) = pixel2(:,2) - offset2(2);


% sim = size(intersect(pixel1,pixel2,'rows'),1)/size(union(pixel1,pixel2,'rows'),1);

%this is faster than the line above, but only works if the arrays of pixels
%have no repeated rows (that is to say, if the arrays are row sets)

sim             = 0; %#ok<NASGU>

%sort rows of both trees
sortedrows      = [pixel1; pixel2];
[ignore,ind]    = sort(sortedrows(:,  2)); %#ok<ASGLU>
[ignore,ind2]   = sort(sortedrows(ind,1)); %#ok<ASGLU>
ind             = ind(ind2);
sortedrows      = sortedrows(ind,:);
%find matching entries
matching        = all(sortedrows(1:end-1,:) == sortedrows(2:end,:), 2);
sizeintersec    = sum(matching);
sim             = sizeintersec/(size(pixel1,1)+size(pixel2,1)-sizeintersec);

