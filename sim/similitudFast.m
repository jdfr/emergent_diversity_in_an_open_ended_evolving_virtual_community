function sim = similitudFast(pixel1, pixel2)


% sim = size(intersect(pixel1,pixel2,'rows'),1)/size(union(pixel1,pixel2,'rows'),1);

%this is faster than the line above, but only works if the arrays of pixels
%have no repeated rows (that is to say, if the arrays are row sets)

sim             = 0; %#ok<NASGU>

%sort rows of both trees
sortedrows      = [pixel1; pixel2];
[ind,ind]    = sort(sortedrows(:,  2)); %#ok<ASGLU>
[ind2,ind2]   = sort(sortedrows(ind,1)); %#ok<ASGLU>
ind             = ind(ind2);
sortedrows      = sortedrows(ind,:);
%find matching entries
matching        = all(sortedrows(1:end-1,:) == sortedrows(2:end,:), 2);
sizeintersec    = sum(matching);
sim             = sizeintersec/(size(pixel1,1)+size(pixel2,1)-sizeintersec);

