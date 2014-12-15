function [Z sim numGroups raster individualsSpecies dists] = testClustering(P)
%this was part of the effor to convince myself that my code had no bugs,
%in spite of Rene and Paco's opinion, when implementing Rene's way to
%determine the number of clusters. I WAS RIGHT  

npixsoil = P.counts(:,3);
raster   = P.raster(npixsoil==0);
off      = P.offset(npixsoil==0,:);

dists = [];

k=1;
for a=1:numel(raster)-1
  for b=a+1:numel(raster)
    dists(k) = 1 - similitud2(raster{a},off(a,:),raster{b},off(b,:));
    k=k+1;
    %dists(b,a) = dists(a,b);
  end
end

Z = linkage(dists);
dists = squareform(dists); 

    sZ = size(Z,1);
    c = 1;
    simi = 1;
    while c<=sZ+1 && simi >0

        % groups = cluster(Z,'cutoff',c);
        groups = cluster(Z,'maxclust',c,'criterion','distance');
        for i = 1 : c
            eachgroup{c}{i}= find(groups == i);
            disrups=zeros(numel(eachgroup{c}{i}),1);
            for w=1:numel(eachgroup{c}{i})
              disrups(w) = mean(dists(eachgroup{c}{i}(w), eachgroup{c}{i}));
            end
            [meanDisrup{c}{i}]= min(disrups);
            individualsSpecies{c}{i} = eachgroup{c}{i};
        end


        sim(c) = sum([meanDisrup{c}{:}])/c;
        simi = sim(c);

        c = c + 1;
    end
numGroups = calculateNumGroups(sim, struct('oldWayToFindPrototype', false, 'thresholdMode', 'absolute',   'threshold', 0.01));



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

