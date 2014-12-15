function leaf = leafTree(ramas)

sr = size(ramas,1);
k=1;
for i = 1 : sr
    r = ramas(i,:);
    x1 = r(2);
    y1 = r(4);
    x0equalx1 = find(ramas(:,1)==x1);
    if isempty(x0equalx1)
        leaf(k,:) = [x1 y1];
        k = k + 1;
    else
        y0equaly1 = find(ramas(x0equalx1,3)==y1);
        if isempty(y0equaly1)
            leaf(k,:) = [x1 y1];
            k = k + 1;
        end
    end
end
