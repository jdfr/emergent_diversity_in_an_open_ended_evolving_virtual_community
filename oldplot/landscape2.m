function f=landscape2(ramas,xmax,xmin,ymax,ymin,nomdir,Generation)

xmax = xmax*100;
xmin = xmin*100;
ymax = ymax*100;
ymin = ymin*100;

% axis([xmin xmax ymin ymax])
sR = size(ramas,2);
% maxyoy1 = zeros(1,sR);
nramas = zeros(1,sR);
N = floor((sR*(xmax+abs(xmin)))/20);
values = ymin*ones(1,N);
label = zeros(1,N);
h = figure(1);
clf
% h = figure('visible','off');
for i = 1 : sR
    i
    if ~isempty(ramas{i})
        r = ramas{i}(:,1:4);
        if ~isempty(r)
            r = floor(r*100);
%             maxyoy1(i) = max(max(r(:,3:4)));
            minx0x1 = min(min(r(:,1:2)));
            maxx0x1 = max(max(r(:,1:2)));
            x = N + 1 - maxx0x1;
            while x + maxx0x1 > N
                randN = randint(1,1,[1,N]); % situo el arbol en algun natural n >= 1
                x = randN - minx0x1;
            end
            r(:,1:2) = r(:,1:2)+x;
            t = r;
            t(:,5) = ramas{i}(:,5);
            ran = [N,0,ymax,ymin];
            drawtree(t,ran,1);
            hold on
            sr = size(r,1);
            nramas(i) = sr;
            for j = 1 : sr
                x0 = r(j,1);
                x1 = r(j,2);
                y0 = r(j,3);
                y1 = r(j,4);
                if x0 == x1
                    maxy0y1 = max(y0,y1);
                    if values(x0) == maxy0y1 && label(x0) ~= i
                        values(x0) = ymin;
                        label(x0) = -1;
                    elseif maxy0y1 > values(x0)
                        values(x0) = maxy0y1;
                        label(x0) = i;
                    end
                elseif x0<x1
                    for k = x0:x1
                        y = y0 + (((k-x0)*(y1-y0))/(x1-x0));
                        if values(k) == y && label(k) ~= i
                            values(k) = ymin;
                            label(k) = -1;
                        elseif y > values(k)
                            values(k) = y;
                            label(k) = i;
                        end

                    end
                else
                    for k = x1:x0
                        y = y0 + (((k-x0)*(y1-y0))/(x1-x0));
                        if values(k) == y && label(k) ~= i
                            values(k) = ymin;
                            label(k) = -1;
                        elseif y > values(k)
                            values(k) = y;
                            label(k) = i;
                        end

                    end

                end
            end
        end
    end
end

saveas(h,[nomdir,'entorno' num2str(Generation, '%03g')],'png')
clf

label=label(label ~= -1); % quitar los -1 (donde se pisan los grafos no tienen luz)
label=label(label ~= 0); % quitar los 0 (las x que no tienen y)
% v = sort(label)';
v = label;
z = zeros(1,sR);

z(v(find([1 diff(v)])))=[diff(find([1 diff(v)])) length(find(v==v(end)))];

% dividiendo por el numero de ramas distinguibles del arbol *|G|
for i = 1 : sR
    if nramas(i) > 0
        f(i)=z(i)/(nramas(i)*10); % |G|= 10
    else
        f(i)= 0;
    end
end



% dividiendo por la longitud del arbol
% for i = 1 : sR
%     if maxyoy1(i) > 0
%         f(i)=z(i)./maxyoy1(i);
%     else
%         f(i)= 0;
%     end
% end




