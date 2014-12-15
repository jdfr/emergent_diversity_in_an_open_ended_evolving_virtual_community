function h = showEvolutionOfComplexity(basedir,numGenerations, comEspeciesAnte, comEspeciesSig)

h = figure('Visible', 'on');

colores = rand(1,3);
for i = 1 : numGenerations-2
    if i == 173
        disp('hola')
    end

    puntos = size(comEspeciesSig{i},2);
    g = 1;
    r = colores;
    colores = [];
    colores2 = [];
    t = [];
    for j = 1 : puntos

        if length(comEspeciesSig{i}{j}) == 1
            line([i-1 i],[comEspeciesAnte{i}{j} comEspeciesSig{i}{j}],'color',r(j,:))
            hold on
            colores(g,:) = r(j,:);
            g = g+1;
            %             pause
        elseif length(comEspeciesSig{i}{j}) > 1

            total = length(comEspeciesSig{i}{j});
            t = rand(total,3);
            for k = 1 : total
                line([i-1 i],[comEspeciesAnte{i}{j} comEspeciesSig{i}{j}(k)],'color',t(k,:))
                hold on
                colores(g,:) = t(k,:);
                g = g+1;
                %                 pause
            end

        end
    end
    for j = 1 : length([comEspeciesAnte{i+1}{:}])
        a = [comEspeciesAnte{i+1}{:}];
%         fe = find([comEspeciesSig{i}{:}]==a(j));
try
        fe = find(vertcat(comEspeciesSig{i}{:})==a(j));
catch
        fe = find(horzcat(comEspeciesSig{i}{:})==a(j));
end
        colores2(j,:) = colores(fe(1),:);
    end
    colores = colores2;
    %     pause
end
hold on
% [x1 y1 x21 y21 r2]= paintComplexityProto(basedir,numGenerations);
% plot([x1 x21],[y1 y21],'k','LineWidth',2)
% hold on
% plot([x1 numGenerations],[y1 r2],'k','LineWidth',2)

[stdt fit1 fit2 fit3]=paintComplexityProto2(basedir,1,numGenerations);
hold on
% line(1:numGenerations,stdt(:,1)','color', 'r','LineWidth',2);
% line(1:numGenerations,stdt(:,2)','color', 'g','LineWidth',2);
% line(1:numGenerations,stdt(:,3)','color', 'b','LineWidth',2);
line(1:numGenerations,polyval(fit1,1:numGenerations),'color', 'k','LineWidth',2);
line(1:numGenerations,polyval(fit2,1:numGenerations),'color', 'k','LineWidth',2);
% line(1:numGenerations,polyval(fit3,1:numGenerations),'color', 'k','LineWidth',2);

% [x1 pol1 pol2]=paintComplexityProto2(basedir,1,numGenerations);
% plot(x1,pol1,'k','LineWidth',2);
% hold on
% plot(x1,pol2,'k','LineWidth',2);

xlabel('generations')
ylabel('evolution of species complexity')

%saveas(h, [basedir filesep 'relacionarSpecies.fig'], 'fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x1 y1 x21 y21 r2]=paintComplexityProto(basedir,numGenerations)

for i = 1 : numGenerations
    load(sprintf('%s/comLZfigureLongitud%03d.mat',basedir,i));
    comLZfigureLongitud = eval(sprintf('comLZfigureLongitud%03d', i));
    desviacionTipica(i) = std(comLZfigureLongitud);
end

mediaDes = mean(desviacionTipica);
[maxDes pos]= max(desviacionTipica);
x1 = 0;
y1 = 0;
x21 = i;
x22 = pos;
y21 = mediaDes;
y22 = maxDes;
m2 = y22/x22;
r2 = numGenerations*m2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stdt fit1 fit2 fit3]=paintComplexityProto2(basedir,GenI,GenF)

stdt = zeros(GenF,3);
for i = GenI : GenF
    if i == 33
        disp('hola')
    end
    load(sprintf('%s/comLZfigureLongitud%03d.mat',basedir,i));
    comLZfigureLongitud = eval(sprintf('comLZfigureLongitud%03d', i));
    if length(comLZfigureLongitud) == 1
        T{i}=1;

    else
        Y = distancia(comLZfigureLongitud);
        Z = linkage(Y,'single');
        T{i} = cluster(Z,'maxclust',3);
    end
    a = comLZfigureLongitud(find(T{i}==1));
    mediaa = mean(a);
    b = comLZfigureLongitud(find(T{i}==2));
    mediab = mean(b);
    c = comLZfigureLongitud(find(T{i}==3));
    mediac = mean(c);

    medias = [mediaa mediab mediac];
    [ordenar pos] = sort(medias);

    %if ~isempty(find(T{i}==pos(1)))
    if ~isempty(find(T{i}==pos(2)))

        stdt(i,1) = max(comLZfigureLongitud(find(T{i}==pos(1))));
        %stdt(i,1) = min(comLZfigureLongitud(find(T{i}==pos(2))));
    else
        stdt(i,1) = 0;
    end

    %if ~isempty(find(T{i}==pos(2)))
    if ~isempty(find(T{i}==pos(3)))

        stdt(i,2) = max(comLZfigureLongitud(find(T{i}==pos(2))));
        %stdt(i,2) = min(comLZfigureLongitud(find(T{i}==pos(3))));
    else
        stdt(i,2) = 0;
    end
    if ~isempty(find(T{i}==pos(3)))

        stdt(i,3) = max(comLZfigureLongitud(find(T{i}==pos(3))));
    else
        stdt(i,3) = 0;
    end
end

fit1 = polyfit(1:GenF,stdt(:,1)',1);
fit2 = polyfit(1:GenF,stdt(:,2)',1);
fit3 = polyfit(1:GenF,stdt(:,3)',1);


y1 = stdt(:,1)';
x1 = 1:GenF;
Y1 = log(y1);
M1 = [length(x1),sum(x1);sum(x1),sum(x1.^2)];
u1 = [sum(Y1);sum(x1.*Y1)];
sol1 = M1\u1;
a1 = exp(sol1(1));
b1 = sol1(2);

pol1 = a1*exp(b1*x1);

y2 = stdt(:,2)';
noceros = find(y2~=0);
y2 = y2(noceros);
x2 = x1(noceros);
Y2 = log(y2);
M2 = [length(x2),sum(x2);sum(x2),sum(x2.^2)];
u2 = [sum(Y2);sum(x2.*Y2)];
sol2 = M2\u2;
a2 = exp(sol2(1));
b2 = sol2(2);

pol2 = a2*exp(b2*x1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dist=distancia(X)

long = length(X);
k = 1;
for i = 1 : long-1
    for j = (i+1): long
        dist(k) = abs(X(j) - X(i));
        k = k + 1;
    end
end

