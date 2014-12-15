function [especies comEspeciesAnte comEspeciesSig]= relacionarSpecies(arg1,numGenerations)

if ischar(arg1)
  basedir = arg1;
  poph = loadSimulation(basedir);
elseif iscell(arg1)
  poph = arg1{1};
  basedir = arg1{2};
else
  error('arg1 has incorrect type!!!!');
end

numIndividuals = zeros(1,numGenerations+1);
for i = 1 : numGenerations+1
    numIndividuals(i) = sum(poph.generation == i-1); % numIndividuals(1) = numero de individuos en la generacion 0
end

indicesPredecesores = poph.tree.iA;

for i = 1 : numGenerations-1
    i
    predecesores = indicesPredecesores(sum(numIndividuals(1:i+1))+1:sum(numIndividuals(1:i+2)));

    load(sprintf('%s/numGroups%03d.mat',basedir,i));
    numGroupsA = eval(sprintf('numGroups%03d', i));

    load(sprintf('%s/individualsSpecies%03d.mat',basedir,i));
    individualsSpeciesA = eval(sprintf('individualsSpecies%03d', i));

    load(sprintf('%s/comLZfigureLongitud%03d.mat',basedir,i));
    comLZfigureLongitudA = eval(sprintf('comLZfigureLongitud%03d', i));

    load(sprintf('%s/numGroups%03d.mat',basedir,i+1));
    numGroupsS = eval(sprintf('numGroups%03d', i+1));

    load(sprintf('%s/individualsSpecies%03d.mat',basedir,i+1));
    individualsSpeciesS = eval(sprintf('individualsSpecies%03d', i+1));

    load(sprintf('%s/comLZfigureLongitud%03d.mat',basedir,i+1));
    comLZfigureLongitudS = eval(sprintf('comLZfigureLongitud%03d', i+1));

    for j = 1 : numGroupsA
        longEspecie = length(individualsSpeciesA{numGroupsA}{j});
        cont = cell(longEspecie,1);
        cont2 = cell(longEspecie,1);
        for k = 1 : longEspecie
            t = 1;
            f = find(predecesores == individualsSpeciesA{numGroupsA}{j}(k));
            for q = 1 : length(f)
                for p = 1 : numGroupsS
                    if ~isempty(find(individualsSpeciesS{numGroupsS}{p}==f(q)))
                        cont{k}{t} = p;
                        t = t + 1;
                    end
                end
            end
            if ~isempty(cont{k})
                cont2{k}=unique([cont{k}{:}]);
            else
                cont2{k}=[];
            end
        end
        especies{i}{j}=unique([cont2{:}]);
        comEspeciesSig{i}{j} = comLZfigureLongitudS(especies{i}{j});
        comEspeciesAnte{i}{j} = comLZfigureLongitudA(j);
    end
end

tosave1 = ['especies' num2str(numGenerations, '%03g')];
tosave2 = ['comEspeciesSig' num2str(numGenerations, '%03g')];
tosave3 = ['comEspeciesAnte' num2str(numGenerations, '%03g')];

eval([tosave1 ' = especies;']);
eval([tosave2 ' = comEspeciesSig;']);
eval([tosave3 ' = comEspeciesAnte;']);

save([basedir filesep tosave1],tosave1);
save ([basedir filesep tosave2],tosave2);
save ([basedir filesep tosave3],tosave3);

h = figure('Visible', 'off');

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
            plot([i-1 i],[comEspeciesAnte{i}{j} comEspeciesSig{i}{j}],'color',r(j,:))
            hold on
            colores(g,:) = r(j,:);
            g = g+1;
            %             pause
        elseif length(comEspeciesSig{i}{j}) > 1

            total = length(comEspeciesSig{i}{j});
            t = rand(total,3);
            for k = 1 : total
                plot([i-1 i],[comEspeciesAnte{i}{j} comEspeciesSig{i}{j}(k)],'color',t(k,:))
                hold on
                colores(g,:) = t(k,:);
                g = g+1;
                %                 pause
            end

        end
    end
    for j = 1 : length([comEspeciesAnte{i+1}{:}])
        a = [comEspeciesAnte{i+1}{:}];
        fe = find([comEspeciesSig{i}{:}]==a(j));
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

[fit1 fit2]=paintComplexityProto2(basedir,1,numGenerations);
plot(1:numGenerations,polyval(fit1,1:numGenerations),'k','LineWidth',2);
hold on
plot(1:numGenerations,polyval(fit2,1:numGenerations),'k','LineWidth',2);

% [x1 pol1 pol2]=paintComplexityProto2(basedir,1,numGenerations);
% plot(x1,pol1,'k','LineWidth',2);
% hold on
% plot(x1,pol2,'k','LineWidth',2);

xlabel('generations')
ylabel('evolution of species complexity')

saveas(h, [basedir filesep 'relacionarSpecies.fig'], 'fig');

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

function [fit1 fit2]=paintComplexityProto2(basedir,GenI,GenF)

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

    if ~isempty(find(T{i}==pos(1)))

        stdt(i,1) = max(comLZfigureLongitud(find(T{i}==pos(1))));
    else
        stdt(i,1) = 0;
    end

    if ~isempty(find(T{i}==pos(2)))

        stdt(i,2) = max(comLZfigureLongitud(find(T{i}==pos(2))));
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

