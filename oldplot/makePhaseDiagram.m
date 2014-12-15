function makePhaseDiagram(basedirs)

C1 = Gray;
C2 = Spring;
C3 = Summer;
color = [C1(1,:);C1(20,:);C1(40,:);C1(60,:);C2(1,:);C2(10,:);C2(20,:);C2(30,:);C3(30,:);C3(40,:);C3(50,:);C3(60,:)];
numDirect = size(basedirs,2);

for k = 1 : numDirect
    varianceAltura = [];
    varianceArea = [];
    varianceRamas = [];
    varianceHojas = [];
    varianceFit = [];
    
    altura = [];
    anchura = [];
    numleafsOnTop = [];
    numRamas = [];
    numleafs = [];
    fitness = [];

    
    poph = loadSimulation(basedirs{k});

    altura = poph.height;
    anchura = poph.width;
    numleafsOnTop = poph.nleafsOnTop;
    numRamas = poph.nbranches;
    numleafs = poph.nleafs;
    fitness = poph.fitness;

    numGenerations = max(poph.generation);

    numIndividuals = zeros(1,numGenerations+1);
    for j = 1 : numGenerations+1
        numIndividuals(j) = sum(poph.generation == j-1); % numIndividuals(1) = numero de individuos en la generacion 0
    end

    varianceAltura(1,:) = var(altura(1:numIndividuals(1)));
    varianceArea(1,:) = var(altura(1:numIndividuals(1)).*anchura(1:numIndividuals(1)));
     varianceRamas(1,:) = var(numRamas(1:numIndividuals(1)));
    varianceHojas(1,:) = var(numleafsOnTop(1:numIndividuals(1))./numleafs(1:numIndividuals(1)));
    varianceFit(1,:) = mean(fitness(1:numIndividuals(1)));

    for i = 2 : numGenerations+1
        varianceAltura(i,:) = var(altura(sum(numIndividuals(1:i-1))+1:sum(numIndividuals(1:i))));
        varianceArea(i,:) = var(altura(sum(numIndividuals(1:i-1))+1:sum(numIndividuals(1:i))).*anchura(sum(numIndividuals(1:i-1))+1:sum(numIndividuals(1:i))));
        varianceRamas(i,:) = var(numRamas(sum(numIndividuals(1:i-1))+1:sum(numIndividuals(1:i))));
        varianceHojas(i,:) = var(numleafsOnTop(sum(numIndividuals(1:i-1))+1:sum(numIndividuals(1:i)))./numleafs(sum(numIndividuals(1:i-1))+1:sum(numIndividuals(1:i))));
        varianceFit(i,:) = mean(fitness(sum(numIndividuals(1:i-1))+1:sum(numIndividuals(1:i))));
    end
    
    figure(1);
    plot(varianceAltura,numIndividuals,'color',color(k,:));
    hold on
    figure(2);
    plot(varianceArea,numIndividuals,'color',color(k,:));
    hold on
    figure(3);
    plot(varianceRamas,numIndividuals,'color',color(k,:));
    hold on
    figure(4);
    plot(varianceHojas,numIndividuals,'color',color(k,:));
    hold on
    figure(5);
    plot(varianceFit,numIndividuals,'color',color(k,:));
    hold on

end

figure(1)
legend('T1','T1','T1','T1','T2','T2','T2','T2','T3','T3','T3','T3')
xlabel('varianceHeight')
ylabel('sizePopulation')

figure(2)
legend('T1','T1','T1','T1','T2','T2','T2','T2','T3','T3','T3','T3')
xlabel('varianceArea')
ylabel('sizePopulation')

figure(3)
legend('T1','T1','T1','T1','T2','T2','T2','T2','T3','T3','T3','T3')
xlabel('varianceNumRamas')
ylabel('sizePopulation')

figure(4)
legend('T1','T1','T1','T1','T2','T2','T2','T2','T3','T3','T3','T3')
xlabel('variance(NumLeafsOnTop/NumLeafs)')
ylabel('sizePopulation')

figure(5)
legend('T1','T1','T1','T1','T2','T2','T2','T2','T3','T3','T3','T3')
xlabel('meanFitness')