function meanDiruptionAllGenerations(basedir)

poph = loadSimulation(basedir);

numGenerations = max(poph.generation);

for i = 1 : numGenerations
    numIndividuals(i) = sum(poph.generation == i);
end

disruption = poph.tree.disruption;

for i = 1 : numGenerations-1
    dis = setdiff(disruption(sum(numIndividuals(1:i))+1:sum(numIndividuals(1:i+1))),0);
    meanDisrupt(i) = mean(dis);
%     maxDisrupt(i) = max(disruption(sum(numIndividuals(1:i))+1:sum(numIndividuals(1:i+1))));
%     minDisrupt(i) = min(disruption(sum(numIndividuals(1:i))+1:sum(numIndividuals(1:i+1))));
end

plot(meanDisrupt,'b');
% hold on;
% plot(maxDisrupt,'r');
% hold on;
% plot(minDisrupt,'g');

% legend('mean disruption', 'masimum disruption', 'minimum disruption')
xlabel('generations')
ylabel('mean disruption')

