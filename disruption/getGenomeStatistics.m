function [stats allvals] = getGenomeStatistics(genome)

sgenome = numel(genome);

numGs = zeros(1, sgenome);
numBra = zeros(1, sgenome);
numSillyBra = zeros(1, sgenome);
numAnid = zeros(1, sgenome);
longGenome = zeros(1, sgenome);
for i = 1 : sgenome
    numGs(i) = sum(genome{i} == 'G');
    numBra(i) = sum(genome{i} == '[');
    numSillyBra(i) = sum(genome{i}(1:end-1) == '[' & genome{i}(2:end) == ']'); 
    numAnid(i) = max(cumsum((genome{i} == '[') + (-1 * (genome{i} == ']'))));
    longGenome(i) = length(genome{i});
end

allvals = struct('numGs', numGs, 'numBra', numBra, 'numSillyBra', numSillyBra, 'numAnid', numAnid, 'longGenome', longGenome);

meanGs = mean(numGs);
meanBra = mean(numBra);
meanSillyBra = mean(numSillyBra);
meanAnid = mean(numAnid);
meanLongGenome = mean(longGenome);

stats = struct('meanGs', meanGs, 'meanBra', meanBra, 'meanSillyBra', meanSillyBra, 'meanAnid', meanAnid, 'meanLongGenome', meanLongGenome, ...
  'round', struct('meanGs', round(meanGs), 'meanBra', round(meanBra), 'meanSillyBra', round(meanSillyBra), 'meanAnid', round(meanAnid), 'meanLongGenome', round(meanLongGenome)));
