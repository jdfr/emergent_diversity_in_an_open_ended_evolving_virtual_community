function genomeRandom=generateGenomeRTipo2(numgenomes,stats,meanGs,meanBra,meanAnid,meanLongGenome)
sgenome = numgenomes;
if not(isempty(stats))
  meanGs = stats.round.meanGs;
  meanLongGenome = stats.round.meanLongGenome;
  meanBra = stats.round.meanBra;
  meanAnid = stats.round.meanAnid;
end
genomeRandom = cell(sgenome, 1);

meanUGs = meanGs/meanLongGenome;
meanUBra = meanBra/meanLongGenome;

for i = 1 : sgenome
    sizeGenome = randint(1,1,[meanLongGenome-10,meanLongGenome+10]);
    genome = repmat('x', 1, sizeGenome);
    
    numG = round(meanUGs * sizeGenome);
    %numG = randint(1,1,[meanGs-10, min(meanGs+10, sizeGenome-2)]); % -2 in order to have at least a pair of []
    permutaciones = randperm(sizeGenome);
    posG = permutaciones(1:numG);
    
    numBra = round(2 * meanUBra * sizeGenome);
    %numBra = randint(1,1,[2 sizeGenome-numG]);
    posBra = permutaciones(numG+1:min(sizeGenome,numG+numBra));
    if(mod(length(posBra),2)~=0)
        posBra(end) = [];
        numBra = numBra - 1;
    end
    posBra = sort(posBra);
    %posBra = reshape(sort(reshape(posBra, 2, [])), 1, []);
    genome(posBra(numBra/2)) = '[';
    genome(posBra(numBra/2 + 1)) = ']';
    rp = randperm(numBra/2 - 1);
    leftAnidIndex = rp(1:meanAnid-1);
    rp = randperm(numBra/2 - 1) + numBra/2 + 1;
    rightAnidIndex = rp(1:meanAnid-1);
	genome(posBra(leftAnidIndex)) = '[';	
    genome(posBra(rightAnidIndex)) = ']';
    posBra([numBra/2, numBra/2 + 1, leftAnidIndex, rightAnidIndex]) = [];
    for(j=1:2:length(posBra))
        genome(posBra(j)) = '[';
        genome(posBra(j+1)) = ']';
    end
    
    posSignos = permutaciones(numG+numBra+1:sizeGenome);

    genome(posG) = 'G';
    %genome(posBra(1:2:end-1)) = '['; % genoma Tipo 3
    %genome(posBra(2:2:end)) = ']';
    symbols = ['+','-'];
    genome(posSignos) = symbols(randint(1,length(posSignos),[1,length(symbols)]));


    noPlus = isempty(find(genome == '+', 1));
    noMinus = isempty(find(genome == '-', 1));
    if noPlus && noMinus
        genome(posBra(1)) = '+';
        genome(posBra(2)) = 'G';
    end
    genomeRandom{i} = genome;
    %fprintf('%s\n', genome);
    if(~corchetesBalanceados(genome))
        error('regla con corchetes no balanceados: %s', genome);
    end
    if(isempty(find(genome=='[', 1)))
        error('genoma sin corchetes: %s', genome);
    end
end
