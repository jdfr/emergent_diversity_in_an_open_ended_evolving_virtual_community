function genomeRandom=generateGenomeR(sgenome,meanGs,meanLongGenome)

for i = 1 : sgenome
    sizeGenome = randint(1,1,[meanLongGenome-10,meanLongGenome+10]);

    numG = randint(1,1,[meanGs-10,meanGs+10]);
    numBra = randint(1,1,[1 sizeGenome-numG]);
    if mod(numBra,2)~=0
        numBra = numBra + 1;
    end
    permutaciones = randperm(sizeGenome);
    posBra = sort(permutaciones(1:numBra));

    posG = randsrc(1,numG,setdiff(1:sizeGenome,posBra));
    posSignos = setdiff(1:sizeGenome,union(posBra,posG));

    genome(posBra(1:2:end-1)) = '[';
    genome(posBra(2:2:end)) = ']';
    genome(posG) = 'G';
    symbols = ['+','-'];
    genome(posSignos) = symbols(randint(1,length(posSignos),[1,length(symbols)]));


    k1 = find(genome == '+');
    k2 = find(genome == '-');
    if isempty(k1) && isempty(k2)
        genome(posBra(1)) = '+';
        genome(posBra(2)) = 'G';
    end
    genomeRandom{i} = genome;
end
