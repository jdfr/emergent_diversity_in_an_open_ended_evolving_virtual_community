%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G indexNewG terclus P parentdata] = evolucion(params,P,generacion,rangeid,generationInRange,terclus,mejorEnConserva)
%
%  range: 2-tuple marking the range which will be replaced. The new
%  individuals do not immediately replace the population in that range, but
%  are stored in G.
%  The indexes of truly new individuals are in indexNewG (the others are
%  just straightforward copies of individuals in P, i.e., they haven't
%  undergone any mutation at all). IMPORTANT: indexNewG is relative to G,
%  not P!!!!
%  Their values 'b' and 'ra' for new individuals are stored in 'bG' and
%  'raG' respectively: the new individuals have initially 0 (as they must
%  be evaluated), but old individuals' values are straightforwardly copied.
%
%  All this stuff is a nuissance, but it enables the code to overlap
%  execution for several subpopulations, thus being able to use the cluster
%  for a subpopulation while applying genetic operators to another
%  subpopulation
%
%   prob : probabilidades de cada operador
%
%      1   elitismo             =  
%      2   op_alteracion        A simbolo  posicion  
%      3   op_dupaleatoria      R segmento posicion
%      4   op_dupnivel          L segmento posicion
%      5   op_dupsecuencia      C segmento posicion
%      6   op_eliminacion       D posicion
%      7   op_insercion         I simbolo  posicion
%
% el patrón utilizado en el fichero arbol.txt es:
% generacion progenitor operador_segmento_posicion
% el patrón operador_segmento_posicion podrá aparecer repetido tantas veces
% como operadores se hayan aplicado sobre el individuo.
% Ejemplo un operador: 1 186 N[FGF-GGF]3
% Ejemplo varios operadores: 1 16 C+8N[+GFF++F]3
% Cuando un operador se aplica al individuo pero no modifica la cadena que
% lo representa se representará por el símbolo '=' y el parámetro posición no se almacena. 
% Ejemplo: 1 95 S=I+8 (generacion progenitor operador_=_operador_segmento_posicion)
% Por último, el patrón para el operador elitista es:
% generacion progenitor =
% Ejemplo: 1 207 =

maxLongsImpl = params.maxLongs;
if ~params.maxLongs.penalizeOperators
  maxLongsImpl.maxInsert  = [];
  maxLongsImpl.maxGs      = [];
  maxLongsImpl.maxString  = [];
end

variablePopSize           = params.variablePopSize;
elitista                  = (rand<=params.prob(1));       % estrategia elitista (mantener mejor)
elitista                  = elitista && (terclus.ranges(rangeid,1)==1); %only if the range starts in the 1st place

if variablePopSize
  if size(terclus.ranges,1)~=1
    error('If using variable population size, there must be only a range!!!');
  end
  alphaF = params.alphaF;
  numSonsMeasure      = P.fitness.^alphaF;
  areImaginary        = imag(numSonsMeasure)~=0;
  numSonsMeasure(areImaginary) = 0;
  %determine the individuals that *won't* leave descend
  toPenalize          = penalizeDescendency(P, params.maxLongs, params);
  switch params.variableNumSonsMode %= %'mixed'; 'round'; 'ceil'; 'floor'
    case 'mixed'
      toReward            = ~toPenalize;
      numToReward         = sum(toReward);
      %minimal guaranteed number of sons
      numSons             = floor(numSonsMeasure);
      %probability of having an additional son
      probAdicional       = numSonsMeasure-numSons;
      %cut low probability values
      probAdicional((probAdicional<params.minimalDescendProbability) & (numSons==0)) = params.minimalDescendProbability;
      %number of sons after lottery
      numSons             = numSons + (rand(size(numSonsMeasure))<=probAdicional);
      %remember to penalize some individuals
      numSons(toPenalize) = 0;
      %if nobody is able to leave descend, stop it right now
      if (~elitista) && ( ((params.minimalDescendProbability==0) && (~any(numSons(toReward)))) || (all(toPenalize)) )
        G = [];
        return;
      end
      n = 1;
      %if we are able to have sons but do not have them yet, try it
      while ~any(numSons)
        numSons(toReward) = numSons + (rand(numToReward,1)<=probAdicional(toReward));
        n                 = n+1;
        if n>1e5
          error('Hey you have set a too low probability for additional sons, I am going to play this game no more!!!');
        end
      end
      clear toReward numToReward probAdicional n;
    case {'round', 'ceil', 'floor'}
      numSonsMeasure(toPenalize)   = 0;
      if isempty(params.variableNumSonsInterval)
        numSons                           = feval(params.variableNumSonsMode, numSonsMeasure);
        numSons(toPenalize)               = 0;
      else
        randshift                         = diff(params.variableNumSonsInterval)*rand(size(numSonsMeasure)) + ...
                                                 params.variableNumSonsInterval(1);
        numSonsMeasure                    = numSonsMeasure+randshift;
        numSonsMeasure(numSonsMeasure<0)  = 0;
        if params.preserveUniqueTree && (numel(numSonsMeasure)==1) && (feval(params.variableNumSonsMode, numSonsMeasure)<1)
          %this is to not abort simulations too easy
          numSons                         = 1;
        else
          numSons                         = feval(params.variableNumSonsMode, numSonsMeasure);
          numSons(toPenalize)             = 0;
        end
      end
      clear randshift;
    otherwise
      error('variableNumsonsNode value <%s> not understood!!!', any2str(params.variableNumSonsMode));
  end
  %prepare indexes of parents to be copied to the next population
  idxInP                  = zeros(sum(numSons),1);
  endR                    = cumsum(numSons);
  startR                  = endR-numSons+1;
  for k=1:numel(P.fitness)
    idxInP(startR(k):endR(k)) = k;
  end
  clear numSonsMeasure areImaginary numSons endR startR toPenalize;
  if elitista
    idxInP                = [1; idxInP];
  end
  %remake the settings for terclus
  terclus.ranges          = [1 numel(idxInP)];
  terclus.rangeSizes      =    numel(idxInP);
  range                   = terclus.ranges(rangeid,:);
  rangesExpanded          = ones(size(P.genome));
  rangoIsCalculated       = false;
else
  
  range                   = terclus.ranges(rangeid,:);

  rangesExpanded          = zeros(size(P.genome));
  for k=1:size(terclus.ranges,1)
    rangesExpanded(terclus.ranges(k,1):terclus.ranges(k,2)) = k;
  end

  %classical way: stocastically select individuals according to their
  %fitnes
  rango                   = calculateRangoSeleccion(P,P.fitness,params.presion,params.maxLongs,params);
  rangoIsCalculated       = true;
  numG                    = range(2)-range(1)+1;
  idxInP                  = seleccion(rango,numG);
end

G           = getSubStruct(P, idxInP, 'extensive');
parentdata.width = G.dimensions(:,2);
%if params.save.disp.any
  parentdata.idxInP = idxInP;
%end
% G           = repmatStruct(params.indivStruct, size(P.genome));
% idxInP      = zeros(size(G.genome));
maskNew     = false(size(G.genome));
mutated     = false(size(G.genome));
if (elitista)
  G.genome{1} = mejorEnConserva.data.genome{1};
  parentdata.width(1) = mejorEnConserva.data.dimensions(1,2);
  maskNew(1)  = true; %false;
  mutated(1)  = false; %false;
  idxInP(1)   = mejorEnConserva.num;
  %fprintf(params.farb,'%d %d %d = %d %d %d    ANTEPASADO = %s\n',generacion,rangeid,range(1),mejorEnConserva.gen, mejorEnConserva.rind,mejorEnConserva.num,mejorEnConserva.data.genome{1});
  %fprintf(params.farb,'%d %d %d = %d %d %d 0 0\n',generacion,rangeid,1,mejorEnConserva.gen, mejorEnConserva.rind,mejorEnConserva.num);
  fprintf(params.farb,'%d %d %d = %d %d %d\n',generacion,rangeid,1,mejorEnConserva.gen, mejorEnConserva.rind,mejorEnConserva.num);
end
Pgenome = P.genome;

prob_ind = params.prob;
probIndexes   = 2:(numel(params.prob));
justOne = params.mutationOrder(1)=='j';
switch params.mutationOrder(1)
  case {'n', 'j'} %'normal' 'justone'
    permutation = 1:numel(probIndexes);
  case 'r' %'random'
    permutation = randperm(numel(probIndexes));
  otherwise
    error('Parameter mutationOrder <%s> not understood!!!', any2str(params.mutationOrder));
end
randomize = (isfield(params.T, 'randomize') && params.T.randomize);
for j=elitista+1:length(G.genome)      % crear individuos a partir de la anterior
    i    = idxInP(j);
%     i    = seleccion(rango,1);
%     G.genome{j} = P.genome{i};
%     idxInP(j) = i;
    prob = mutar(prob_ind); %#ok<NASGU> %ESTO HACE FALTA PARA MANTENER REPLICABILIDAD
    
    cad = '';
    
    if justOne
      %if only a genetic operator is to be applied, select just one
      prob_ind = rand(size(params.prob))<=params.prob;
      prob_ind(1) = 0;
      if sum(prob_ind)>1
        indexes_prob_ind = find(prob_ind);
        indexes_prob_ind = indexes_prob_ind(randperm(numel(indexes_prob_ind)));
        prob_ind         = zeros(size(prob_ind));
        prob_ind(indexes_prob_ind(1)) = 1;
      end
    end
    %apply genetic operators. All of them are fairly similar, so a loop is
    %used to keep code tidier
    for zz=1:(numel(params.prob)-1)
      z = permutation(zz);
      if (rand<=prob_ind(probIndexes(z)))
        switch probIndexes(z)
          case 8 %only 8-th operator requires this parameter to be set
            if ~rangoIsCalculated
              rango             = calculateRangoSeleccion(P,P.fitness,params.presion,params.maxLongs,params);
              rangoIsCalculated = true;
            end
            cad2                = Pgenome{seleccion(rango,1)};
          otherwise
            cad2 = '';
        end
        %apply the operator coded by probIndexes(z)
        [G.genome{j},seg,n,operatorChar] = cod_implicita(probIndexes(z),G.genome{j},cad2, maxLongsImpl, params);
%         [G.genome{j},seg,n,operatorChar,disruption] = cod_implicita(probIndexes(z),G.genome{j},cad2, maxLongsImpl, params);
        %record tree
        if (seg~='=')
            cad = [cad,operatorChar,seg,num2str(n)]; %#ok<AGROW>
            maskNew(j) = randomize || not(strcmp(G.mingnm{j}, sanitizeString(G.genome{j})));
            mutated(j) = true;
        else
            cad = [cad,operatorChar,seg]; %#ok<AGROW>
        end
      end
    end
    %pause
    indexNewInd = j+range(1)-1;

    ranOld = rangesExpanded(i);
    genOld = generationInRange(ranOld);
    if ~isempty(cad)
        %añadir al fichero de texto: generación, indice del nuevo
        %individuo, indice del viejo individuo, cambios
        fprintf(params.farb,'%d %d %d %s %d %d %d\n',generacion,rangeid,indexNewInd,cad,genOld,ranOld,i);
%         fprintf(params.farb,'%d %d %d %s %d %d %d %d %f\n',generacion,rangeid,indexNewInd,cad,genOld,ranOld,i,n,disruption);
    else
        %añadir al fichero de texto
        fprintf(params.farb,'%d %d %d = %d %d %d\n',generacion,rangeid,indexNewInd,genOld,ranOld,i);
%         fprintf(params.farb,'%d %d %d = %d %d %d 0 0\n',generacion,rangeid,indexNewInd,genOld,ranOld,i);
    end
end

%if variablePopSize && not(params.save.disp.any)
%    P = [];
%end

%get the indexes for new individuals
indexNewG       = find(maskNew);
%if params.save.disp.any
  parentdata.maskNew = maskNew;
  parentdata.mutated = mutated;
%end

%put default values in all new individuals
fieldsNoGenome  = fieldnames(G);
fieldsNoGenome  = fieldsNoGenome([2, 4:end]); %do not assign the genome field nor the randN (xpos) field
% G = assignSubStruct(G, ~maskNew, 'extensive', P, idxInP(~maskNew), 'extensive', fieldsNoGenome);
G               = assignSubStruct(G, maskNew, 'extensive', params.indivStruct, 1, 'extensive', fieldsNoGenome);

if elitista
  G      = assignSubStruct(G, 1, 'extensive', mejorEnConserva.data, 1, 'extensive');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rango = calculateRangoSeleccion(P,f,presion,maxLongs,params)
sf=length(f);
rango = zeros(sf,2);
M = max(f);
m = min(f);

if (m~=M)
    rango(1:sf,2)       = 1+presion*(((f(1:sf)'-m)/(M-m))-1);

    toPenalize          = penalizeDescendency(P, maxLongs, params);
    
    rango(toPenalize,2) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function toPenalize = penalizeDescendency(P, maxLongs, params)
toPenalize = false(size(P.genome));
if maxLongs.penalizeFitness
  if ~isempty(maxLongs.maxGs)
    numGs         = cellfun(@(x)sum(x=='G'), P.genome);
    toPenalize    = toPenalize | (numGs>maxLongs.maxGs);
  end
  if ~isempty(maxLongs.maxString)
    longs         = cellfun(@numel, P.genome);
    toPenalize    = toPenalize | (longs>maxLongs.maxString);
  end
  if ~isempty(params.heightT)
    heights       = cellfun(@safeMax, P.maxiy);
    toPenalize    = toPenalize | (heights>params.heightT);
  end
  if (params.negYThreshold<inf)
    negYCounts    = P.counts(:,3);
    toPenalize    = toPenalize | (negYCounts>params.negYThreshold);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = safeMax(array)
if isempty(array);  m=int32(0);  else  m=max(array); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inds = seleccion(rango, ntimes)
%this function selects an amount "ntimes" of individuals with a probability
%directly proportional to its fitness
%
% rango: matrix num_individualsX2, each row of the for [0..num], expressing
% the probability of selecting each individual.

if ntimes==1
  r = randnumber(rango);%rango(:,1) + ( rand(size(rango,1),1).*(rango(:,2) - rango(:,1)) );
else
  %generates a matrix of size "size(rango,1)Xntimes" such that the i-th row
  %is an array of "ntimes" random numbers in the range [rango(i,1)...rango(i,1)]  
  r = bsxfun(@plus, rango(:,1), bsxfun(@times, rand(size(rango,1),ntimes), (rango(:,2) - rango(:,1))) );
end

[nevermind,inds] = max(r); %#ok<ASGLU>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OG = mutar(p)
% funcion para decidir que operadores se aplican
%       p contiene las probabilidades absolutas de cada operador: pDS, pDN,
%       pDA, pI, pE, pA
%       OG es un vector binario que especifica con 1 que operador debe
%       aplicarse

OG = rand(size(p))<=p;