function [varargout] = trackPopulation(poph, individuals, startGen, endGen, varargin)
%inputs:
%   -poph: population history, loaded with loadSimulation
%   -startGen: starting generation
%   -endGen: ending generation
%   -individuals: indexes of individuals of generation startGen. They
%                 will be tracked until generation endGen.
%                 Note that endGen can be lower or greater than startGen:
%                    -if greater, descendants will be tracked
%                    -if lower, ancestors will be tracked
%   -outputs: This is a list of char arguments, that can be in any order:
%        -'all': output indexes of all tracked individuals
%        -'generations': a cell array containing the indexes of the
%                        individuals, each cell corresponding to a
%                        generation 
%        -'amounts': the number of tracked indivudals in each generation
%        -'parents': return the list containing the parent's index for each
%                    node (0 for roots)
%        -'sons': return a cell array, each cell containing the sons'
%                 indexes (empty list for leaves)
%        -'plot': show a plot of the amount of tracked individuals per
%                 generation
%The list of output arguments will correspond to the list of outputs, so
%for example:
%  [all amounts parents generations] = trackPopulation(poph, 47, 10, 50, 'all', 'plot', 'amounts', 'parents', 'generations');  
%That is to say, the list of output arguments will match the list of chars
%modulo 'plot'

if numel(varargin)==0
  error('You have not specified any output!!!');
end

if floor(startGen)~=startGen
  error('startGen must be integer!!!');
end
if floor(endGen)~=endGen
  error('endGen must be integer!!!');
end

showplot       = false;
outamounts     = 0;
outparents     = 0;
outsons        = 0;
outgenerations = 0;
outall         = 0;
z=1;
for k=1:numel(varargin)
  switch varargin{k}
    case 'plot'
      showplot       = true;
    case 'parents'
      if outparents;     error('You have specified ''parents'' twice!!!'); end
      outparents     = z;
      z              = z+1;
    case 'sons'
      if outsons;        error('You have specified ''sons'' twice!!!'); end
      outsons        = z;
      z              = z+1;
    case 'amounts'
      if outamounts;     error('You have specified ''amounts'' twice!!!'); end
      outamounts     = z;
      z              = z+1;
    case 'generations'
      if outgenerations; error('You have specified ''generations'' twice!!!'); end
      outgenerations = z;
      z              = z+1;
    case 'all'
      if outall;         error('You have specified ''all'' twice!!!'); end
      outall         = z;
      z              = z+1;
  end
end

if outgenerations
  generations = cell(abs(endGen-startGen)+1,1);
end
if showplot || outamounts
  amounts     = zeros(abs(endGen-startGen)+1,1);
end

if startGen<endGen
  comp = @le;
  adv  = @plus;
  
else
  comp = @ge;
  adv  = @minus;
end

if isfield(poph.tree, 'parents')
  modeIndirect = false;
  parents = poph.tree.parents;
%   idxD    = poph.tree.idxD;
%   idxA    = poph.tree.idxA;
  gD      = poph.tree.gD;
  iD      = poph.tree.iD;
%   gA      = poph.tree.gA;
else
  modeIndirect = true;
  allowed = (poph.tree.gD>=min(startGen,endGen)) & (poph.tree.gD<=max(startGen,endGen));
  idxD    = poph.tree.idxD(allowed);
  idxA    = poph.tree.idxA(allowed);
  gD      = poph.tree.gD(allowed);
  iD      = poph.tree.iD(allowed);
%   gA      = poph.tree.gA(allowed);
  allowed = find(allowed);
%   try
%     sparseIndexes = sparse(1,idxD+1,1:numel(idxD),1, max(idxD)+1, numel(idxD));
%     parents       = full(sparseIndexes(idxA+1)); 
%     clear sparseIndexes;
%   catch ME %#ok<NASGU>
%     if strcmpi(ME.identifier, 'MATLAB:nomem')
      [nevermind, parents] = ismember(idxA, idxD); %#ok<ASGLU>
      parents = parents';
      clear nevermind;
%     else
%       rethrow(ME);
%     end
%   end
  clear idxD idxA;
end

%translate [generation, individuals] to indexes
indexes  = find(gD==startGen);
indexes  = indexes(ismember(iD(indexes), individuals));
%get 'sons' information
if (startGen<endGen) || outsons
  if isfield(poph.tree, 'sons')
    sons = poph.tree.sons;
  else
    sons = cell(size(parents));
    for z=1:numel(parents)
      pz = parents(z);
      if pz>0
        sons{pz}(end+1,1) = z;
      end
    end
  end
end

%do tracking
indexes         = reshape(unique(indexes),[],1);
allindexes      = indexes;
goon            = ~isempty(indexes);
while goon
  if (startGen>=endGen)
    %find parents
    news        = parents(indexes)';
    %wipe out parents outside bounds
    news        = news(news>0);
    news        = news((gD(news)>=endGen));
  else %(startGen<endGen)
    %fin sons
    news        = sons(indexes);
    news        = vertcat(news{:});
    %wipe out sons outside bounds
    news        = news(gD(news)<=endGen);
  end
  %end tracking when no more individuals can be tracked
  goon          = ~isempty(news);
  %add them!
  if goon
    allindexes  = unique(vertcat(allindexes, news));
  end
  indexes       = news;
end

%compile statistics by generation
if showplot || outamounts || outgenerations
  k= startGen;
  kn = 1;
  while comp(k, endGen)
    thisGen = gD(allindexes)==k;
    if outamounts || showplot
      amounts(kn) = sum(thisGen);
    end
    if outgenerations
      generations{kn} = allindexes(thisGen);
    end
    k=adv(k,1);
    kn=kn+1;
  end
end

%do output
if outparents
  varargout{outparents}            = parents;
  if modeIndirect
    varargout{outparents}          = allowed(varargout{outparents});
  end
end
if outsons
  varargout{outsons}               = sons;
  if modeIndirect
    varargout{outsons}             = allowed(varargout{outsons});
  end
end
if outgenerations
  varargout{outgenerations}        = generations;
  if modeIndirect
    varargout{outgenerations}      = cellfun(@(x)allowed(x), varargout{outgenerations}, 'uniformoutput', false);
  end
end
if outamounts
  varargout{outamounts}            = amounts;
end
if outall
  varargout{outall}                = allindexes;
  if modeIndirect
    varargout{outall}              = allowed(varargout{outall});
  end
end
if showplot
  figure;
  if startGen<=endGen
    gens = startGen:endGen;
    plot(gens,amounts);
  else
    gens = startGen:-1:endGen;
    plot(gens,amounts);
    set(gca, 'XDir', 'reverse');
  end
  xlabel('generations');
  ylabel('amount of individuals');
  if (startGen>=endGen)
    str =  'ascendants';
  else
    str = 'descendants';
  end
  title(sprintf('tracking from generation %d to generation %d (tracking individuals'' %s)', startGen, endGen, str));
  grid on;
end