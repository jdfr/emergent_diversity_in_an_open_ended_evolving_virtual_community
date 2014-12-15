function axiom = sanitizeStringNEW(axiom, randomExpression)
%remove all redundancies from a string, by means of regular expressions (it
%cannot find all redundancies)
persistent regularExpressions1;
if isempty(regularExpressions1)
  regularExpressions1 = {...
    ... %to remove trailing +/signs
    {false, '([+-]+)\]', ']'}, ... 
    ... %to remove unnecessary brackets (brackets that are neededlessly nested)
    {false, '\[(?<open>\[+)(?<inner>[^\[\]]*)(??@repmat(''\]'', size($<open>)))\]', '[$<inner>]'}, ... 
    ... %to scorch out void brackets (spared by the previous expression)
    {false,  '\[\]', ''}...
    ... %again, to remove trailing +/signs, that migh become exposed after void brackets are deleted
    {false, '([+-]+)\]', ']'}, ... 
    {true, '(\+\-)|(\-\+)', ''}, ...
    ... %to remove simple repeated brackets [CAD][CAD] and [CAD]CAD that do not contribute anything
        ... %to remove simple expressions of the form [CAD]CAD that are redundant [G][G-]
        ...%{true, '\[([^\[\]]+)\]\[\1\]', '[$1]'}, {true, '\[([^\[\]]+)\]\1', '$1'}, ...
    {true, '\[([^\[\]]+)\]((\[([^\[\]]+)\])*)((\[\1\])|\1)', '$2$3'}, ...
    };
end
if nargin<2
  randomExpression = false;
end

if randomExpression
  error('Currently, we do not support sanitization of strings to be randomly expressed!!!!');
end

%remove all redundancies from a string
goon       = true;
while goon
  antNum   = numel(axiom);
  goon_inner = true;
  while goon_inner
    antNumInner = numel(axiom);
    axiom    = simplifyAxiom(axiom,regularExpressions1);
    goon_inner   = numel(axiom)~=antNumInner;
%     fprintf('AFTER REGEXP:   %s\n', axiom);
  end
  chg = false;
  if (~isempty(axiom)) % && (numel(axiom)~=antNum)
    t = toTree(axiom);
    [t chg] =  simplifyTree(t, false);
    axiom  = toAxiom(t);
    t=[]; %#ok<NASGU>
%     fprintf('AFTER SIMPLIFY: %s\n', axiom);
  end
 goon   = chg || (numel(axiom)~=antNum);
end
end

% function changeDir = minimizeChangeDir(changeDir)
% totalChange = sum(changeDir=='+')*2-numel(changeDir);
% changeDir   = [repmat('+', 1, max(0, totalChange)) repmat('-', 1, max(0, -totalChange))];
% end
% 
function axiom = simplifyAxiom(axiom, regularExpressions)

for m=1:numel(regularExpressions)
  if regularExpressions{m}{1}
    goon      = true;
    while goon
      lastNum = numel(axiom);
      axiom   = regexprep(axiom, regularExpressions{m}{2:3});
      goon    = (numel(axiom)~=lastNum);
    end
  else
    axiom     = regexprep(axiom, regularExpressions{m}{2:3});
  end
end
end

function root     = toTree(axiom)
tokens            = regexp(axiom, '(\[)|(|])|([^\[\]]+)', 'match');
[starts, ends]    = cellfun(@(x)deal(x(1)=='[', x(1)==']'), tokens);
others            = ~(starts | ends);
nestingLevel      = cumsum(starts-ends);

root              = populateSubTree(0,     1, numel(tokens));
  function tree   = populateSubTree(level, s, e)
  if s==e
    tree          = tokens(s);
  else
    thisLevel     = find( (starts(s:e)&(nestingLevel(s:e)==(level+1))) | (others(s:e)&(nestingLevel(s:e)==level)) );

    tree          = cell(1,numel(thisLevel));
    for k=1:numel(tree)
      current     = s+thisLevel(k)-1;
      if others(current)
        tree{k}   = tokens{current};
      else
        next      = current+1;
        endToken  = next+find(ends(next:e) & (nestingLevel(next:e)==level), 1)-2;
        tree{k}   = populateSubTree(level+1, next, endToken);
      end
    end
  end
  end
end


% function [tree n] = toTree(axiom, n)
% %convert the axiom from string to tree form
% tree = {};
% while n<=numel(axiom)
%   switch axiom(n)
%     case {'+', '-', 'G', 'g'}
%       tree{end+1} = axiom(n); %#ok<AGROW>
%     case '['
%       [subtree n] = toTree(axiom, n+1);
%       tree{end+1} = subtree; %#ok<AGROW>
%     case ']'
%       return
%   end
%   n=n+1;
% end
% end

function [tree changed] = simplifyTree(tree, notRoot)
changed = false;
%simplify an axiom in tree form
if iscell(tree)
  superCells = find(areSuperCells(tree));
  %integrate supercells in lower level, from last one to first one, to
  %preserve index integrity in the process
  for k=numel(superCells):-1:1
    idx = superCells(k);
    tree = [tree(1:idx-1) tree{idx} tree(idx+1:end)];
    changed = true;
  end
end

[tree chg] = changeBranches(tree);
changed = changed || chg;

if notRoot %do not alter the last one if this is the root node, since the results are critically different when the rule is applied 3 times
  %this is already done by regular expressions, but they cannot de-nest
  %complex expressions
  while iscell(tree) && (numel(tree)==1) && iscell(tree{1})
    tree = tree{1};
  end
  
  if iscell(tree) && (numel(tree)>1) 
    arecells = cellfun(@iscell, tree);
    
    if arecells(1)
      %create a supercell to be wiped in the next pass, using the first
      %cells and wrapping the rest. This subsumes the previous code
      notcell = find(not(arecells),1);
      if isempty(notcell)
        error('This should not happen!!!!');
      end
        tree = [tree(1:(notcell-1)) {tree(notcell:end)}];
        changed = true;
    elseif arecells(end)
      %de-nest the last node if necessary
      tree = [tree(1:end-1) tree{end}];
      changed = true;
    end
  end
end
% if notRoot %do not de-nest the last one if this is the root node, since
% the results are critically different when the rule is applied 3 times
%   %de-nest the last node if necessary
%   if iscell(tree) && (numel(tree)>1) && iscell(tree{end})
%     tree = [tree(1:end-1) tree{end}];
%   end
% end
if iscell(tree)
  %recursively apply simplifyTree
  for k=1:numel(tree)
    if iscell(tree{k})
      [tree{k} chg] = simplifyTree(tree{k}, true);
      changed = changed || chg;
    end
  end
end
end

function scs = areSuperCells(tree)
%supercells are nodes whose all sons are not leafs. They can safely be
%brought up one level: for example, [[G][+G]] can safely become [G][+G]
if iscell(tree)
  scs = cellfun(@iscell, tree);
  for k=1:numel(tree)
    if scs(k)
      subtree = tree{k};
      for m=1:numel(subtree)
        scs(k) = iscell(subtree{m});
        if ~scs(k)
          break;
        end
      end
    end
  end
else
  scs = false;
end
end

function axiom = toAxiom(tree)
%convert an axiom from tree form to string
axiom = cell(size(tree));
for k=1:numel(tree)
  if iscell(tree{k})
    axiom{k} = ['[' toAxiom(tree{k}) ']'];
  else
    axiom{k}  = tree{k};
  end
end
axiom = [axiom{:}];
end

function [tree changed]= changeBranches(tree)
changed = false;
arecells = cellfun(@iscell, tree);
indscells = find(arecells);
ni = numel(indscells);
if ni>1
  df = diff(indscells);
  consec = df==1;
  init=indscells(1);
  nc = ni-1;
  for k=1:nc
    if consec(k) && (k<nc)
      continue;
    end
    if consec(k)
      fin = indscells(end);
    else
      fin = indscells(k);
    end
    if fin>init
      [tree(init:fin) chg] = removeCommon(tree(init:fin), tree, fin);
      changed = changed || chg;
    end
    init = indscells(k+1);
  end
end
end

function [tree changed] = removeCommon(tree, alltree, fintree)
changed = false;

%this function operates over a subtree that is a set of consecutive cells

arenull = cellfun('prodofsize', tree)==0;
tn = tree(arenull);

tree = tree(not(arenull));

firsts = cellfun(@(x)x{1}, tree, 'uniformoutput', false);

%only operate on those having initial characters, not subcells
arechar = cellfun(@ischar, firsts);

chars = firsts(arechar);

if numel(chars)<2
  return
end

subno = tree(not(arechar));

sub = tree(arechar);

lens =cellfun('prodofsize', chars);

[lens idx] = sort(lens, 'descend');

sub = sub(idx);
chars = chars(idx);

%this is not optimal: we should not look for perfect matches, but for the
%longest common prefix. However, this might turn out to have non-optimal
%implications. Let's settle on look for whole prefixes
for k=1:numel(idx)-1
  if isempty(sub{k}) || isempty(sub{k}{1}) || not(ischar(sub{1}{1}))
    continue;
  end
  for z=k+1:numel(idx)
    if (numel(sub{z})>0) && not(isempty(sub{z}{1})) && ischar(sub{z}{1}) && not(isempty(strmatch(sub{z}{1}, sub{k}{1}))) %strmatch(chars{z}, chars{k})
      if lens(z)==lens(k)
        %if the strings are equal, we can do this very easily
        sub{k} = [sub{k}(1) {sub{z}(2:end)} sub{k}(2:end)];
        sub{z}={};
        changed = true;
      else
        %this case is more difficult to process in general. We only process
        %it if it is easy
        if numel(sub{z})==1 %do this only if the the string "z" is the only content of sub{z}
          %in this case, we can safely remove this, since its only contents
          %are the substring
          sub{z}={};
          changed = true;
        end
      end
    end
  end
end

if not(changed) && (fintree==(numel(alltree)-1)) && ischar(alltree{end})
  for k=1:numel(sub)
    if (numel(sub{k})==1) && not(isempty(strmatch(sub{k}{1}, alltree{end})))
      sub{k}={};
      changed=true;
    end
  end
end

tree = [tn subno sub];

end
