function axiom = sanitizeString(axiom, randomExpression)
%remove all redundancies from a string, by means of regular expressions (it
%cannot find all redundancies)
persistent regularExpressions1 regularExpressions2;
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
  end
  if (~isempty(axiom)) % && (numel(axiom)~=antNum)
    axiom  = toAxiom(simplifyTree(toTree(axiom),false));
  end
 goon   = numel(axiom)~=antNum;
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

function tree = simplifyTree(tree, notRoot)
%simplify an axiom in tree form
if iscell(tree)
  superCells = find(areSuperCells(tree));
  %integrate supercells in lower level, from last one to first one, to
  %preserve index integrity in the process
  for k=numel(superCells):-1:1
    idx = superCells(k);
    tree = [tree(1:idx-1) tree{idx} tree(idx+1:end)];
  end
end
if notRoot %do not alter the last one if this is the root node, since the results are critically different when the rule is applied 3 times
  %this is already done by regular expressions, but they cannot de-nest
  %complex expressions
  while iscell(tree) && (numel(tree)==1) && iscell(tree{1})
    tree = tree{1};
  end
  if iscell(tree) && (numel(tree)>1) 
    if all(cellfun(@iscell, tree(1:end-1))) && (~iscell(tree{end}))
      %wrap the last node, creating a superCell to be wiped in the next
      %pass
      tree{end} = {tree{end}};
    elseif iscell(tree{end})
      %de-nest the last node if necessary
      tree = [tree(1:end-1) tree{end}];
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
      tree{k} = simplifyTree(tree{k}, true);
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