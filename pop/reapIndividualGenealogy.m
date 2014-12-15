function g = reapIndividualGenealogy(basedir, generation, rangeid, individual)
%reap the genealogy of a individual specified by its generation, range and
%index in the population. basedir is the directory of the simulation

g = reapGenealogyData(basedir);

indexes   = zeros(0,1);
continuar = true;
ind       = '';
while continuar
  idx = find((g.gDescendant==generation) & (g.iDescendant==individual));%g.idxDescendant==g.idxDF(generation, rangeid, individual), 1);
  continuar = ~isempty(idx);
  if continuar
    generation = g.gAncestor(idx);
    %rangeid    = g.rAncestor(idx);
    individual = g.iAncestor(idx);
    if ~strcmp(ind, g.ancestor(idx))
      indexes  = [indexes; idx]; %#ok<AGROW>
      ind      = g.ancestor(idx);
    end
  end
end

%g  = rmfield(g, {'idxDF'});%, 'idxDescendant'});
fs = fieldnames(g);
for k=1:numel(fs)
  g.(fs{k}) = g.(fs{k})(indexes);
end

function g = reapGenealogyData(basedir)
%load in structure 'g' the contents of arbol.txt 

%given a list of values, this function provides the minimum power of ten 
%that is greater than the maximum value
nextDigit = @(values) power(10, (1+ceil(log10(max(values)))));

farb = fopen([basedir filesep 'arbol.txt'], 'r');
datos = textscan(farb, '%d %d %d %s %d %d %d %*s %*s %s');
fclose(farb);
g = struct('gDescendant', [], ...
           'rDescendant', [], ...
           'iDescendant', [], ...
           ...%'idxDescendant', [], ...
           ...%'idxDF', [], ...
           'change', [], ...
           'gAncestor', [], ...
           'rAncestor', [], ...
           'iAncestor', [], ...
           ...%'idxAncestor', [], ...
           ...%'idxAF', [], ...
           'ancestor', []...
           );
g.gDescendant   = double(datos{1}); datos{1} = [];
g.rDescendant   = double(datos{2}); datos{2} = [];
g.iDescendant   = double(datos{3}); datos{3} = [];
%unique value identifying descendant
% powerRD         = nextDigit(g.iDescendant);
% g.idxDescendant = g.iDescendant   + powerRD * g.rDescendant;
% powerGD         = nextDigit(g.idxDescendant);
% g.idxDescendant = g.idxDescendant + powerGD * g.gDescendant;
g.change        = datos{4};         datos{4} = [];
g.gAncestor     = double(datos{5}); datos{5} = [];
g.rAncestor     = double(datos{6}); datos{6} = [];
g.iAncestor     = double(datos{7}); datos{7} = [];
% %unique value identifying ancestor
% powerRA         = nextDigit(g.iAncestor);
% g.idxAncestor   = g.iAncestor   + powerRA * g.rAncestor;
% powerGA         = nextDigit(g.idxAncestor);
% g.idxAncestor   = g.idxAncestor + powerGA * g.gAncestor;
g.ancestor      = datos{8};         datos{8} = []; %#ok<NASGU>

clear basedir farb datos nextDigit;
% g.idxDF         = @(gen, rind, idx) idx+powerRD*rind+powerGD*gen;
% g.idxAF         = @(gen, rind, idx) idx+powerRA*rind+powerGA*gen;


