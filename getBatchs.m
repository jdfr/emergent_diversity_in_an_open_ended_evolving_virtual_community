function out = getBatchs(bd, batchs, alphs)

if not(exist('bd', 'var')) || (isempty(bd) && not(ischar(bd)))
  bd = 'lsystemdani\dataset\resumenok\test';
end

if not(exist('batchs', 'var')) || isempty(batchs)
  batchs = {'LONG', 'LONGG', 'LONGOLD', 'LONGOLD2'};
end

if ischar(batchs)
  batchs = {batchs};
end

batchs   = cellfunc(@(x)[bd x], batchs);
%batchs   = cellfunc(@(x)[bd x], {'110NB'});

allreses = cellfunc(@(x)load([x filesep 'allres.mat']), batchs);
allreses = cellfunc(@(x)x.allres, allreses);
nisd = not(cellfun(@(x)isfield(x, 'disps'), allreses));
if any(nisd) && not(all(nisd))
  nisd = find(nisd);
  for k=1:numel(nisd)
    allreses{nisd(k)}(1).disps = [];
  end
end
allreses = vertcat(allreses{:});

alphas = [allreses.alphaG]';
if not(exist('alphs', 'var')) || isempty(alphs)
  alphs = unique(alphas);
end
out = allreses(ismember(alphas, alphs));
%all075 = allreses;
