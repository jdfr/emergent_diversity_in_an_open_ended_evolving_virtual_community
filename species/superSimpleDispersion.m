% function [daniel rene, md] = superSimpleDispersion(basedir, prefix)
function [daniel rene variance] = superSimpleDispersion(basedir, prefix)

file = [basedir filesep prefix 'matrixDist.mat'];
matrixDist = load(file);
matrixDist = matrixDist.matrixDist;

for k=1:numel(matrixDist)
  if isempty(matrixDist{k})
    matrixDist{k} = 0;
  end
end

[daniel rene variance] = cellfun(@measures, matrixDist);

function [daniel rene variance] = measures(matrixDist)

if not(size(matrixDist,1)==size(matrixDist,2))
  error('not square!!!');
end

if isempty(matrixDist)
  daniel = 0;
  rene = 0;
  return
end

%daniel measure: mean of the values of the matrix distance
daniel = mean(matrixDist(:));

variance = var(matrixDist(:));

sums = sum(matrixDist);
[mn mn] = min(sums);
%rene measure
rene = sums(mn)/numel(sums);

