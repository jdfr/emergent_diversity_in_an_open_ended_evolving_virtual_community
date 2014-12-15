function [P colors params]=fitnessSimilarity(P,params,varargin)

sR = numel(P.genome);

P.randN(1:end) = 0;

colors = rand(sR, 3)*0.8;

if isempty(params.treeToMatch)
  error('The parameter treeToMatch cannot be empty!!!');
elseif ischar(params.treeToMatch)
  [raster, offset, counts, dimension, isLeaf, maxiy]=ls2(params.treeToMatch, params.T);
  params.treeToMatch    = {raster, offset, counts, dimension, isLeaf, maxiy, params.treeToMatch};
  clear raster offset counts dimension isLeaf maxiy;
end

%prepare pixels of reference tree
if iscell(params.treeToMatch{1})
  pixelsReferencia      = [params.treeToMatch{1}{1:2}];
  pixelsReferencia(:,1) = pixelsReferencia(:,1)-params.treeToMatch{2}(1);
  pixelsReferencia(:,2) = pixelsReferencia(:,2)-params.treeToMatch{2}(2);
  params.treeToMatch    = [{pixelsReferencia} params.treeToMatch];
else
  pixelsReferencia      = params.treeToMatch{1};
end

for k=1:sR
  %prepare pixels of current tree
  pixelsTree      = [P.raster{k}{1:2}];
  pixelsTree(:,1) = pixelsTree(:,1)-P.offset(k,1);
  pixelsTree(:,2) = pixelsTree(:,2)-P.offset(k,2);
  %fitness is similitud
  P.fitness(k)    = similitudFast(pixelsReferencia, pixelsTree);
end
  
