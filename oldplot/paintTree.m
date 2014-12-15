% <<<<<<< .mine
function array = paintTree(axiom, mode, varargin)
% =======
% function paintTree(axiom, mode, alsoLeafs)
% >>>>>>> .r201

%paintTree(axiom,struct('shadowing','transparentBranches','angle',22))

if ischar(axiom)
  if ~exist('mode', 'var')
    mode = 'canonical';
  end
% <<<<<<< .mine
%  [raster, offset, counts, dimension, isLeaf, maxiy, leaves]=ls2(axiom, mode);
  [raster, offset, counts, dimension]=ls2(axiom, mode);
% =======
%   [raster, offset, counts, dimension, isLeaf, maxiy, leafs]=ls2(axiom, mode);
% >>>>>>> .r201
elseif iscell(axiom)
  raster    = axiom{1};
  dimension = axiom{2};
  offset    = axiom{3};
  if numel(axiom)>3
    leafs = axiom{4};
  end
else
  error('axiom not understood!!!');
end

% <<<<<<< .mine
[a array] = drawtree(raster,dimension,offset,0,1,[0 0 1; 1 0 0],varargin{:});
% drawtreeleafed(raster, leaves, dimension,offset,0,1,[0 0.6 0; 1 0 0],varargin{:});

% =======
% if exist('alsoLeafs', 'var') && alsoLeafs && exist('leafs', 'var')
%   drawtree(raster,dimension,offset,0,1,[0 0 0; 1 0 0], 'fit', leafs);
% else
%   drawtree(raster,dimension,offset,0,1,[0 0 0; 1 0 0], 'fit');
% end
% >>>>>>> .r201
drawnow