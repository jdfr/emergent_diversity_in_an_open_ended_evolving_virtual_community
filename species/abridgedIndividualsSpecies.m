function varargout = abridgedIndividualsSpecies(basedir, prefix, saveIt, varargin)
%abridgedIndividualsSpecies({'rene/dataset/G0.5F1/uno/1=0.001_1', 'rene/dataset/G0.5F1/uno/1=0.001_2', 'rene/dataset/G1F1/uno/1=0.001_1'}, {'E1_', 'E2_', 'MP_', 'NM_'}, true);
%makeAllSpeciations('rene/dataset/G0.5F1/uno/1=0.001_1', 1, 79, 'rene/src');

if iscell(basedir)
  for k=1:numel(basedir)
    abridgedIndividualsSpecies(basedir{k}, prefix, saveIt, varargin{:});
  end
  return
end
if iscell(prefix)
  for k=1:numel(prefix)
    abridgedIndividualsSpecies(basedir, prefix{k}, saveIt, varargin{:});
  end
  return
end

if nargin<4
  f1 = [basedir filesep prefix 'individualsSpecies.mat'];
  f2 = [basedir filesep prefix 'numGroups.mat'];
  if not(exist(f1, 'file')) 
    fprintf('REQUIRED FILE %s not found!!!!\n', f1);
    return
  end
  if not(exist(f2, 'file'))
    fprintf('REQUIRED FILE %s not found!!!!\n', f2);
    return
  end
  individualsSpecies = load(f1);
  if not(isfield(individualsSpecies, 'individualsSpecies'))
    individualsSpecies
    error('jarllll!!!!');
  end
  individualsSpecies = individualsSpecies.individualsSpecies;
  numGroups = load(f2);
else
  individualsSpecies = varargin{1};
  numGroups = varargin{2};
end

if numel(numGroups)~=numel(individualsSpecies)
  error('Both vars must be the same length (%d, %d)!!!!', numel(numGroups), numel(individualsSpecies));
end

individualsSpeciesJust = cell(size(individualsSpecies));
for k=1:numel(individualsSpecies)
  individualsSpeciesJust{k} = individualsSpecies{k}{numGroups{k}};
end

if saveIt
  save([basedir filesep prefix 'individualsSpeciesJust.mat'], 'individualsSpeciesJust');
end

if nargout>0
  varargout{1} = individualsSpeciesJust;
end
