function data = loadSpeciation(basedir, maxGeneration, matFileForGroupSons)

%this is the file containing a cell array with a cell for each generation.
%In the cell N, there is a subcell. Each array in position "i" in the
%subcell "j" indicates the indexes of the groups in the generation "j+1"
%that descend from the group "i" in the generation "j"
if ~exist('matFileForGroupSons', 'var')
  matFileForGroupSons = 'especies';
end

load([basedir filesep matFileForGroupSons]); groupSons = eval(matFileForGroupSons); clear(matFileForGroupSons);

protos             = cell(maxGeneration, 1);
individualsByProto = cell(maxGeneration, 1);
complexities       = cell(maxGeneration, 1);
species            = cell(maxGeneration, 1);
numIndividuals     = cell(maxGeneration, 1);
individualsByGroup = cell(maxGeneration, 1);


%"species" contains, for each group in each generation, its species
%identificator
species{1} = 1:numel(groupSons{1});
%"speciesTree" contains, for each species, the species that descend from it
speciesTree = cell(size(species{1}));
newSpecies = numel(groupSons{1})+1;
for k=2:numel(groupSons);
  sons = groupSons{k-1};
  numSons = cellfun(@numel, sons);
  speciesParents = species{k-1};
  species{k} = zeros(numel(groupSons{k}), 1);
  
  groupsBeingSons = 1:numel(species{k});
  appearances = arrayfun(@(groupBeingSon)cellfun(@(group)any(group==groupBeingSon), sons), groupsBeingSons, 'uniformoutput', false);
  
  zeroGroups = find(cellfun(@(appears)sum(appears)==0, appearances)); %#ok<EFIND>
  if ~isempty(zeroGroups)
    error('data in %s are corrupted!!!', matFileForGroupSons);
  end
  
  %groups that descend from just one group in the previous generation
  straightGroups = find(cellfun(@(appears)sum(appears)==1, appearances));
  if ~isempty(straightGroups)
    parentsOfStraightGroups = cellfun(@find, appearances(straightGroups));
      %groups whose parent has just one child: continue the species
      parentsOfStraightGroupsWith1Son = numSons(parentsOfStraightGroups)==1;
      if any(parentsOfStraightGroupsWith1Son)
        species{k}(straightGroups(parentsOfStraightGroupsWith1Son)) = speciesParents(parentsOfStraightGroups(parentsOfStraightGroupsWith1Son));
      end
      %groups whose parent has several childs: start new species
      parentsOfStraightGroupsWithMoreSons = numSons(parentsOfStraightGroups)>1;
      if any(parentsOfStraightGroupsWithMoreSons)
        parentsOfStraightGroupsHavingSeveralSons = unique(parentsOfStraightGroups(parentsOfStraightGroupsWithMoreSons));
        numNewSpecies = numel(straightGroups(parentsOfStraightGroupsWithMoreSons));
        species{k}(straightGroups(parentsOfStraightGroupsWithMoreSons)) = newSpecies:(newSpecies+numNewSpecies-1);
        newSpecies = newSpecies+numNewSpecies;
        speciesTree{newSpecies-1} = [];
        speciesParentsWithMoreSons = speciesParents(parentsOfStraightGroupsHavingSeveralSons);
        sonsOfParentsWithMoresons  = sons(parentsOfStraightGroupsHavingSeveralSons);
        for m=1:numel(speciesParentsWithMoreSons)
          theseSpecies = species{k}(sonsOfParentsWithMoresons{m});
          speciesTree{speciesParentsWithMoreSons(m)} = theseSpecies(theseSpecies~=0);
        end
      end
  end
  
  %groups that descend from several groups in the previous generation
  conflictiveGroups = find(cellfun(@(appears)sum(appears)>1, appearances));
  if ~isempty(conflictiveGroups)
    parentsOfConflictiveGroups = cellfun(@find, appearances(conflictiveGroups), 'uniformoutput', false);
    numParentsOfConflictiveGroups = cellfun(@numel, parentsOfConflictiveGroups);
      %start new species
      numNewSpecies = numel(conflictiveGroups);
      species{k}(conflictiveGroups) = newSpecies:(newSpecies+numNewSpecies-1);
      newSpecies = newSpecies+numNewSpecies;
      speciesTree{newSpecies-1} = [];
      for m=1:numel(conflictiveGroups)
        conflictiveSpecies = species{k}(conflictiveGroups(m));
        for n=1:numParentsOfConflictiveGroups(m)
          speciesOfThisParent = speciesParents(parentsOfConflictiveGroups{m}(n));
          speciesTree{speciesOfThisParent} = [ ...
                speciesTree{speciesOfThisParent}; ...
                conflictiveSpecies ...
            ];
        end
      end
  end
    
  
%   for m=1:numel(sons)
%     switch numel(sons{m})
%       case 0
%       case 1
%         species{k}(sons{m}) = speciesParents(m);
%       otherwise
%         numNew              = numel(sons{m});
%         species{k}(sons{m}) = newSpecies:(newSpecies+numNew-1);
%         speciesTree{speciesParents(m)}   = species{k}(sons{m});
%         speciesTree(species{k}(sons{m})) = cell(size(sons{m}));
%         newSpecies          = newSpecies+numNew;
%     end
%   end
end

names = {'numGroups', 'individualsSpecies', 'comLZfigureLongitud', 'proto', 'P'};
loadComplexities = ismember('comLZfigureLongitud', names);
loadProtoGenomes = ismember('individualsSpecies', names);
loadProtos       = ismember('proto', names);

str = '';

for k=1:maxGeneration
  if numel(str)>0
    fprintf(repmat('\b', 1, numel(str)));
  end
  str = sprintf('Amos por la generacion %d, pisha', k);
  fprintf(str);
  
  for m=1:numel(names)
    gname = sprintf('%s%03d', names{m}, k);
    load([basedir filesep gname]);
    eval(sprintf('%s = %s; clear %s;', names{m}, gname, gname));
  end
  
  if loadProtos
    proto = proto{numGroups};
  end
  if loadProtoGenomes
    individualsSpecies = individualsSpecies{numGroups};
  end
  
  if loadComplexities
    complexities{k}       = comLZfigureLongitud;
  end
  if loadProtos
    protos{k}             = proto;
  end
  
  if loadProtos
    individualsByProto{k} = cell(size(proto));
  end    
  
  if loadProtoGenomes
    numIndividuals{k}     = cellfun(@numel, individualsSpecies);
  end    
  
  if loadProtos && loadProtoGenomes
    individualsByGroup{k} = individualsSpecies;
    for m=1:numel(proto)
      candidatisimos = individualsSpecies{m}(cellfun(@(x)isequal(proto{m}, x), P.raster(individualsSpecies{m})));
      candidatosos = struct;
      snames = fieldnames(P);
      for n=1:numel(snames)
        candidatosos.(snames{n}) = P.(snames{n})(candidatisimos(1),:);
      end
      candidatosos.genomes = P.genome(candidatisimos);
      candidatosos.randNs  = P.randN(candidatisimos);
      candidatosos.individuals = candidatisimos;
      individualsByProto{k}{m} = candidatosos;
    end
  end
  
  clear cadidatisimos cadidatosos;
  
  cellfun(@clear, names);
end

if numel(str)>0
  fprintf(repmat('\b', 1, numel(str)));
end

%"especies" contains, for each species, the generations it has been alive,
%the groups adscribed to that species in each corresponding generation, and
%the complexity of the prototype in each of these groups
especiesGroups      = cell(newSpecies-1, 1);
especiesGenerations = cell(newSpecies-1, 1);
especiesComplexity  = cell(newSpecies-1, 1);
for k=1:numel(species)
  speciesk = species{k};
  for m=1:numel(speciesk)
    especiesGenerations{speciesk(m)} = [especiesGenerations{speciesk(m)}; k];
    especiesGroups{speciesk(m)}      = [especiesGroups{speciesk(m)};      m];
    if loadComplexities && (k<=numel(complexities))
      especiesComplexity{speciesk(m)}  = [especiesComplexity{speciesk(m)};  complexities{k}(m)];
    else
      especiesComplexity{speciesk(m)}  = [especiesComplexity{speciesk(m)};  NaN];
    end
  end
end

data = struct('groups', [], 'species', []);
data.species = struct('tree', {speciesTree}, 'generations', {especiesGenerations}, 'groups', {especiesGroups}, 'complexities', {especiesComplexity});
data.groups  = struct('protos', {protos}, 'complexities', {complexities}, 'individualsByProto', {individualsByProto}, 'numIndividuals', {numIndividuals}, 'individualsByGroup', {individualsByGroup});
