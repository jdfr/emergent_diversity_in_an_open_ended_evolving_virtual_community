function processSpeciation(data, poph)

fig = figure; %#ok<NASGU>

colors = rand(numel(data.species.generations), 3)*0.8;

pointArgs       = {'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 3}; %#ok<NASGU>
lineArgs        = {'LineWidth', 2, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 3};
connectLineArgs = {'LineStyle', '--', 'Color', 'k'};

dots  = cell(size(data.species.generations)); %#ok<NASGU>
lines = cell(size(data.species.generations));

str = '';
for k=1:numel(data.species.generations)
  if numel(str)>0
    fprintf(repmat('\b', 1, numel(str)));
  end
  str = sprintf('Amos por la especie %d/%d, pisha', k,numel(data.species.generations));
  fprintf(str);
  %dots{k}  = line(data.species.generations{k}([1 end]), data.species.complexities{k}([1 end]), pointArgs{:});
  lines{k} = line(data.species.generations{k}, data.species.complexities{k}, 'Color', colors(k,:), ...
    'ButtonDownFcn', makeCallback(fig, k, data.species.generations{k}, data.species.complexities{k}, data.species.groups{k}, data, poph), ...
    lineArgs{:});
  for m = 1:numel(data.species.tree{k})
    hija = data.species.tree{k}(m);
    line( ...
        [data.species.generations{k}(end),  data.species.generations{hija}(1)], ...
        [data.species.complexities{k}(end), data.species.complexities{hija}(1)], ...
        connectLineArgs{:});
  end
end
if numel(str)>0
  fprintf(repmat('\b', 1, numel(str)));
end

    
grid on;



function cb = makeCallback(fig, species, generations, complexities, groups, data, poph)
cb = @(src, event)doCallback(fig, species, generations, complexities, groups, data, poph);

function doCallback(fig, species, generations, complexities, groups, data, poph) %#ok<INUSD>
ax = get(fig, 'CurrentAxes');
currentPoint= get(ax, 'CurrentPoint'); 
cx            = currentPoint(1,1);
cy            = currentPoint(1,2);

errs          = abs([generations(:)-cx, complexities(:)-cy]);
[ignore idx]  = min(errs(:,1)+errs(:,2));%#ok<ASGLU>

generation     = generations(idx);
group          = groups(idx);
complexity     = complexities(idx);
inds           = data.groups.individualsByProto{generation}{group};
numIndividuals = data.groups.numIndividuals{generation}(group);
genomes        = unique(inds.genomes);
genomes        = cellfun(@sanitizeString, genomes, 'uniformoutput', false);
genomes        = unique(genomes);
numGenomes     = numel(genomes);
genomesSize    = cellfun(@numel, genomes);
[minsz idxmin] = min(genomesSize);

numProtos = numel(data.groups.individualsByProto{generation}{group}.individuals);
genomedata = data.groups.individualsByProto{generation}{group}.genomes{1};
individualToSearch = data.groups.individualsByProto{generation}{group}.individuals(1); %data.groups.individualsByGroup{generation}{group}(1);
idxpoph = find(poph.generation==generation);
idxpoph = idxpoph(poph.individual(idxpoph)==individualToSearch);
genomepoph = poph.genome{idxpoph};
idxidxpoph = poph.idx(idxpoph);
idxpophtree = find(poph.tree.idxD==idxidxpoph);
changepoph  = poph.tree.change{idxpophtree};
idxpophancestor = poph.tree.idxA(idxpophtree);
if idxpophancestor>0
  gAncestor = poph.tree.gA(idxpophtree);
  iAncestor = poph.tree.iA(idxpophtree);
  indancestor = find((poph.generation==gAncestor) & (poph.individual==iAncestor) );%find(poph.idx==idxpophancestor);
  genomeAncestor = poph.genome{indancestor};
  strAncestor = sprintf('\nAncestor''s generation: %d, individual %d, genome: %s', gAncestor, iAncestor, genomeAncestor);
else
  strAncestor = '';
end

% figure;
% [rasterpoph, offsetpoph, countspoph, dimensionpoph]=ls2(genomepoph, poph.ACIParams.T);
% drawtree(rasterpoph,dimensionpoph,offsetpoph,0);
% title('THIS IS FROM POPH');

titleTXT = sprintf('You have clicked on species number %d, group %d in generation %d, complexity %s\ndimensions: %s,  nbranches: %d, npixels: %d', species, group, generation, mat2str(complexity), mat2str(inds.dimensions), inds.counts(1), numel(inds.raster{1}{1}));
labelTXT = sprintf('%d individuals in that group and generation, %d equal to the prototype (with %d distinct genomes), minimal purified genome (length %d):\n%s\nFor the first individual in Gema''s speciation for this species:\ngeneration %d, individual %d, CHANGE IN GENOME: %s\nGENOME IN POPH: %s\nGENOME IN DATA: %s%s', numIndividuals, numProtos, numGenomes, minsz, genomes{idxmin}, generation, individualToSearch, changepoph, genomepoph, genomedata, strAncestor);

title(ax,  titleTXT);
xlabel(ax, labelTXT);

TRUE_SCALE = false;

sel_typ = get(fig,'SelectionType');
switch sel_typ 
 %case 'normal'; case 'extend'
 case 'extend'
   fprintf('Individuals for species %d at generation %d: %s\n', species, generation, mat2str(data.groups.individualsByGroup{generation}{group}));
 case 'alt'
   figure;
   drawtree(inds.raster{1},inds.dimensions,inds.offset,0,1,[0 0 0], ~TRUE_SCALE);
   title(titleTXT);
   xlabel(labelTXT);
%    figure;
%    proto = data.groups.protos{generation}{group};
%    drawtree(proto,[max(proto{1}), max(proto{2})], 1,0);
end


