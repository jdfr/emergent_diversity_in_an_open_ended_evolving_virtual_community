function makeDiversityPlots(arg1, basedirClassification, toFile, override, paperSize)
%make phase plots about diversity (number of species)

if ~exist('paperSize', 'var')
  paperSize = 1000;
end

if ~exist('toFile', 'var')
  toFile = true;
end

if ~exist('override', 'var')
  override = false;
end

speciesMAT = [basedirClassification, filesep 'speciesFacts.mat'];

if (~override) && exist(speciesMAT, 'file')
  load(speciesMAT);
else
      if ischar(arg1)
        poph = loadSimulation(arg1);
      else
        poph = arg1;
        clear arg1;
      end

      d = dir(basedirClassification);

      minGen = min(poph.generation);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %get numSpecies
      toUse = find(arrayfun(@(x)~isempty(strfind(x.name, 'numSpecies')), d));
      if isempty(toUse)
        error('directory <%s> does''nt harbor any file named numSpecies*!!!',  basedirClassification);
      elseif numel(toUse)>1
        fprintf('Several files numSpecies* in directory <%s>!!! Using <%s>!!!!\n', basedirClassification, d(toUse(1)).name);
      end
      fname   = d(toUse(1)).name;
      varname = fname(1:end-4);
      load([basedirClassification filesep fname], varname);
      eval(sprintf('numSpecies = %s; clear %s;', varname, varname));

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %get individualsSpecies
      toUse = find(arrayfun(@(x)~isempty(strfind(x.name, 'individualsSpecies')), d));
      if isempty(toUse)
        error('directory <%s> does''nt harbor any file named individualsSpecies*!!!',  basedirClassification);
      elseif numel(toUse)>1
        fprintf('Several files individualsSpecies* in directory <%s>!!! Using <%s>!!!!\n', basedirClassification, d(toUse(1)).name);
      end
      fname   = d(toUse(1)).name;
      varname = fname(1:end-4);
      load([basedirClassification filesep fname], varname);
      eval(sprintf('individualsSpecies = %s; clear %s;', varname, varname));

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %compute number of species by generation, total population by generation
      numSpeciesByGen = cellfun(@(x)sum(x~=0), numSpecies); %#ok<USENS>

      longGenome         = cellfun(@numel, poph.genome);

      populationByGen    = zeros(size(numSpeciesByGen));
      meanLongGenomeByGen= zeros(size(numSpeciesByGen));
      meanHeightByGen    = zeros(size(numSpeciesByGen));
      meanWidthByGen     = zeros(size(numSpeciesByGen));
      meanNBranchesByGen = zeros(size(numSpeciesByGen));
      meanFitnessByGen   = zeros(size(numSpeciesByGen));
      
      maxLongGenomeByGen= zeros(size(numSpeciesByGen));
      maxHeightByGen    = zeros(size(numSpeciesByGen));
      maxWidthByGen     = zeros(size(numSpeciesByGen));
      maxNBranchesByGen = zeros(size(numSpeciesByGen));
      maxFitnessByGen   = zeros(size(numSpeciesByGen));

      currentGen = minGen;

      for k=1:numel(populationByGen)
        thisGen = find(poph.generation==currentGen);
        populationByGen(k)     = numel(thisGen);
        meanLongGenomeByGen(k) = mean(longGenome(thisGen));
        meanHeightByGen(k)     = mean(poph.height(thisGen));
        meanWidthByGen(k)      = mean(poph.width(thisGen));
        meanNBranchesByGen(k)  = mean(poph.nbranches(thisGen));
        maxLongGenomeByGen(k)  =  max(longGenome(thisGen));
        maxHeightByGen(k)      =  max(poph.height(thisGen));
        maxWidthByGen(k)       =  max(poph.width(thisGen));
        maxNBranchesByGen(k)   =  max(poph.nbranches(thisGen));
        auxfit                 = poph.fitness(thisGen);
        auxfit(auxfit<0)       = 0;
        meanFitnessByGen(k)    = mean(auxfit);
        maxFitnessByGen(k)     =  max(auxfit);
        currentGen             = currentGen+1;
        clear auxfit thisGen;
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %amount of individuals and means per species (across
      %all its timesteps)


      speciesNumIndividuals = zeros(size(numSpecies{end}));

      speciesMeanLongGenome = zeros(size(numSpecies{end}));
      speciesMeanHeight     = zeros(size(numSpecies{end}));
      speciesMeanWidth      = zeros(size(numSpecies{end}));
      speciesMeanFitness    = zeros(size(numSpecies{end}));
      speciesMeanNBranches  = zeros(size(numSpecies{end}));

      speciesStartT         = zeros(size(numSpecies{end}));
      speciesEndT           = zeros(size(numSpecies{end}));

      speciesIndividual = zeros(size(longGenome));
      currentGen = minGen;
      for k=1:numel(individualsSpecies) %#ok<USENS>
        thisGen = find(poph.generation==currentGen);
        for m=1:numel(individualsSpecies{k})
          speciesIndividual(thisGen(individualsSpecies{k}{m})) = m;
        end
        maxGen                 = currentGen;
        currentGen             = currentGen+1;
      end

      generations = minGen:maxGen;
      
      for k=1:numel(speciesNumIndividuals)
        indexes = find(speciesIndividual==k);
        speciesNumIndividuals(k) = numel(indexes);
        speciesMeanLongGenome(k) = sum(longGenome(indexes));
        speciesMeanHeight(k)     = sum(poph.height(indexes));
        speciesMeanWidth(k)      = sum(poph.width(indexes));
        auxfit                   = poph.fitness(indexes);
        auxfit(auxfit<0)         = 0;
        speciesMeanFitness(k)    = sum(auxfit);
        clear auxfit;
        speciesMeanNBranches(k)  = sum(poph.nbranches(indexes));
        auxgen                   = poph.generation(indexes);
        speciesStartT(k)         = min(auxgen);
        speciesEndT(k)           = max(auxgen);
        clear auxgen;
      end
      speciesSpanT          = speciesEndT - speciesStartT + 1;
      speciesMeanLongGenome = speciesMeanLongGenome./speciesNumIndividuals;
      speciesMeanHeight     = speciesMeanHeight./speciesNumIndividuals;
      speciesMeanWidth      = speciesMeanWidth./speciesNumIndividuals;
      speciesMeanFitness    = speciesMeanFitness./speciesNumIndividuals;
      speciesMeanNBranches  = speciesMeanNBranches./speciesNumIndividuals;
      
      save(speciesMAT, ...
        'numSpecies', 'individualsSpecies', 'numSpeciesByGen', 'populationByGen', ...
        'speciesNumIndividuals', 'speciesMeanLongGenome', 'speciesMeanHeight', ...
        'speciesMeanWidth', 'speciesMeanFitness', 'speciesMeanNBranches', ...
        'speciesStartT', 'speciesEndT', 'speciesSpanT', ...
        'meanLongGenomeByGen', 'meanHeightByGen', 'meanWidthByGen', ...
        'meanNBranchesByGen', 'meanFitnessByGen', ...
        'maxLongGenomeByGen', 'maxHeightByGen', 'maxWidthByGen', ...
        'maxNBranchesByGen', 'maxFitnessByGen', 'generations');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lineargsScatter = {'LineStyle', 'none', 'Marker', 'd'};%, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'};
lineargsPhase   = {'LineStyle', 'none', 'Marker', 'd'};%, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'};
% lineargsPhase   = {};%, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'};

if toFile
  figargs = {'Visible', 'off', 'Units', 'pixels', 'PaperUnits', 'points', 'PaperPosition', [0 0 paperSize paperSize]};
else
  figargs = {};
end

cmap = 'cool';

% make3DPhasePlotColor(basedirClassification, figargs, lineargsPhase, cmap, toFile, ...
%   sprintf('phase plot: diversity by generation VS. population size VS. mean Branch Count\nColor axis denotes generation'), ...
%   'diversityVSpopSizeVSmeanBranchCount', ...
%   numSpeciesByGen,    'number of distinct species by generation', ...
%   populationByGen,    'population size', ...
%   meanNBranchesByGen, 'mean Branch Count', ...
%   generations);

makePhasePlot = @makePhasePlotColor;
%makePhasePlot = @makePhasePlotPlain;

coloraxis = '\nColor axis denotes generation';

makePhasePlot(basedirClassification, figargs, lineargsPhase, cmap, toFile, ...
  sprintf(['phase plot: diversity by generation VS. population size' coloraxis]), ...
  'diversityVSpopSize', ...
  numSpeciesByGen, 'number of distinct species by generation', ...
  populationByGen, 'population size', ...
  generations);

templatePhasePlot = @(name, values) makePhasePlot(basedirClassification, figargs, lineargsPhase, cmap, toFile, ...
  sprintf(['phase plot: diversity by generation VS. ' name ' by generation' coloraxis]), ...
  ['diversityVS' name(name~=' ')], ...
  numSpeciesByGen, 'number of distinct species by generation', ...
  values,          [name ' by generation'], ...
  generations);

templatePhasePlot('mean Height',        meanHeightByGen);
templatePhasePlot('mean Width',         meanWidthByGen);
templatePhasePlot('mean Branch Count',  meanNBranchesByGen);
templatePhasePlot('mean Fitness',       meanFitnessByGen);
templatePhasePlot('mean Genome Length', meanLongGenomeByGen);
templatePhasePlot( 'max Height',        maxHeightByGen);
templatePhasePlot( 'max Width',         maxWidthByGen);
templatePhasePlot( 'max Branch Count',  maxNBranchesByGen);
templatePhasePlot( 'max Fitness',       maxFitnessByGen);
templatePhasePlot( 'max Genome Length', maxLongGenomeByGen);


coloraxis = '\nColor axis denotes species'' first generation';

templateScatterPlot = @(name, values) makeScatterPlot(basedirClassification, figargs, lineargsScatter, toFile, ...
  sprintf(['scatter plot: for each species, mean genome length VS. species'' amount of individuals' coloraxis]), ...
  [name(name~=' ') 'VSNumIndiv'], ...
  values,                [name ' (across all timesteps)'], ...
  speciesNumIndividuals, 'amount of individuals (across all timesteps)', ...
  speciesStartT);

templateScatterPlot('mean Genome Length',     speciesMeanLongGenome);
templateScatterPlot('mean Tree Height',       speciesMeanHeight);
templateScatterPlot('mean Tree Width',        speciesMeanWidth);
templateScatterPlot('mean Fitness',           speciesMeanFitness);
templateScatterPlot('mean Branch Count',      speciesMeanNBranches);
templateScatterPlot('mean Species Time Span', speciesSpanT);

makeScatterPlot(basedirClassification, figargs, lineargsScatter, toFile, ...
  sprintf(['scatter plot: for each species, mean genome length VS. species'' amount of individuals' coloraxis]), ...
  'meanLongGenomeVSNumIndiv', ...
  speciesMeanLongGenome, 'mean genome length (across all timesteps)', ...
  speciesNumIndividuals, 'amount of individuals (across all timesteps)', ...
  speciesStartT);

makeScatterPlot(basedirClassification, figargs, lineargsScatter, toFile, ...
  sprintf(['scatter plot: for each species, mean tree height VS. species'' amount of individuals' coloraxis]), ...
  'meanHeightVSNumIndiv', ...
  speciesMeanHeight,     'mean tree height (across all timesteps)', ...
  speciesNumIndividuals, 'amount of individuals (across all timesteps)', ...
  speciesStartT);

makeScatterPlot(basedirClassification, figargs, lineargsScatter, toFile, ...
  sprintf(['scatter plot: for each species, mean tree width VS. species'' amount of individuals' coloraxis]), ...
  'meanWidthVSNumIndiv', ...
  speciesMeanWidth,      'mean tree width (across all timesteps)', ...
  speciesNumIndividuals, 'amount of individuals (across all timesteps)', ...
  speciesStartT);

makeScatterPlot(basedirClassification, figargs, lineargsScatter, toFile, ...
  sprintf(['scatter plot: for each species, mean tree fitness VS. species'' amount of individuals' coloraxis]), ...
  'meanFitnessVSNumIndiv', ...
  speciesMeanFitness,    'mean tree fitness (across all timesteps)', ...
  speciesNumIndividuals, 'amount of individuals (across all timesteps)', ...
  speciesStartT);

makeScatterPlot(basedirClassification, figargs, lineargsScatter, toFile, ...
  sprintf(['scatter plot: for each species, mean amount of branches VS. species'' amount of individuals' coloraxis]), ...
  'meanNBranchesVSNumIndiv', ...
  speciesMeanNBranches,  'mean amount of branches (across all timesteps)', ...
  speciesNumIndividuals, 'amount of individuals (across all timesteps)', ...
  speciesStartT);

makeScatterPlot(basedirClassification, figargs, lineargsScatter, toFile, ...
  sprintf(['scatter plot: for each species, time span VS. species'' amount of individuals' coloraxis]), ...
  'timeSpanVSNumIndiv', ...
  speciesSpanT,          'species'' time span (across all timesteps)', ...
  speciesNumIndividuals, 'amount of individuals (across all timesteps)', ...
  speciesStartT);

function makePhasePlotPlain(basedir, figargs, lineargs, tofile, titulo, filename, x, xlab, y, ylab, varargin) %#ok<DEFNU>
fig = figure(figargs{:});
line(x,y, lineargs{:});
grid on; xlabel(xlab); ylabel(ylab); title(titulo);
if tofile
  drawnow;
  saveas(fig,[basedir filesep 'species_' filename '.png'],'png');
  close(fig);
end

function makePhasePlotColor(basedir, figargs, lineargs, cmap, tofile, titulo, filename, x, xlab, y, ylab, generations) %#ok<DEFNU>
fig = figure(figargs{:});
cols = generations-min(generations)+1;
cmap3=feval(cmap, (max(cols)));
set(gcf, 'Colormap', cmap3);
hc = colorbar('Location', 'EastOutside');%, 'YTick', numgens, 'YTickLabel', strgens);
freezeColors(hc);

line(x,y, 'Color', 'k');
%for k=1:numel(x)
for k=numel(x):-1:1
  color = cmap3(cols(k),:);
  line(x(k),y(k), lineargs{:}, 'MarkerFaceColor', color, 'MarkerEdgeColor', color);
end
grid on; xlabel(xlab); ylabel(ylab); title(titulo);
if tofile
  drawnow;
  saveas(fig,[basedir filesep 'species_' filename '.png'],'png');
  close(fig);
end

function makeScatterPlot(basedir, figargs, lineargs, tofile, titulo, filename, x, xlab, y, ylab, tStart) %#ok<DEFNU>
fig = figure(figargs{:});
cols = tStart-min(tStart)+1;
cmap3=jet(max(cols));
% [uniquecols idxs] = unique(cols);
% tims              = arrayfun(@mat2str, tStart(idxs), 'uniformoutput', false);
colormap(cmap3);
hc = colorbar('Location', 'EastOutside');%, 'YTick', uniquecols, 'YTickLabel', tims);
freezeColors(hc);

for k=1:numel(x)
%for k=numel(x):-1:1
  color = cmap3(cols(k),:);
  line(x(k),y(k), lineargs{:}, 'MarkerFaceColor', color, 'MarkerEdgeColor', color);
end
grid on; xlabel(xlab); ylabel(ylab); title(titulo);
if tofile
  drawnow;
  saveas(fig,[basedir filesep 'species_' filename '.png'],'png');
  close(fig);
end

function make3DPhasePlotColor(basedir, figargs, lineargs, cmap, tofile, titulo, filename, x, xlab, y, ylab, z, zlab, generations) %#ok<DEFNU>
fig = figure(figargs{:});
cols = generations-min(generations)+1;
cmap3=feval(cmap, (max(cols)));
set(gcf, 'Colormap', cmap3);
hc = colorbar('Location', 'EastOutside');%, 'YTick', numgens, 'YTickLabel', strgens);
freezeColors(hc);

line(x,y,z, 'Color', 'k');
%for k=1:numel(x)
for k=numel(x):-1:1
  color = cmap3(cols(k),:);
  line(x(k),y(k),z(k), lineargs{:}, 'MarkerFaceColor', color, 'MarkerEdgeColor', color);
end
grid on; xlabel(xlab); ylabel(ylab); zlabel(zlab); title(titulo);
view(3);
if tofile
  drawnow;
  saveas(fig,[basedir filesep 'species_' filename '.png'],'png');
  close(fig);
end

