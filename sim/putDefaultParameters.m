function params = putDefaultParameters(params)

if ~isfield(params, 'randmethod')
    params.randmethod = 'twister';
end

if ~isfield(params, 'factorG');
  params.factorG = 10;
end

if ~isfield(params, 'initialGeneration')
  params.initialGeneration = 0;
end

if ~isfield(params, 'N')
  %mimic the assignment in main.m
  params.N = params.tam_poblacion*15; if params.N<3000; params.N = 3000; end
end

if ~isfield(params, 'negYThreshold')
  params.negYThreshold = inf;
end

if ~isfield(params, 'treePlacement')
  params.treePlacement = 'random';
end

if ~isfield(params, 'tooTallArePlaced')
  params.tooTallArePlaced = true;
end

if ~isfield(params, 'treePlacementOffset')
 params.treePlacementOffset = 0;
end

if ~isfield(params, 'variablePopSize')
  params.variablePopSize = false;
end

if ~isfield(params, 'alphaF')
  params.alphaF = 0.5;
end

if ~isfield(params, 'plotDendroTree')
  params.plotDendroTree = false;
end

if ~isfield(params, 'plottingSample') 
  params.plottingSample = false;
end

if ~isfield(params, 'variableNumSonsMode')
  params.variableNumSonsMode = 'mixed';
end

if isfield(params, 'op_insercion_one_level') && (~isfield(params, 'mode_op_insercion'))
  if params.op_insercion_one_level
    params.mode_op_insercion = 'one_level';
  else
    params.mode_op_insercion = 'normal';
  end
end

if (~isfield(params, 'mode_op_insercion'))
  params.mode_op_insercion = 'normal';
end

if ~isfield(params, 'op_eliminacion_extended')
  params.op_eliminacion_extended = false;
end

if ~isfield(params, 'fitnessFunction')
  params.fitnessFunction = 'fitness2';
end

if ~isfield(params, 'plotSampleMode')
  params.plotSampleMode = 'firstIndexes';
end

if ~isfield(params, 'mutationOrder')
  params.mutationOrder = 'normal';
end

if ~isfield(params, 'insertOnlyNested')
  params.insertOnlyNested = false;
end

if ~isfield(params, 'stopIfTooTall')
  params.stopIfTooTall = [];
end

if ~isfield(params, 'selectSymbolForMutation')
  params.selectSymbolForMutation = 'bySymbol';
end

if ~isfield(params, 'initialPos')
  params.initialPos = [];
end

if ~isfield(params, 'plottingFactor')
  params.plottingFactor = 1;
end

if ~isfield(params, 'tooLongEvTime')
  params.tooLongEvTime = inf;
end

if ~isfield(params, 'randomPlacementForInitialGeneration')
  params.randomPlacementForInitialGeneration = true;
end

% if ~isfield(params, 'MAINSUB')
%   params.MAINSUB = 'resultados';
% end

if ~isfield(params, 'allowPrefixChangeDir')
  params.allowPrefixChangeDir = true;
end

if not(isfield(params, 'preserveUniqueTree'))
  params.preserveUniqueTree = true;
end

if not(isfield(params, 'correc'))
  params.correc = struct;
end

if not(isfield(params.correc, 'zeroHeightOK'))
  params.correc.zeroHeightOK = false;
end

if not(isfield(params, 'save'))
  params.save = struct;
end

if not(isfield(params.save, 'disp'))
  params.save.disp = struct;
end

if not(isfield(params.save.disp, 'gnm'))
  params.save.disp.gnm = false;
end

if not(isfield(params.save.disp, 'mng'))
  params.save.disp.mng = false;
end

if not(isfield(params.save.disp, 'pht'))
  params.save.disp.pht = false;
end

