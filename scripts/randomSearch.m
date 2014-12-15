function randomSearch
% main('embedded', 'newsimSeeded', '[G]', 'prob', {[0 0 0 0 0.05 0.05 0.05 0]}, 'alphaF', {1/3}, 'tam_poblacion', {1}, 'variablePopSize', {true}, 'stopIfTooTall', {15000}, 'variableNumSonsMode', {'round#[-0.25 0.25]'});
% main('embedded', 'newsimSeeded', 'G', 'prob', {[0 0.05 0.05 0.05 0.05 0.05 0.05 0]}, 'alphaF', {1/3}, 'tam_poblacion', {1}, 'variablePopSize', {true}, 'stopIfTooTall', {15000}, 'variableNumSonsMode', {'round#[0 0.6]'}, 'alphaG', {0.75}, 'treePlacementParam', {0}, 'treePlacementOffset', {300});
%main('embedded', 'newsimSeeded', 'G', 'prob', {[0 0.05 0 0 0 0.05 0.05 0]}, 'PRESIONES', {0.0001}, 'tam_poblacion', {500}, 'variablePopSize', {false}, 'stopIfTooTall', {15000}, 'alphaG', {0.3}, 'treePlacementParam', {0}, 'treePlacementOffset', {300}, 'mode_op_insercion', {{'+', '-', 'G'}});

%main('embedded', {['resultados' filesep 'T2_A20'], 'newsimSeeded'}, 'G',   'prob', {[0 0.05 0.05 0.05 0.05 0.05 0.05 0]}, 'factorG', {1}, 'alphaF', {1/3}, 'tam_poblacion', {1}, 'initialPos', {0}, 'N', {inf}, 'variablePopSize', {true}, 'stopIfTooTall', {15000}, 'variableNumSonsMode', {'round#[0 0.6]'}, 'alphaG', {0.75}, 'treePlacementParam', {0}, 'treePlacementOffset', {100}, 'stopIfTooTall', {15000}, 'generaciones', {5000}, 'ensayos', {10}, 'plottingFactor', {10}, 'T.angle', {22});
%main('embedded', {['resultados' filesep 'T3_A20'], 'newsimSeeded'}, '[G]', 'prob', {[0 0    0    0    0.05 0.05 0.05 0]}, 'factorG', {1}, 'alphaF', {1/3}, 'tam_poblacion', {1}, 'initialPos', {0}, 'N', {inf}, 'variablePopSize', {true}, 'stopIfTooTall', {15000}, 'variableNumSonsMode', {'round#[0 0.6]'}, 'alphaG', {0.75}, 'treePlacementParam', {0}, 'treePlacementOffset', {100}, 'stopIfTooTall', {15000}, 'generaciones', {5000}, 'ensayos', {10}, 'plottingFactor', {10}, 'T.angle', {22});
%main('embedded', {['resultados' filesep 'T2_A05'], 'newsimSeeded'}, 'G',   'prob', {[0 0.05 0.05 0.05 0.05 0.05 0.05 0]}, 'factorG', {1}, 'alphaF', {1/3}, 'tam_poblacion', {1}, 'initialPos', {0}, 'N', {inf}, 'variablePopSize', {true}, 'stopIfTooTall', {15000}, 'variableNumSonsMode', {'round#[0 0.6]'}, 'alphaG', {0.75}, 'treePlacementParam', {0}, 'treePlacementOffset', {100}, 'stopIfTooTall', {15000}, 'generaciones', {5000}, 'ensayos', {10}, 'plottingFactor', {10});
%main('embedded', {['resultados' filesep 'T3_A05'], 'newsimSeeded'}, '[G]', 'prob', {[0 0    0    0    0.05 0.05 0.05 0]}, 'factorG', {1}, 'alphaF', {1/3}, 'tam_poblacion', {1}, 'initialPos', {0}, 'N', {inf}, 'variablePopSize', {true}, 'stopIfTooTall', {15000}, 'variableNumSonsMode', {'round#[0 0.6]'}, 'alphaG', {0.75}, 'treePlacementParam', {0}, 'treePlacementOffset', {100}, 'stopIfTooTall', {15000}, 'generaciones', {5000}, 'ensayos', {10}, 'plottingFactor', {10});

%main('embedded', {['resultados' filesep 'T2_A20_F10'], 'newsimSeeded'}, 'G',   'prob', {[0 0.05 0.05 0.05 0.05 0.05 0.05 0]}, 'factorG', {1}, 'alphaF', {0.1}, 'tam_poblacion', {1}, 'initialPos', {0}, 'N', {inf}, 'variablePopSize', {true}, 'stopIfTooTall', {15000}, 'variableNumSonsMode', {'round#[0 0.6]'}, 'alphaG', {0.75}, 'treePlacementParam', {0}, 'treePlacementOffset', {100}, 'stopIfTooTall', {15000}, 'generaciones', {5000}, 'ensayos', {10}, 'plottingFactor', {10}, 'T.angle', {22});
%main('embedded', {['resultados' filesep 'T3_A20_F10'], 'newsimSeeded'}, '[G]', 'prob', {[0 0    0    0    0.05 0.05 0.05 0]}, 'factorG', {1}, 'alphaF', {0.1}, 'tam_poblacion', {1}, 'initialPos', {0}, 'N', {inf}, 'variablePopSize', {true}, 'stopIfTooTall', {15000}, 'variableNumSonsMode', {'round#[0 0.6]'}, 'alphaG', {0.75}, 'treePlacementParam', {0}, 'treePlacementOffset', {100}, 'stopIfTooTall', {15000}, 'generaciones', {5000}, 'ensayos', {10}, 'plottingFactor', {10}, 'T.angle', {22});
%main('embedded', {['resultados' filesep 'T2_A05_F10'], 'newsimSeeded'}, 'G',   'prob', {[0 0.05 0.05 0.05 0.05 0.05 0.05 0]}, 'factorG', {1}, 'alphaF', {0.1}, 'tam_poblacion', {1}, 'initialPos', {0}, 'N', {inf}, 'variablePopSize', {true}, 'stopIfTooTall', {15000}, 'variableNumSonsMode', {'round#[0 0.6]'}, 'alphaG', {0.75}, 'treePlacementParam', {0}, 'treePlacementOffset', {100}, 'stopIfTooTall', {15000}, 'generaciones', {5000}, 'ensayos', {10}, 'plottingFactor', {10});
%main('embedded', {['resultados' filesep 'T3_A05_F10'], 'newsimSeeded'}, '[G]', 'prob', {[0 0    0    0    0.05 0.05 0.05 0]}, 'factorG', {1}, 'alphaF', {0.1}, 'tam_poblacion', {1}, 'initialPos', {0}, 'N', {inf}, 'variablePopSize', {true}, 'stopIfTooTall', {15000}, 'variableNumSonsMode', {'round#[0 0.6]'}, 'alphaG', {0.75}, 'treePlacementParam', {0}, 'treePlacementOffset', {100}, 'stopIfTooTall', {15000}, 'generaciones', {5000}, 'ensayos', {10}, 'plottingFactor', {10});

%main('embedded', {['resultados' filesep 'T2_A20_F10G50'], 'newsimSeeded'}, 'G',   'prob', {[0 0.05 0.05 0.05 0.05 0.05 0.05 0]}, 'factorG', {1}, 'alphaF', {0.1}, 'tam_poblacion', {1}, 'initialPos', {0}, 'N', {inf}, 'variablePopSize', {true}, 'stopIfTooTall', {15000}, 'variableNumSonsMode', {'round#[0 0.6]'}, 'alphaG', {0.5}, 'treePlacementParam', {0}, 'treePlacementOffset', {100}, 'stopIfTooTall', {15000}, 'generaciones', {5000}, 'ensayos', {10}, 'plottingFactor', {10}, 'T.angle', {22});
%main('embedded', {['resultados' filesep 'T2_A20_F15G50'], 'newsimSeeded'}, 'G',   'prob', {[0 0.05 0.05 0.05 0.05 0.05 0.05 0]}, 'factorG', {1}, 'alphaF', {1/15}, 'tam_poblacion', {1}, 'initialPos', {0}, 'N', {inf}, 'variablePopSize', {true}, 'stopIfTooTall', {15000}, 'variableNumSonsMode', {'round#[0 0.6]'}, 'alphaG', {0.5}, 'treePlacementParam', {0}, 'treePlacementOffset', {100}, 'stopIfTooTall', {15000}, 'generaciones', {5000}, 'ensayos', {10}, 'plottingFactor', {10}, 'T.angle', {22});
%main('embedded', {['resultados' filesep 'T2_A20_F20G50'], 'newsimSeeded'}, 'G',   'prob', {[0 0.05 0.05 0.05 0.05 0.05 0.05 0]}, 'factorG', {1}, 'alphaF', {1/20}, 'tam_poblacion', {1}, 'initialPos', {0}, 'N', {inf}, 'variablePopSize', {true}, 'stopIfTooTall', {15000}, 'variableNumSonsMode', {'round#[0 0.6]'}, 'alphaG', {0.5}, 'treePlacementParam', {0}, 'treePlacementOffset', {100}, 'stopIfTooTall', {15000}, 'generaciones', {5000}, 'ensayos', {10}, 'plottingFactor', {10}, 'T.angle', {22});

%555555555555555555555555
% msub = 'T2_A20';     subsub = '2009_Jun_23_14_55_16'; subsubn='1=0.0001_1'; main('embedded', {['resultados' filesep msub], 'continue'}, [subsub filesep subsubn], ['continua_guai_' subsub], 0.0001, 10000, 'stopIfTooTall', {5000}, 'ensayos', {1});
% msub = 'T2_A20_F10'; subsub = '2009_Jun_24_15_27_25'; subsubn='1=0.0001_1'; main('embedded', {['resultados' filesep msub], 'continue'}, [subsub filesep subsubn], ['continua_guai_' subsub], 0.0001, 10000, 'stopIfTooTall', {5000}, 'ensayos', {1});
% msub = 'T2_A20';     subsub = '2009_Jun_24_12_31_04'; subsubn='1=0.0001_1'; main('embedded', {['resultados' filesep msub], 'continue'}, [subsub filesep subsubn], ['continua_guai_' subsub], 0.0001, 10000, 'stopIfTooTall', {5000}, 'ensayos', {1});
% msub = 'T2_A20';     subsub = '2009_Jun_24_12_31_01'; subsubn='1=0.0001_1'; main('embedded', {['resultados' filesep msub], 'continue'}, [subsub filesep subsubn], ['continua_guai_' subsub], 0.0001, 10000, 'stopIfTooTall', {5000}, 'ensayos', {1});
% msub = 'T3_A20';     subsub = '2009_Jun_24_12_31_10'; subsubn='1=0.0001_1'; main('embedded', {['resultados' filesep msub], 'continue'}, [subsub filesep subsubn], ['continua_guai_' subsub], 0.0001, 10000, 'stopIfTooTall', {5000}, 'ensayos', {1});
% msub = 'T3_A20';     subsub = '2009_Jun_24_12_31_25'; subsubn='1=0.0001_4'; main('embedded', {['resultados' filesep msub], 'continue'}, [subsub filesep subsubn], ['continua_guai_' subsub], 0.0001, 10000, 'stopIfTooTall', {5000}, 'ensayos', {1});
% 
% main('embedded', {['resultados' filesep 'T2_A20'], 'newsimSeeded'}, 'stopIfTooTall', {5000});
% 
% main('embedded', {['resultados' filesep 'T3_A20'], 'newsimSeeded'}, '[G]', 'prob', {[0 0    0    0    0.05 0.05 0.05 0]}, 'factorG', {1}, 'alphaF', {1/3}, 'tam_poblacion', {1}, 'initialPos', {0}, 'N', {inf}, 'variablePopSize', {true}, 'stopIfTooTall', {15000}, 'variableNumSonsMode', {'round#[0 0.6]'}, 'alphaG', {0.75}, 'treePlacementParam', {0}, 'treePlacementOffset', {100}, 'stopIfTooTall', {15000}, 'generaciones', {5000}, 'ensayos', {20}, 'plottingFactor', {10}, 'T.angle', {22}, 'tooLongEvTime', {1500}, 'stopIfTooTall', {5000});
% main('embedded', {['resultados' filesep 'T2_A20'], 'newsimSeeded'}, 'G', 'prob',   {[0 0.05 0.05 0.05 0.05 0.05 0.05 0]}, 'factorG', {1}, 'alphaF', {1/3}, 'tam_poblacion', {1}, 'initialPos', {0}, 'N', {inf}, 'variablePopSize', {true}, 'stopIfTooTall', {15000}, 'variableNumSonsMode', {'round#[0 0.6]'}, 'alphaG', {0.75}, 'treePlacementParam', {0}, 'treePlacementOffset', {100}, 'stopIfTooTall', {15000}, 'generaciones', {5000}, 'ensayos', {20}, 'plottingFactor', {10}, 'T.angle', {22}, 'tooLongEvTime', {1500}, 'stopIfTooTall', {5000});

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %T2 con DA+I+E+A
% main('embedded', {'resultados/T2_DA+I+E+A', 'newsimSeeded'}, 'G',   'prob',   {[0 0.05 0.05 0 0.00 0.05 0.05 0]},                                                                  'factorG', {1}, 'alphaF', {1}, 'tam_poblacion', {1}, 'initialPos', {0}, 'N', {inf}, 'variablePopSize', {true}, 'variableNumSonsMode', {'round#[0 0.5]'}, 'alphaG', {0.75}, 'treePlacementParam', {0}, 'treePlacementOffset', {100}, 'generaciones', {5000}, 'ensayos', {5}, 'plottingFactor', {10}, 'T.angle', {22}, 'tooLongEvTime', {1500}, 'stopIfTooTall', {5000});
% %T2 con I+E+A
% main('embedded', {'resultados/T2_I+E+A',    'newsimSeeded'}, 'G',   'prob',   {[0 0.05 0.00 0 0.00 0.05 0.05 0]},                                                                  'factorG', {1}, 'alphaF', {1}, 'tam_poblacion', {1}, 'initialPos', {0}, 'N', {inf}, 'variablePopSize', {true}, 'variableNumSonsMode', {'round#[0 0.5]'}, 'alphaG', {0.75}, 'treePlacementParam', {0}, 'treePlacementOffset', {100}, 'generaciones', {5000}, 'ensayos', {5}, 'plottingFactor', {10}, 'T.angle', {22}, 'tooLongEvTime', {1500}, 'stopIfTooTall', {5000});
% %T3 con DS+I+E+A
% main('embedded', {'resultados/T3_DS+I+E+A', 'newsimSeeded'}, '[G]', 'prob',   {[0 0.05 0.00 0 0.05 0.05 0.05 0]}, 'mode_op_insercion', {{'G', '+', '-'}},                          'factorG', {1}, 'alphaF', {1}, 'tam_poblacion', {1}, 'initialPos', {0}, 'N', {inf}, 'variablePopSize', {true}, 'variableNumSonsMode', {'round#[0 0.5]'}, 'alphaG', {0.75}, 'treePlacementParam', {0}, 'treePlacementOffset', {100}, 'generaciones', {5000}, 'ensayos', {5}, 'plottingFactor', {10}, 'T.angle', {22}, 'tooLongEvTime', {1500}, 'stopIfTooTall', {5000});
% %T3 con Im+E+A
% main('embedded', {'resultados/T3_Im+E+A',   'newsimSeeded'}, '[G]', 'prob',   {[0 0.05 0.00 0 0.00 0.05 0.05 0]}, 'mode_op_insercion', {'one_level'}, 'insertOnlyNested', {false}, 'factorG', {1}, 'alphaF', {1}, 'tam_poblacion', {1}, 'initialPos', {0}, 'N', {inf}, 'variablePopSize', {true}, 'variableNumSonsMode', {'round#[0 0.5]'}, 'alphaG', {0.75}, 'treePlacementParam', {0}, 'treePlacementOffset', {100}, 'generaciones', {5000}, 'ensayos', {5}, 'plottingFactor', {10}, 'T.angle', {22}, 'tooLongEvTime', {1500}, 'stopIfTooTall', {5000});

%%%%%%%%%%%%%%%%%%%%%%%%%%%

stableParams = {...
  'tam_poblacion', {1}, 'variablePopSize', {true}, 'N', {10000}, 'initialPos', {5000}, ...
  'variableNumSonsMode', {'round#[0 0.6]'}, 'factorG', {1} 'treePlacementParam', {0}, ...
  'plottingFactor', {10}, 'ensayos', {1}, 'generaciones', {500}, 'stopIfTooTall', {15000}, 'tooLongEvTime', {1000}, ...
  'variableNumSonsMode', {'round#[0 0.6]'} ...
  };

otherParams = {@angle, @treePlacementOffset, @alphaF, @alphaG};

for k=1:1000
  pop = poptype();
  op  = cellfun(@(x)x(), otherParams, 'uniformoutput', false);
  op  = [op{:}];
  try
    params = main('embedded', pop{:}, stableParams{:}, op{:});
    makeAllGraphicsRec(params.nomdir, 'nonono.png');
  catch ME
    fprintf('Uuuups, ERROR!!!:\n%s', showError(ME));
  end
end

function vnsm = variableNumSonsMode
%different distributions to shift number of sons
vnsm = {'variableNumSonsMode', {['round#' any2str((rand*0.5-0.25)+[0 ((1+rand)/2)])]}};

function a = angle
%angles between 1º and 45º
a = {'T.angle', {((45-1)*rand+1)}};

function tpo = treePlacementOffset
tpo = {'treePlacementOffset', {500*rand+10}};

function ag = alphaG
ag = {'alphaG', {(1-0.1)*rand+0.1}};

function af = alphaF
af = {'alphaF', {(0.5-0.05)*rand+0.05}};

function pt = poptype
if rand<0.5
  pt = {'newsimSeeded', 'G', 'prob', {[0 0.05 0.05 0.05 0.05 0.05 0.05 0]}};
else
  pt = {'newsimSeeded', '[G]', 'prob', {[0 0 0 0 0.05 0.05 0.05 0]}};
end