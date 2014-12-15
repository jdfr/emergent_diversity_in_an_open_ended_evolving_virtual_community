function exhaustiveSearchEdit(Gvalues, Fvalues, ensayos)

stableParams = {...
  'tam_poblacion', {1}, 'variablePopSize', {true}, 'N', {inf}, 'initialPos', {5000}, ...
  'factorG', {1} 'treePlacementParam', {0}, 'treePlacementOffset', {100}, ...
  'plottingFactor', {10}, 'ensayos', {ensayos}, 'generaciones', {500}, 'stopIfTooTall', {5000}, 'tooLongEvTime', {1500}, ...
  'variableNumSonsMode', {'round#[0 0.6]'}, 'T.angle', {22}, ...
  'prob', {[0 0.05 0 0 0 0.05 0.05 0]}, ...
  };

for m=1:numel(Gvalues)
  for n=1:numel(Fvalues)
    params = {'alphaG', {Gvalues(m)}, 'alphaF', {Fvalues(n)}};
    subdir = ['resultados' filesep 'systematicEdit' filesep 'G' mat2str(Gvalues(m)) 'F' mat2str(Fvalues(n))];
    for k=1:numel(ensayos)
      main('embedded', {subdir, 'newsimSeeded'}, 'G', stableParams{:}, params{:});
    end
  end
end

% exhaustiveSearch(0.1, 0.1:0.1:0.5, 5); n0002
% exhaustiveSearch(0.1, 0.6:0.1:1.0, 5); n0002
% exhaustiveSearch(0.2, 0.1:0.1:0.5, 5); n0003
% exhaustiveSearch(0.2, 0.6:0.1:1.0, 5); n0003
% exhaustiveSearch(0.3, 0.1:0.1:0.5, 5); n0004
% exhaustiveSearch(0.3, 0.6:0.1:1.0, 5); n0004
% exhaustiveSearch(0.4, 0.1:0.1:0.5, 5); n0005
% exhaustiveSearch(0.4, 0.6:0.1:1.0, 5); n0005
% exhaustiveSearch(0.5, 0.1:0.1:0.5, 5); n0006
% exhaustiveSearch(0.5, 0.6:0.1:1.0, 5); n0006
% exhaustiveSearch(0.6, 0.1:0.1:0.5, 5); n0007
% exhaustiveSearch(0.6, 0.6:0.1:1.0, 5); n0007
% exhaustiveSearch(0.7, 0.1:0.1:0.5, 5); n0008
% exhaustiveSearch(0.7, 0.6:0.1:1.0, 5); n0008
% exhaustiveSearch(0.8, 0.1:0.1:0.5, 5); n0009
% exhaustiveSearch(0.8, 0.6:0.1:1.0, 5); n0009
% exhaustiveSearch(0.9, 0.1:0.1:0.5, 5); n0011
% exhaustiveSearch(0.8, 0.6:0.1:1.0, 5); n0011
% exhaustiveSearch(1.0, 0.1:0.1:0.5, 5); n0012
% exhaustiveSearch(1.0, 0.6:0.1:1.0, 5); n0012