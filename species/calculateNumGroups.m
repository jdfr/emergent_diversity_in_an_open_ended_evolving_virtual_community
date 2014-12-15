function numGroups = calculateNumGroups(sim, options)

if sim == 0
    derivada = 0;
else
  switch options.thresholdMode(1:6)
    case {'absolu', 'abslas', 'reldif'}
      derivada = abs(diff(sim));
      if options.thresholdMode(1)=='r'
        derivada = derivada./max(derivada);
      end
    case 'relval'
      ss = sim./max(sim);
      derivada = abs(diff(ss));
  end
end

% h1 = figure('visible','off');
% bar(sim);
% xlabel('number of groups');
% ylabel('dispersion');
% saveas(h1,[basedir filesep 'dispersion' num2str(Generacion, '%03g')],'png');
% close(h1);
% 
% h2 = figure('visible','off');
% bar(derivada);
% xlabel('number of groups');
% ylabel('derivative of the dispersion');
% saveas(h2,[basedir filesep 'derivadaDispersion' num2str(Generacion, '%03g')],'png');
% close(h2);

switch options.thresholdMode
  case {'absolute', 'reldiff', 'relval'}
    f = find(derivada <= options.threshold, 1, 'first'); %classical mode
  case {'abslast', 'reldifflast', 'relvallast'}
    below = derivada<=options.threshold;
    if numel(derivada)<2
      f = find(below, 1, 'first');
    else
      %first time such that all subsequent values are below threshold
      if not(below(end))
        f = length(derivada)+1;
      else
        f = find(diff(below)==1, 1, 'last')+1;
        if isempty(f)
          f = find(below, 1, 'first');
        end
      end
    end
  otherwise
      error('options options.thresholdMode=<%s> not understood!!!', options.thresholdMode);
end
if numel(f) == 0
    f = length(derivada)+1;
end
numGroups = f(1);

