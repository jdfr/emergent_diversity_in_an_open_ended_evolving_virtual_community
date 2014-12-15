function showSlanted

types = {'normal' 'fewtrees', 'slanted'};
type_testLONG      = [2 1 1 3 1 3 3 1 1 1 1 2 1 1 3 1 1 2 3 1];
type_test110NB     = [1 1 2 1 1 1 1 3 3 3];
type_test110NB_old = [1 1 1 1 1 2 1 2 3 1];
type_testOKPRFX    = [3 1 3 3 3 3 3 3 3 2 3 3 3 2 3 3 1 3 3 1];

names = {'testLONG', 'test110NB', 'test110NB_old', 'testOKPRFX'};
tits = names;

figure;
for k=1:numel(names)
  subplot(2, 2, k);
  v = eval(['type_' names{k}]);
  hist(v, 1:3);
  title(sprintf('%s (%d)', tits{k}, numel(v)));
  set(gca, 'XTick', 1:3, 'XTickLabel', types);
  rotateticklabel(gca,90);
end
  