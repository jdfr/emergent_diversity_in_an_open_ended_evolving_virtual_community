function makeComparisonGenomes(mild, harsh, veryharsh, filename)
%mild = analyzeGenomeReduced(mildpoph, [], false);
%etc. Esta hecho en el cluster para 3:
%harsh = load('lsystemdani\figurasrene\reducedGenomeStudyG0.75F1_uno_1.mat')
%harsh = harsh.res;
%veryharsh = load('lsystemdani\figurasrene\reducedGenomeStudyG1F1_uno_1.mat')
%veryharsh = veryharsh.res;
%mild = load('lsystemdani\figurasrene\reducedGenomeStudyG0.5F1_dos_1.mat')
%mild = mild.res;

res   = { mild,   harsh,   veryharsh}';
names = {'mild', 'harsh', 'very harsh'}';
nm1 = cellfun(@(x)[x ' (nonreduced)'],         names, 'uniformoutput', false);
nm2 = cellfun(@(x)[x ' (reduced)'], names, 'uniformoutput', false);
nm = {};
legendargs = { 'FontSize', 14};
axesargs = {'FontName','times','FontSize',18};
lw = 2;
lw2 = 4;
linefun = @stairs;%@line; %@stairs;

if not(exist('filename', 'var'))
  filename = '';
end

figargs = {'Color', 'w'};

if not(isempty(filename))
  figargs = [figargs {'visible', 'off'}];
end

h = figure(figargs{:});
  
x=axes;
set(x,axesargs{:}, 'YLim', [0 300]);
hold on
cols = {repmat(0, 1, 3), repmat(0.4, 1, 3), repmat(0.7, 1,3)};
for k=1:numel(res)
  linefun(res{k}.generations, res{k}.nominimized.mean,  'Color', cols{k}, 'LineWidth', lw);
  nm{end+1,1} = nm1{k};
end
for k=1:numel(res)
  linefun(res{k}.generations, res{k}.minimized.mean,  'Color', cols{k}, 'LineWidth', lw2, 'LineStyle', '-');%':');
  nm{end+1,1} = nm2{k};
end
lh = legend(nm,'Location','NorthWest');
set(lh, legendargs{:});
xlabel('generations')
ylabel('genome length')

if not(isempty(filename))
  img = imcapture(h, 'all', 600);
  imwrite(img, filename);
  close(h);
  clear img;
end
  
nm1 = cellfun(@(x)['mean ratio between genome lengths (' x ')'],         names, 'uniformoutput', false);
nm2 = cellfun(@(x)[x ' (reduced)'], names, 'uniformoutput', false);
nm3 = cellfun(@(x)[x ' (nonreduced)'], names, 'uniformoutput', false);
nm =  {};
figure('Color', 'w');
x=axes;
set(x,axesargs{:});
hold on
% for k=1:numel(res)
%   linefun(res{k}.generations, res{k}.compression.mean,  'Color', cols{k}, 'LineWidth', lw);
%   nm{end+1,1} = nm1{k};
% end
for k=1:numel(res)
  linefun(res{k}.generations, res{k}.nonminT3,  'Color', cols{k}, 'LineWidth', lw, 'LineStyle', '-');
  nm{end+1,1} = nm3{k};
end
for k=1:numel(res)
  linefun(res{k}.generations, res{k}.minT3,  'Color', cols{k}, 'LineWidth', lw2, 'LineStyle', '-');%':');
  nm{end+1,1} = nm2{k};
end
lh = legend(nm,'Location','NorthWest');
set(lh, legendargs{:});
xlabel('generations')
ylabel('ratio')
  
