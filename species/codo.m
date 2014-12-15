function codo(imgname, sim, max2show, ejes, threshold)
if not(exist('threshold', 'var'))
  threshold = 1;
end
if ischar(sim)
  S = load(name);%sim290
  fnames = fieldnames(S);

  a = S.(fnames{1});
else
  a = sim;
end
b = abs(diff(a));
bb=diff(b);

c = 1;
if not(exist('max2show', 'var'))
  max2show=25;
end
d = min(max2show, numel(a));

h = figure;
x=axes;
fons = 14;
set(x,'FontName','times','FontSize',fons)
hold on
h2 = stairs((c+1:d-1),(b(c:d-2)), 'color', [0.5 0.5 0.5], 'LineWidth', 2);
h1 = stairs(c:d-1,a(c:d-1),        'color', [0.0 0.0 0.0], 'LineWidth', 2);
%h3 = stairs(c+1:d-2,bb(c+1:d-2),        'color', [1.0 0.0 0.0], 'LineWidth', 2);
% bar(c:d-1,a(c:d-1),'w')
% bar(c+1:d-1,abs(b(c:d-2)),'FaceColor',[0.7 0.7 0.7])
% % bar(c+1:d-1,a(c+1:d-1),'w')
% % bar(c+1:d-1,abs(b(c:d-2)),'FaceColor',[0.7 0.7 0.7])
xlabel('number of clusters','FontSize',fons)
ylabel('dispersion', 'FontSize', fons);%'$R_{m}$,  $\triangle R_{m}$', 'Interpreter', 'latex');
%bar(2:length(sim290),abs(diff(sim290)),'k','facecolor','w')
if  exist('ejes', 'var') && not(isempty(ejes))
  axis(ejes)
end
line([0 40],[threshold threshold],'color','k','linestyle','--')
%ylabel('mean intra-class dispersion','FontName','times','FontSize',12)
legend([h1 h2], {sprintf('mean dispersion $R_{m}$'),'decrease in mean dispersion $\triangle R_{m}$'}, 'Interpreter', 'latex', 'FontSize', fons)
legend(x, 'boxoff');

if not(isempty(imgname))
  img = imcapture(h, 'all', 600);
  imwrite(img, imgname, 'png');
  close(h);
end
