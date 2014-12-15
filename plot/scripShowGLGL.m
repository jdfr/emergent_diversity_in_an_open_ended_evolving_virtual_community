function scripShowGLGL(pophm, pophh, pophvh, fname)

% colmap = gray(256);
% colmap=colmap(end:-1:1,:);
colmap = jet(255);
colmap = [1 1 1; colmap];

datos = {pophm pophh pophvh };
names = {'mild', 'harsh', 'very harsh'};
files = {'M', 'H', 'VH'};
topng =  exist('fname', 'var') && not(isempty(fname));
if topng
  filtemp = [fname '%s.png'];
end


for k=1:3
  poph = datos{k};
  g1 = cellfun('prodofsize', poph.genome);
  g2 = cellfun('prodofsize', poph.minGenome);
  ming1 = min(g1);
  maxg1 = max(g1);
  ming2 = min(g2);
  maxg2 = max(g2);
  gs1 = ming1:maxg1;
  gs2 = ming2:maxg2;
  h=myhist3(g1, 'genomic length', gs1, [], g2, 'reduced genomic length', gs2, [], true, colmap);
  tit = sprintf('relationship between genome lengths in a %s simulation', names{k});
  title(tit);
  postpros(h);
  if topng
    img = imcapture(h, 'all', 600);
    imwrite(img, sprintf(filtemp, files{k}), 'png');
    close(h);
  end
end


function postpros(h)
figargs =  {'Color', 'w'};
axesargs = {'FontName','times','FontSize',14};
labargs = {'FontName','times','FontSize',14};
a = get(h, 'currentaxes');
set(h, figargs{:});
set(a, axesargs{:});
set(get(a, 'XLabel'), labargs{:});
set(get(a, 'YLabel'), labargs{:});%, 'string', 'complexity (log scale)');
set(get(a, 'Title'),  labargs{:});

