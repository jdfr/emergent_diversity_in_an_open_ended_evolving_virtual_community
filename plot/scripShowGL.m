function scripShowGL(pophm, pophh, pophvh, tipog, fname)

% colmap = gray(256);
% colmap=colmap(end:-1:1,:);
colmap = jet(255);
colmap = [1 1 1; colmap];

switch tipog
  case 1
    tittemp = 'Evolution of genomic length in a %s simulation';
    fld = 'genome';
    nm = 'genome length';
    gap = 1;15;
  case 2
    tittemp = 'Evolution of reduced genomic length in a %s simulation';
    fld = 'minGenome';
    nm = 'reduced genome length';
    gap = 1;4;
end


datos = {pophm pophh pophvh };
names = {'mild', 'harsh', 'very harsh'};
files = {'M', 'H', 'VH'};
topng =  exist('fname', 'var') && not(isempty(fname));
if topng
  filtemp = [fname '%s.png'];
end


for k=1:3
  poph = datos{k};
  ming = min(poph.generation);
  maxg = max(poph.generation);
  gs = ming:maxg;
  generations = poph.generation;
  gnm = cellfun('prodofsize',poph.(fld));
  minc = min(gnm);
  maxc = max(gnm);
  h=myhist3(generations, 'generation', gs, [], gnm, nm, minc:gap:maxc, [], true, colmap);
  tit = sprintf(tittemp, names{k});
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

