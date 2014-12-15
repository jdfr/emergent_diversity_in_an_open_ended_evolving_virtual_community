function scripShowCompressionVsGnmLen(pophm, compressionm, pophh, compressionh, pophvh, compressionvh, tipog, fname)

% colmap = gray(256);
% colmap=colmap(end:-1:1,:);
colmap = jet(255);
colmap = [1 1 1; colmap];

tittemp = 'Evolution of phenotypic complexity in a %s simulation';

pophs = {pophm pophh pophvh};
comps = {compressionm.compression compressionh.compression compressionvh.compression};

names = {'mild', 'harsh', 'very harsh'};
files = {'M', 'H', 'VH'};
topng =  exist('fname', 'var') && not(isempty(fname));
if topng
  filtemp = [fname '%s.png'];
end

useLogY=true;

for k=1:3
  switch tipog
    case 1
      gnm = pophs{k}.genome;
      xlab = 'genome length';
      gnm = cellfun('prodofsize', gnm);
      xdivs = min(gnm):10:max(gnm);
    case 2
      gnm = pophs{k}.minGenome;
      xlab = 'reduced genome length';
      gnm = cellfun('prodofsize', gnm);
      xdivs = min(gnm):1:max(gnm);
  end
  tit = sprintf(tittemp, names{k});
  h = showLogHist3(gnm, xlab, xdivs, comps{k}, 'complexity', 100, tit, colmap, useLogY);
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
set(get(a, 'YLabel'), labargs{:}, 'string', 'complexity (log scale)');
set(get(a, 'Title'),  labargs{:});

