function scripShowCompression(pophm, compressionm, pophh, compressionh, pophvh, compressionvh, fname)

% colmap = gray(256);
% colmap=colmap(end:-1:1,:);
colmap = jet(255);
colmap = [1 1 1; colmap];

tittemp = 'Evolution of phenotypic complexity in a %s simulation';


datos = {{pophm compressionm.compression} {pophh compressionh.compression} {pophvh compressionvh.compression}};
names = {'mild', 'harsh', 'very harsh'};
files = {'M', 'H', 'VH'};
topng =  exist('fname', 'var') && not(isempty(fname));
if topng
  filtemp = [fname '%s.png'];
end


for k=1:3
  tit = sprintf(tittemp, names{k});
  h = showCompression(datos{k}{:}, tit, colmap, true);
  postpros(h);
  if topng
    img = imcapture(h, 'all', 600);
    imwrite(img, sprintf(filtemp, files{k}), 'png');
    close(h);
  end
end

% tit = sprintf(tittemp, 'mild');
% h = showCompression(pophm, compressionm.compression, tit, colmap, true);
% postpros(h);
% 
% tit = sprintf(tittemp, 'harsh');
% h = showCompression(pophh, compressionh.compression, tit, colmap, true);
% postpros(h);
% 
% tit = sprintf(tittemp, 'veryharsh');
% h = showCompression(pophvh, compressionvh.compression, tit, colmap, true);
% postpros(h);
% 

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

