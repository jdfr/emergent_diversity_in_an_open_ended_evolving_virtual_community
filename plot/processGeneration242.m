function processGeneration242(basedir, P, minGenome, proto, numGroups, comLZfigureLongitud, lgenome, lmingenome)
if not(exist('basedir', 'var')) || isempty(basedir)
  basedir = 'lsystemdani\figurasrene\genomecomp\a\1';
end
if not(exist('P', 'var')) || isempty(P)
  P = load([basedir filesep 'P242.mat']);
  P = P.P242;
end
if not(exist('minGenome', 'var')) || isempty(minGenome)
  minGenome = cellfun(@sanitizeString, P.genome, 'uniformoutput', false);
end
if not(exist('proto', 'var')) || isempty(proto)
  proto = load([basedir filesep 'proto242.mat']);
  proto = proto.proto242;
end
if not(exist('numGroups', 'var')) || isempty(numGroups)
  numGroups = load([basedir filesep 'numGroups242.mat']);
  numGroups = numGroups.numGroups242;
end
if not(exist('comLZfigureLongitud', 'var')) || isempty(comLZfigureLongitud)
  comLZfigureLongitud = load([basedir filesep 'comLZfigureLongitud242.mat']);
  comLZfigureLongitud = comLZfigureLongitud.comLZfigureLongitud242;
end
if not(exist('lgenome', 'var')) || isempty(lgenome)
  lgenome = load([basedir filesep 'lgenome242.mat']);
  lgenome = lgenome.lgenome242;
end
if not(exist('lmingenome', 'var')) || isempty(lmingenome)
  lmingenome = load([basedir filesep 'lmingenome242.mat']);
  lmingenome = lmingenome.lmingenome242;
end

proto = proto{numGroups};
if numel(lgenome)~=numel(proto)
  error('jarrlll %d/%d', numel(lgenome), numel(proto));
end
if numel(lmingenome)~=numel(proto)
  error('jarrlll %d/%d', numel(lmingenome), numel(proto));
end

lenGenome  = cellfun('prodofsize', P.genome);
lenMGenome = cellfun('prodofsize', minGenome);

sameRaster     = false(numel(P.raster), numel(proto));
sameLenGenome  = false(numel(P.raster), numel(proto));
sameLenMGenome = false(numel(P.raster), numel(proto));

raster = P.raster;

fprintf('CALCULATING SAMEMATRICES\n');
for a=1:size(sameRaster,1)
  for b=1:size(sameRaster,2)
    sameRaster(a,b) = isequal(raster{a}, proto{b});
    sameLenGenome(a,b) = lenGenome(a)==lgenome(b);
    sameLenMGenome(a,b) = lenMGenome(a)==lmingenome(b);
  end
end

sameAll = sameRaster & sameLenGenome;% & sameLenMGenome;

[idx numproto] =find(sameAll);

[numproto pi]=unique(numproto);
idx = idx(pi);

fprintf('GENERATING IMAGES...\n');
imgs = cell(numel(numproto),1);
dms = P.dimensions(idx,:);
off =P.offset(idx,:);
rs = P.raster(idx);
for zn=1:numel(numproto)
  img = ones(dms(zn,:));
  img(sub2ind(dms(zn,:), rs{zn}{:})) = 0;
  img = img(end:-1:1,:);
  imgs{zn}=img;
  clear img;
end

comps = comLZfigureLongitud(numproto)';
lens  = lenGenome(idx);
mlens = lenMGenome(idx);

%zz=[array2cell([numproto, idx, comLZfigureLongitud(numproto)', lgenome(numproto)', lenGenome(idx), lenMGenome(idx)]), minGenome(idx), P.genome(idx)];

rescale = true;

if rescale
  fprintf('RESCALING...\n');
  firstDownScale = 0.8;
  for k=1:numel(imgs)
    imgs{k} = imresize(imgs{k}, firstDownScale);
    imgs{k}(imgs{k}<0) = 0;
    imgs{k}(imgs{k}>1) = 1;
    imgs{k} = cat(3, imgs{[k k k]});
  end
  fprintf('RESCALED...\n');
end

indexesSortByComplexity = false;
letras = 'a':'z';
if indexesSortByComplexity
  [idxToShow idxToShow] = sort(comps);
  idxToShow = array2cell(letras(idxToShow));
else
  idxToShow = array2cell(letras(1:numel(numproto)));
end
paint = true;
showText = false;not(paint);
show1 = true;
show2 = true;
imscale = 0.0002;
xe = 140;
xlab = 'genotype length';
ye = 27000;
xtcks = 0:20:140;
ylab = 'phenotype complexity';
ytcks = (0:0.2:2)*1e4;
offy=300;
offi=[...
  01 0000 0000 double('a'); ...
  02 0000 0000 double('b'); ...
  03 0000 0000 double('c'); ...
  04 0000 0000 double('d'); ...
  05   -1 0000 double('e'); ...
  06    1 0000 double('f'); ...
  07 0000 0000 double('g'); ...
  08 0000 0000 double('h'); ...
  09 0003 0000 double('i'); ...
  10 0000 0000 double('j'); ...
  11 0003 0000 double('k'); ...
  12    5 -600 double('l'); ...
  13 0000 0000 double('m'); ...
  14   -2 0000 double('n'); ...
  15 0003 0000 double('o'); ...
  16 0000 0000 double('p'); ...
  17 0000 0000 double('q'); ...
  18 0000 0000 double('r'); ...
  19   -4 0000 double('s'); ...
  20 0005 0000 double('t'); ...
  ];
offt=[...
  01 0000 0000 double('a'); ...
  02 0000 0000 double('b'); ...
  03 0000 0000 double('c'); ...
  04 0000 0000 double('d'); ...
  05 0000 0000 double('e'); ...
  06 0000 0000 double('f'); ...
  07 0000 0000 double('g'); ...
  08 0000 0000 double('h'); ...
  09 0000 0000 double('i'); ...
  10 0000 0000 double('j'); ...
  11 0000 0000 double('k'); ...
  12 0000 0000 double('l'); ...
  13 0000 0000 double('m'); ...
  14 0000 0000 double('n'); ...
  15 0000 0000 double('o'); ...
  16 0000 0000 double('p'); ...
  17 0000 0000 double('q'); ...
  18 0000 0000 double('r'); ...
  19 0000 0000 double('s'); ...
  20 0000 0000 double('t'); ...
  ];

if show1
  img1 = makeFigure(lens,  comps, imgs, off, dms, xe, xtcks, xlab, ye, ytcks, ylab, offy, offi, offt, imscale, paint, showText, idxToShow);
end

xe = 60;
xtcks = 0:10:60;
xlab = 'reduced genotype length';
offi=[...
  01 0000 0000 double('a'); ...
  02 0000 0000 double('b'); ...
  03 0000 0000 double('c'); ...
  04 0000 0000 double('d'); ...
  05 01.0 -600 double('e'); ...
  06 01.0 -600 double('f'); ...
  07   -1 0000 double('g'); ...
  08 00.4 0000 double('h'); ...
  09 0000 0000 double('i'); ...
  10 0000 0000 double('j'); ...
  11 0000 0000 double('k'); ...
  12 01.8 -600 double('l'); ...
  13 0000 0000 double('m'); ...
  14 -0.1 0000 double('n'); ...
  15 0000 0000 double('o'); ...
  16 0000 0000 double('p'); ...
  17 0000 0000 double('q'); ...
  18 0000 0000 double('r'); ...
  19 -0.8 -600 double('s'); ...
  20  3.5 -600 double('t'); ...
  ];
offt=[...
  01 0000 0000 double('a'); ...
  02 0000 0000 double('b'); ...
  03 0000 0000 double('c'); ...
  04 0000 0000 double('d'); ...
  05 0000 0000 double('e'); ...
  06 0000 0000 double('f'); ...
  07 0000 0000 double('g'); ...
  08 0000 0000 double('h'); ...
  09 0000 0000 double('i'); ...
  10 0000 0000 double('j'); ...
  11 0000 0000 double('k'); ...
  12 0000 0000 double('l'); ...
  13 0000 0000 double('m'); ...
  14 0000 0000 double('n'); ...
  15 0000 0000 double('o'); ...
  16 0000 0000 double('p'); ...
  17 0000 0000 double('q'); ...
  18 0000 0000 double('r'); ...
  19 0000 0000 double('s'); ...
  20 0000 0000 double('t'); ...
  ];

if show2
  img2 = makeFigure(mlens, comps, imgs, off, dms, xe, xtcks, xlab, ye, ytcks, ylab, offy, offi, offt, imscale, paint, showText, idxToShow);
end

if paint
  clear img imgs;
  if show1 && show2
    recortaA1 = 600;
    recortaB1 = 600;
    recortaA2 = 600;
    recortaB2 = 600;
    img1 = img1(:, recortaA1:end-recortaB1,:);
    img2 = img2(:, recortaA2:end-recortaB2,:);
    img = [img1, img2];
  elseif show1
    img = img1;
  elseif show2
    img = img2;
  else
    img = [];
  end
  clear img1 img2;
  f = 'lsystemdani\figurasrene\prueba.png';
  imwrite(img, f, 'png');
end
  

function img = makeFigure(posx, posy, imgs, off, dms, xe, xtcks, xlab, ye, ytcks, ylab, offy, offs, offt, imscale, paint, showText, idxToShow)

axesargs = {'FontName','times','FontSize',14, 'XLim', [0 xe], 'YLim', [0 ye], 'XTick', xtcks, 'YTick', ytcks};

h = figure;
a = axes;

if paint
  v = 'off';
else
  v = 'on';
end

set(h, 'Visible', v, 'Colormap', [0 0 0; 1 1 1; jet(62)]);%, 'PaperUnits', 'points', 'PaperSize', 4*[2000 2000], 'Position', 4*[0 0 2000 2000]);
set(a, axesargs{:});

hold on;

for k=1:numel(imgs)
  lx = (dms(k,2))*imscale*xe;
  offx = off(k,2)*imscale*xe;
  ly = dms(k,1)*imscale*ye;
  xd = offs(k,2)+posx(k)-offx+[0,lx];
  yd = offs(k,3)+posy(k)+ly+offy+[0,-ly];
  img = double(imgs{k});
  szim = size(img);
  img2 = zeros(szim([1 2]));
  img2(img(:,:,1)<1)=1;
  image('CData', img, 'XData', xd, 'YData', yd, 'AlphaDataMapping', 'none', 'AlphaData', img2);
end

line(posx,  posy, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
% xl = [lens mlens nan(size(lens))]';
% yl = [comps comps nan(size(lens))]';
% line(xl,    yl,    'LineStyle', '-', 'Color', 'k');
%line(lens,  comps, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
%line(mlens, comps, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');

if showText
  for k=1:numel(imgs)
    text(posx(k)+offt(k,2), posy(k)-300+offt(k,3), idxToShow{k}, 'FontName','times','FontSize',14);
  end
end

xlabel(xlab);
ylabel(ylab);

axis square;

%saveas(h, 'lsystemdani\figurasrene\prueba.png', 'png');

if paint
  img = imcapture(h, 'all', 600);
%   imwrite(img, fname, 'png');
  close(h);
else
  img = [];
  grid on;
end
