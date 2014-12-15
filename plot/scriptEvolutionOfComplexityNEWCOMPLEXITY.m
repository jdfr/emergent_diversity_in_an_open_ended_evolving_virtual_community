function scriptEvolutionOfComplexityNEWCOMPLEXITY

basedir = 'gemabestia\G0.7F0.5\2009_Jul_05_17_08_32\1=0.0001_4';

shownProtos = [...
  ... GENERACION COMPLEJIDAD
  53 304; ...
  74 1169; ...
  90 1122; ...
  117 6415; ...
  129 10070; ...
  166 18279; ...
  241 24292; ...
  241 7519; ...
  ...%241 5421; ...
   ];

maxgen = 243;


level = 9;

namefil = [basedir filesep 'NEWCOMPLEXITY.mat'];

if exist(namefil, 'file')
  res = load(namefil);
  comEspeciesAnte = res.res.comEspeciesAnte;
  comEspeciesSig = res.res.comEspeciesSig;
  structcomplex = res.res.complexityProto;
else
  name2 = [basedir filesep 'TEMPCOMPLEXITY.mat'];
  if exist(name2, 'file')
    structcomplex = load(name2);
    structcomplex = structcomplex.structcomplex;
  else
    fprintf('CALCULATING COMPLEXITY...\n');
    structcomplex = complexityProtoNEWCOMPLEXITY(basedir,1,maxgen, level);
    save(name2, 'structcomplex');
  end
  fprintf('RELATING SPECIES...\n');
  [especies comEspeciesAnte comEspeciesSig]= relacionarSpeciesNEWCOMPLEXITY(basedir,maxgen, structcomplex);
  res = struct('complexityProto', {structcomplex}, 'especies', {especies}, 'comEspeciesAnte', {comEspeciesAnte}, 'comEspeciesSig', {comEspeciesSig}); %#ok<NASGU>
  fprintf('SAVING RESULTS...');
  save(namefil, 'res');
  if exist(name2, 'file')
    delete(name2);
  end
end
% state = rand('twister');
% rand('twister', 0);
% h = showEvolutionOfComplexityNEWCOMPLEXITY(basedir, maxgen, comEspeciesAnte, comEspeciesSig, structcomplex);
% a = get(h, 'CurrentAxes');
% set(a, 'YLim', [0 25000]);%, 'YScale', 'log');
% rand('twister', state);
% return

comEspeciesAnte243 = load([basedir filesep 'comEspeciesAnte243.mat']);
comEspeciesAnte243 = comEspeciesAnte243.comEspeciesAnte243;
comEspeciesSig243  = load([basedir filesep 'comEspeciesSig243.mat']);
comEspeciesSig243  = comEspeciesSig243.comEspeciesSig243;

generateIMAGE(basedir, maxgen, comEspeciesAnte, comEspeciesSig, structcomplex, shownProtos);

function generateIMAGE(basedir, maxgen, comEspeciesAnte, comEspeciesSig, structcomplex, shownProtos)

gens = shownProtos(:,1)+1;
comps = shownProtos(:,2);
idxs = zeros(size(gens));
for k=1:numel(comps)
  idxs(k) = find(structcomplex.complexities{gens(k)}==comps(k));
end

fprintf('LOADING PROTOS\n');
tic;
protos = cell(size(gens));
offs = cell(size(gens));
gs = unique(gens);
for k=1:numel(gs)
  g = gs(k);
  nm = sprintf('numGroups%03d', g);
  numGroups = load([basedir filesep nm]);
  numGroups = numGroups.(nm);
  nm = sprintf('proto%03d', gens(k));
  proto = load([basedir filesep nm]);
  proto = proto.(nm);
  proto = proto{numGroups};
  nm = sprintf('offsetproto%03d', g);
  offsetproto = load([basedir filesep nm]);
  offsetproto = offsetproto.(nm);
  offsetproto = offsetproto{numGroups};
  protos(gens==g) = proto(idxs(gens==g));
  offs(gens==g)   = offsetproto(idxs(gens==g));
end
dims = cellfun(@(x)double([max(x{1}) max(x{2})]), protos, 'uniformoutput', false);
fprintf('LOADED PROTOS %s\n', mat2str(toc));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

fprintf('GENERATING IMAGES...\n');
imgs = cell(numel(protos),1);
dms = dims;
off =offs;
rs = protos;
for zn=1:numel(protos)
  img = ones(dms{zn});
  img(sub2ind(dms{zn}, rs{zn}{:})) = 0;
  img = img(end:-1:1,:);
  imgs{zn}=img;
  clear img;
end


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
  idxToShow = array2cell(letras(1:numel(protos)));
end
paint = true;
showText = false;not(paint);
imscale = 0.0002;
xe = 260;
xlab = 'generations';
ye = 26000;
xtcks = 0:50:250;
ylab = 'phenotype complexity';
ytcks = (0:0.2:2.5)*1e4;
offy=300;
offi=[...
  01 0000 0500 double('a'); ...
  02 0000 0500 double('b'); ...
  03 0000 0500 double('c'); ...
  04 0000 0000 double('d'); ...
  05    0 0000 double('e'); ...
  06    0 0000 double('f'); ...
  07 0-10 -4500 double('g'); ...
  08 0010 0000 double('h'); ...
  09 0000 0000 double('i'); ...
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
  ];

state = rand('twister');
rand('twister', 0);
h = showEvolutionOfComplexityNEWCOMPLEXITY(basedir, maxgen, comEspeciesAnte, comEspeciesSig, structcomplex);
a = get(h, 'CurrentAxes');
set(a, 'YLim', [0 25000]);%, 'YScale', 'log');
rand('twister', state);


  img1 = makeFigure(h, gens-1,  comps, imgs, off, dms, xe, xtcks, xlab, ye, ytcks, ylab, offy, offi, offt, imscale, paint, showText, idxToShow);



if paint
  clear img imgs;
    recortaA1 = 600;
    recortaB1 = 600;
    img1 = img1(:, recortaA1:end-recortaB1,:);
  f = 'lsystemdani\figurasrene\evocompnew.png';
  imwrite(img1, f, 'png');
end
  

function img = makeFigure(h, posx, posy, imgs, off, dms, xe, xtcks, xlab, ye, ytcks, ylab, offy, offs, offt, imscale, paint, showText, idxToShow)

axesargs = {'FontName','times','FontSize',14, 'XLim', [0 xe], 'YLim', [0 ye], 'XTick', xtcks, 'YTick', ytcks};

a = get(h, 'CurrentAxes');

if paint
  v = 'off';
else
  v = 'on';
end

set(h, 'Visible', v, 'Colormap', [0 0 0; 1 1 1; jet(62)]);%, 'PaperUnits', 'points', 'PaperSize', 4*[2000 2000], 'Position', 4*[0 0 2000 2000]);
set(a, axesargs{:});

hold on;

for k=1:numel(imgs)
  lx = (dms{k}(2))*imscale*xe;
  offx = off{k}(2)*imscale*xe;
  ly = dms{k}(1)*imscale*ye;
  xd = offs(k,2)+posx(k)-offx+[0,lx];
  yd = offs(k,3)+posy(k)+ly+offy+[0,-ly];
  img = double(imgs{k});
  szim = size(img);
  img2 = zeros(szim([1 2]));
  img2(img(:,:,1)<1)=1;
  image('CData', img, 'XData', xd, 'YData', yd, 'AlphaDataMapping', 'none', 'AlphaData', img2);
end

line(posx,  posy, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15, 'MarkerEdgeColor', 'k');

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


