function getLongPopSizes(modo, varargin)

if not(exist('modo', 'var'))
  modo = 'pop';
end

all075 = getBatchs(varargin{:});

alphas = [all075.alphaG]';
gens   = {all075.gens}';
bios   = {all075.biomass}';
pops   = {all075.popsize}';
genmn  = {all075.gnmlenMEAN}';
genrmn = {all075.gnmlenredMEAN}';
disps  = {all075.disps}';
mxg    = max(vertcat(gens{:}));
mxp    = max(vertcat(pops{:}));
mxb    = max(vertcat(bios{:}));
mxg    = max(vertcat(genmn{:}));
mxr    = max(vertcat(genrmn{:}));

%cols = distinguishable_colors(numel(gens), 'w');

%imgs = cell(numel(gens),1);

switch modo
  case 'genom'
    x = gens; xl = 'generations'; y = genmn; yl = 'mean genome size'; ax = [0 mxg 0 mxg];
    axesprops = {};
  case 'genomr'
    x = gens; xl = 'generations'; y = genrmn; yl = 'mean reduced genome size'; ax = [0 mxg 0 mxr];
    axesprops = {};
  case 'pop'
    x = gens; xl = 'generations'; y = pops; yl = 'population size'; ax = [0 mxg 0 mxp];
    axesprops = {};
  case 'bio'
    x = gens; xl = 'generations'; y = bios; yl = 'biomass'; ax = [0 mxg 0 mxb];
    axesprops = {};
  case 'biol'
    x = gens; xl = 'generations'; y = bios; yl = 'biomass'; ax = [0 mxg 0 mxb];
    axesprops = {'YScale', 'log'};
end    

fs = 24;
h = figure;
sz = 20*[2 1];
for k=1:numel(gens);
  cla;
  line(x{k}, y{k});
  grid on;
  xlabel(xl, 'FontSize', fs);
  ylabel(yl, 'FontSize', fs);
  title(sprintf('%04d generations (%02d, \\alpha=%s)', gens{k}(end), k, mat2str(alphas(k))), ...
       'FontSize', fs, 'units', 'normalized', 'Position', [0.2, 0.9]);
  axis(ax);
  if not(isempty(axesprops))
    set(gca, axesprops{:});
  end
  set(gca, 'FontSize', fs);
  set(h, 'PaperUnits', 'centimeters', 'PaperSize', sz, 'PaperPosition', [0 0 sz]);
  saveas(h, ['lsys\longs\' sprintf('alllong%03d.pdf', k)]);
  %close(h); return;
end;
close(h);

%img = vertcat(imgs{:});
