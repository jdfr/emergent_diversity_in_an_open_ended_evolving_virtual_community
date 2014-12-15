function makePhaseDiagramPopBio2(allres, mode, xax, yax, gen, alphas, axprops, filename)

%allong = getBatchs([], {'110NB'}); makePhaseDiagramPopBio(allong, 'diff', 'pop', 'bio', [0:1000], 0.5:0.05:1, {'XLim', [0 2600], 'YLim', [1 275e6], 'YTick', 10.^(1:8)}, 'phaseplot%04d.pdf');    

cf = @(varargin)cellfun(varargin{:}, 'uniformoutput', false);

af = @(varargin)arrayfun(varargin{:}, 'uniformoutput', false);

dirs = strrep({allres.dir}', '\', '/');
slashes = cf(@(x)find(x=='/'), dirs);
dirs = cf(@(x,y)x(y(end-3)+1:y(end)-1), dirs, slashes);

gs = [allres.alphaG]';

gens = {allres.gens}';
numgs = cellfun(@max, gens);

%fitm = {allres.fitnessMEAN}';

ngens = cellfun(@numel, gens);

if not(exist('filename', 'var'))
  filename = '';
end

interactive = isempty(filename);
if interactive
  figargs = {};
else
  figargs = {'Visible', 'off'};
end

if not(exist('alphas', 'var')) || isempty(alphas)
  alphas = unique(gs);
end

%mymap = distinguishable_colors(numel(alphas), 'w'); mymap2 = mymap; m = 0.025;
mymap = jet(numel(alphas));                         mymap2 = mymap; m = 0.025;
%mymap = jet(numel(alphas));                         mymap2 = jet(256); m = 0;

xax = getax(xax, allres);
yax = getax(yax, allres);

fs1 = {'FontName', 'Times', 'FontSize', 16};
fs2 = {'FontName', 'Times', 'FontSize', 12};

samefig = iscell(gen);
if samefig
  gen = gen{1};
  h = figure(figargs{:});
end

lineargs1 = {'Marker', 'o', 'MarkerSize', 3};
lineargs2 = {'Marker', 's', 'MarkerSize', 4};
tb = 0;
for z=0:numel(gen)
  if tb>0; fprintf(repmat('\b', 1, tb)); end
  if z>0
    tb = fprintf('gen %04d/%04d', gen(z), gen(end));
  end

  if samefig
    clf(h);
  else
    h = figure(figargs{:});
  end

  if z>0
    g = gen(z);
  else
    g=nan;
  end
  
  tx         = xax.d;
  ty         = yax.d;
  tdirs      = dirs;
  tgens      = gens;
  tnumgs     = numgs;
  tgs        = gs;
  
  if z>0
    switch mode
      case 'all'
        innerloop(tgs, mymap, alphas, g, false, lineargs1, 0, tgens, tx, ty, tdirs, interactive);
      case {'exclude', 'diff'}
        todelete = tnumgs<g;
        switch mode
          case 'diff'
            tx         = xax.d(todelete);
            ty         = yax.d(todelete);
            tdirs      = dirs(todelete);
            tgens      = gens(todelete);
            tnumgs     = numgs(todelete);
            tgs        = gs(todelete);
            innerloop(tgs, mymap, alphas, g, true, lineargs2, -1, tgens, tx, ty, tdirs, interactive);
        end
        toremain   = not(todelete);
        tx         = xax.d(toremain);
        ty         = yax.d(toremain);
        tdirs      = dirs(toremain);
        tgens      = gens(toremain);
        tnumgs     = numgs(toremain);
        tgs        = gs(toremain);
        innerloop(tgs, mymap, alphas, g, true, lineargs1, 0, tgens, tx, ty, tdirs, interactive);
      otherwise
        error('unknown mode <%s>!!!', mode);
    end
  end

  a = get(h, 'CurrentAxes');
  set(a, fs1{:});%, 'YScale', 'log', 'YLim', [2000 2.5e8], 'YTick', [1e3 1e4 1e5 1e6 1e7 1e8]);%, 'XLim', [1 1500], 'XTick', [1 1e2 1e3 1e4]);
  if exist('axprops', 'var')
    set(a, axprops{:});
  end
  %grid on;

  if false; interactive;
    dcm_obj = datacursormode(h);
    set(dcm_obj, 'Updatefcn', @showdir)
  end

  colormap(mymap2);
  caxis([min(alphas)-m max(alphas)+m]);
  xlabel(xax.l);
  ylabel(yax.l);
  %zlabel('generations');
%   title(sprintf('state of simulations at generation %03d', g));
  %grid on;
  %view(3);
  cbar_axes = colorbar('YTick', 0.5:0.05:1, fs1{:}); 
  ylabel(cbar_axes, '\alpha value', fs1{:});
  
  if z>0
    set(a,         'Visible', 'off');
    set(cbar_axes, 'Visible', 'off');
  end
  
  if not(isempty(filename))
    saveas(h, sprintf(filename, z));
    %export_fig(sprintf(filename, z), '-pdf', '-nocrop', h);
    if not(samefig)
      close(h);
    end
  else
    drawnow;
  end
end

if not(isempty(filename)) && samefig
  close(h);
end

if tb>0; fprintf(repmat('\b', 1, tb)); end
%fprintf('\n');

function txt = showdir(obj, ev)
target = get(ev, 'Target');
if strcmp(get(target, 'type'), 'line')
  txt = get(target, 'UserData');
else
  txt = 'not available';
end
if isempty(txt)
  txt = 'empty';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function innerloop(tgs, mymap, alphas, g, exclude, lineargs, colfac, tgens, tx, ty, tdirs, interactive)
  for k=1:numel(tgs)
    c = mymap(tgs(k)==alphas,:);
    posg = g+1;
    if interactive
      mor = {'UserData', tdirs{k}};
    else
      mor = {};
    end
    %if not(exclude)
      posg = min(numel(tgens{k}), posg);
    %end
  %   bs = log(bios{k}); bs(isinf(bs)) = 0;
  %   modo = 1:5:numel(gens{k});
  %   line(pops{k}(modo), bs(modo), gens{k}(modo), 'Color', c, ... %'Marker', 'o', 'MarkerEdgeColor', c, 'MarkerFaceColor', c, 'MarkerSize', 5, ...
  %                                'UserData', dirs{k});
    line(tx{k}(posg), ty{k}(posg), tgens{k}(posg), 'MarkerEdgeColor', max(0, min(1, c+colfac)), 'MarkerFaceColor', c, lineargs{:}, mor{:});

  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function st = getax(ax, allres)
switch ax
  case 'mpress'
    st.l = 'mean pressure'; st.d = {allres.pressureMEAN}';
  case 'dfit'
    st.l = 'std fitness'; st.d = {allres.fitnessMEAN}';
  case 'fit'
    st.l = 'mean fitness'; st.d = {allres.fitnessMEAN}';
  case 'fitpot'
    st.l = 'mean potential fitness'; st.d = {allres.fitnesspotMEAN}';
  case 'dbio'
    st.l = 'std of tree biomass'; st.d={allres.nbranchesSTD}';
  case 'mbio'
    st.l = 'mean biomass'; st.d={allres.nbranchesMEAN}';
  case 'bio'
    st.l = 'biomass'; st.d = {allres.biomass}';
  case 'mh'
    st.l = 'mean height'; st.d={allres.heightMEAN}';
  case 'pop'
    st.l = 'population size'; st.d = {allres.popsize}';
  case 'genom'
    st.l = 'mean genome length'; st.d = {allres.gnmlenMEAN}';
  case 'genomr'
    st.l = 'mean reduced genome length'; st.d = {allres.gnmlenredMEAN}';
  case 'disp.pht'
    st.l = 'mean dispersion (phenotype distance)'; st.d = arrayfunc(@(x)x.disps.dispersion.pht.mean, allres);
  case 'disp.gnm'
    st.l = 'mean dispersion (genome distance)'; st.d = arrayfunc(@(x)x.disps.dispersion.gnm.mean, allres);
  case 'disp.mng'
    st.l = 'mean dispersion (reduced genome distance)'; st.d = arrayfunc(@(x)x.disps.dispersion.mng.mean, allres);
end  


function c = parsecolor(s)
  if ischar(s)
    c = colorstr2rgb(s);
  elseif isnumeric(s) && size(s,2) == 3
    c = s;
  else
    error('MATLAB:InvalidColorSpec','Color specification cannot be parsed.');
  end


function c = colorstr2rgb(c)
  % Convert a color string to an RGB value.
  % This is cribbed from Matlab's whitebg function.
  % Why don't they make this a stand-alone function?
  rgbspec = [1 0 0;0 1 0;0 0 1;1 1 1;0 1 1;1 0 1;1 1 0;0 0 0];
  cspec = 'rgbwcmyk';
  k = find(cspec==c(1));
  if isempty(k)
    error('MATLAB:InvalidColorString','Unknown color string.');
  end
  if k~=3 || length(c)==1,
    c = rgbspec(k,:);
  elseif length(c)>2,
    if strcmpi(c(1:3),'bla')
      c = [0 0 0];
    elseif strcmpi(c(1:3),'blu')
      c = [0 0 1];
    else
      error('MATLAB:UnknownColorString', 'Unknown color string.');
    end
  end

