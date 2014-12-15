function figuras = makeFigurasStruct(option, rellenar, alsoByInds)
if nargin<2
  rellenar = false;
end
if nargin<3
  alsoByInds = false;
end

switch (option)
  case 1
    veryharsh = 'genomas\porhacer\G1F1\2009_Jul_29_14_22_12\5\1=0.001_1';
    harsh     = 'genomas\porhacer\G0.8F1\2009_Jul_29_12_21_07\2\1=0.001_1';
    mild      = 'genomas\porhacer\G0.5F1\2009_Jul_29_12_16_39\3\1=0.001_1';
    basedirs = {veryharsh, harsh, mild}';
    tipos = {'very harsh', 'harsh', 'mild'};
    colors = {'r', 'b', 'g'};
    indexes = {1, 2, 3};
    prefix = '';
  case 2
    veryharsh = 'newdists\editdist\2009_Jul_29_14_22_12\1=0.001_1';
    harsh     = 'newdists\editdist\2009_Jul_29_12_21_07\1=0.001_1';
    mild      = 'newdists\editdist\2009_Jul_29_12_16_39\1=0.001_1';
    basedirs = {veryharsh, harsh, mild}';
    tipos = {'very harsh', 'harsh', 'mild'};
    colors = {'r', 'b', 'g'};
    indexes = {1, 2, 3};
    prefix = '';
  case 3
    veryharsh = 'newdists\mutdist\2009_Jul_29_14_22_12\1=0.001_1';
    harsh     = 'newdists\mutdist\2009_Jul_29_12_21_07\1=0.001_1';
    mild      = 'newdists\mutdist\2009_Jul_29_12_16_39\1=0.001_1';
    basedirs = {veryharsh, harsh, mild}';
    tipos = {'very harsh', 'harsh', 'mild'};
    colors = {'r', 'b', 'g'};
    indexes = {1, 2, 3};
    prefix = '';
  case 4
    veryharsh = 'newdists\miniedit\2009_Jul_29_14_22_12\1=0.001_1';
    harsh     = 'newdists\miniedit\2009_Jul_29_12_21_07\1=0.001_1';
    mild      = 'newdists\miniedit\2009_Jul_29_12_16_39\1=0.001_1';
    basedirs = {veryharsh, harsh, mild}';
    tipos = {'very harsh', 'harsh', 'mild'};
    colors = {'r', 'b', 'g'};
    indexes = {1, 2, 3};
    prefix = '';
  case 5
%     mild = '\\150.214.56.111\l-system\l-system\codigo\codigo\alphaGalphaF\G0.1F0.5\2009_Jul_04_11_51_30\1=0.0001_3';
%    veryharsh     = '\\150.214.56.111\l-system\l-system\codigo\codigo\alphaGalphaF\G1F0.9\2009_Jul_04_13_25_15\1=0.0001_11';
%     harsh      = '\\150.214.56.111\l-system\l-system\codigo\codigo\alphaGalphaF\G0.7F0.2\2009_Jul_04_03_31_39\1=0.0001_1';
   mild = 'fromlsystem\g0.1f0.5\1';
   harsh = 'fromlsystem\g0.7F0.2\1=0.0001_1';
   veryharsh= 'fromlsystem\g1f0.9\1';
    basedirs = {veryharsh, harsh, mild}';
    tipos = {'very harsh', 'harsh', 'mild'};
    colors = {'r', 'b', 'g'};
    indexes = {1, 2, 3};
    prefix = '';
  otherwise
    error('unknown option!!!');
end
figuras = makeFigurasStructDOIT(rellenar, alsoByInds, basedirs, prefix, tipos, colors, indexes);

function figuras = makeFigurasStructDOIT(rellenar, alsoByInds, basedirs, prefix, tipos, colors, indexes)
% figuras = makeFigurasStruct({'newdists\editdist\2009_Jul_29_12_16_39\1=0.001_1', 'newdists\editdist\2009_Jul_29_14_22_12\1=0.001_1', );
% todoEnoUnEDITDISTANCE();
% todoEnoUnEDITDISTANCE('newdists\editdist\2009_Jul_29_12_21_07\1=0.001_1');

campos = {'nombres'; ...
          'generations'; ...
          'popsizes'; ...
          'numRamas'; ...
          'numSpecies'; ...
          'speciesGenerations'; ...
          'legend'; ...
          'color'; ...
          'order'};
        
if nargin<3
  colors = mat2cell(lines(numel(basedirs)), ones(numel(basedirs),1), 3);
end
if nargin<4
  indexes = array2cell(1:numel(basedirs));
end
if not(iscell(indexes))
  indexes = array2cell(indexes);
end
tabla = cell(numel(campos), numel(basedirs));

tabla(1,:) = basedirs(:)';
tabla(7,:) = tipos(:)';
tabla(8,:) = colors(:)';
tabla(9,:) = indexes(:)';


for k=1:numel(basedirs)
  b = basedirs{k};
  fprintf('Going for %s\n', b);
  fn = [b filesep 'poph.mat'];
  if exist(fn, 'file')
    poph = load(fn);
    poph = poph.poph;
  else
    poph = loadSimulation(b);
  end
  mng = min(poph.generation);
  mxg = max(poph.generation);
  gs = (mng:mxg)';
  tabla{2,k} = gs;
  pops = zeros(size(gs));
  bios = zeros(size(gs));
  fprintf('   pops+bios...\n');
  for z=1:numel(gs)
    thisg = find(poph.generation==gs(z));
    thisg = thisg(poph.nPixelsInSoil(thisg)==0);
    pops(z) = numel(thisg);
    bios(z) = sum(poph.nbranches(thisg));
  end
  tabla{3,k} = pops;
  tabla{4,k} = bios;
    gss = (1:mxg)';
    numSpecies = zeros(size(gss));
    if alsoByInds
      maxSpecies = zeros(size(gss));
    end
    ok = true;
    fprintf('   numSpecies...\n');
    gemafile = [b filesep prefix 'numGroups.mat'];
    if exist(gemafile, 'file')
      especies = load(gemafile);
      especies = especies.especies;
      tabla{5,k} = [especies.numGroups{:}]';
      tabla{6,k} = especies.generations(:);
      tabla{10,k} = cellfun(@numel, especies.individualsSpecies);
      tabla{10,k} = tabla{10,k}(:);
    else
      for z=1:mxg
        nm = sprintf('%snumGroups%03g', prefix, gss(z));
        name = [b filesep nm '.mat'];
        if alsoByInds
          nm2 = sprintf('%sindividualsSpecies%03g', prefix, gss(z));
          name2 = [b filesep nm2 '.mat'];
        end
        ok = exist(name, 'file');
        if ok
          s = load(name);
          numSpecies(z) = s.(nm);
          if alsoByInds
            s2 = load(name2);
            maxSpecies(z) = numel(s2.(nm2));
          end
        else
          if rellenar
            numSpecies(z) = numSpecies(z-1);
            if alsoByInds
              maxSpecies(z) = maxSpecies(z-1);
            end
          else
            break;
          end
        end
      end
      if not(ok)
        gss = gss(1:z-1);
        numSpecies = numSpecies(1:z-1);
      end
      tabla{5,k} = numSpecies;
      tabla{6,k} = gss;
      if alsoByInds
        tabla{10,k} = maxSpecies;
      end
    end
end

figuras = struct('campos', {campos}, ...
                 'tipos', {tipos}, ...
                 'tabla', {tabla});
