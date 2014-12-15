function testSpeciesPredecesorNOPROTO(useNew, poph, thefiles)

maxgen = 442;
umbral = 10;
mode = 'morpho';

basedir1 = 'bestia\variable_pop_retune\afinando1\T2_A20\2009_Jun_24_12_31_04\1=0.0001_1';
basedir2 = 'bestia\variable_pop_retune\afinando1\T2_A20\continua_guai_2009_Jun_24_12_31_04\1=0.0001_1';

if not(exist('poph', 'var')) || isempty(poph)
  poph = loadSimulation(basedir1);
end

% nmg = [basedir2 filesep 'minGenome' mat2str(maxgen) '.mat'];
% if exist(nmg, 'file')
%   minGenome = load(nmg);
%   minGenome = minGenome.minGenome;
% else
%   minGenome = cell(find(poph.generation==maxgen, 1, 'last'),1);
%   gnm = poph.genome;
%   for k=1:numel(minGenome)
%     if mod(k,1000)==0
%       fprintf('%d/%d\n', k, numel(minGenome));
%     end
%     minGenome{k} = sanitizeString(gnm{k});
%   end
%   clear gnm;
%   save(nmg, 'minGenome');
% end
% 
% poph.minGenome = minGenome;

if useNew
  if exist('thefiles', 'var')
    matrixDist = thefiles{1};
    indicesClasificados = thefiles{2};
  else
    if not(exist([basedir2 filesep 'matrixDist.mat'], 'file'))
      HierarchicalClusteringALL(basedir2,poph,'morpho',1,maxgen,false,'', '', struct('justMatrixDist', true, 'alsoClasificados', true));
      compressGemaFiles(basedir2, 1, maxgen, '', 5);
    end
    matrixDist = load([basedir2 filesep 'matrixDist.mat']);
    matrixDist = matrixDist.matrixDist;
    indicesClasificados = load([basedir2 filesep 'indicesClasificados.mat']);
    indicesClasificados = indicesClasificados.indicesClasificados;
  end
  ret= speciesPredecesorNOPROTO(basedir2,poph,umbral,maxgen, mode, matrixDist, indicesClasificados);
  save([basedir2 filesep 'miniSpeciationNUEVANUEVANUEVA.mat'], 'ret');
%     ret= speciesPredecesorNOPROTO(basedir2,poph,umbral,maxgen, mode);%, matrixDistsALL, indicesClasificadosALL);
%     save([basedir2 filesep 'miniSpeciationNUEVA.mat'], 'ret');
else
  [numSpecies individualsSpecies]= speciesPredecesor(basedir2,poph, umbral,maxgen);
  ret.numSpecies = numSpecies;
  ret.individualsSpecies = individualsSpecies;
  save([basedir2 filesep 'miniSpeciationANTIGUA.mat'], 'ret');
end


ret = ret; %#ok<NASGU>

