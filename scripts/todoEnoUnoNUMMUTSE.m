function todoEnoUnoNUMMUTSE(basedir)
% todoEnoUnoNUMMUTSE('newdists\mutdist\2009_Jul_29_12_16_39\1=0.001_1');
% todoEnoUnoNUMMUTSE('newdists\mutdist\2009_Jul_29_14_22_12\1=0.001_1');
% todoEnoUnoNUMMUTSE('newdists\mutdist\2009_Jul_29_12_21_07\1=0.001_1');

if ~exist([basedir filesep 'poph.mat'], 'file')
    poph          = loadSimulation(basedir);
  poph.minGenome  = cellfun(@sanitizeString, poph.genome, 'uniformoutput', false);
  save([basedir filesep 'poph.mat'], 'poph');
else
  fprintf('>>>>ESTOY EN LOADSIMULATION 3<<<<\n');
  load([basedir filesep 'poph.mat'], 'poph');
end
  poph            = rmfield(poph, 'minGenome');

maxGen = max(poph.generation);

% if ~exist([basedir filesep 'P' num2str(maxGen, '%03g') '.mat'], 'file')
%   fprintf('>>>>ESTOY EN CALCULATEPMAT<<<<\n');
%   calculatePmat('embedded', poph, [0 maxGen], basedir);
% end

minGen = 1;

if ~exist([basedir filesep 'indicesClasificados' num2str(maxGen, '%03g') '.mat'], 'file')
  fprintf('>>>>ESTOY EN HierarchicalClustering<<<<\n');
  HierarchicalClusteringNUMMUTS(basedir, poph,minGen,maxGen);
end

if ~exist([basedir filesep 'numGroups' num2str(maxGen, '%03g') '.mat'], 'file')
  fprintf('>>>>ESTOY EN numberOfSpecies<<<<\n');
  numberOfSpecies(basedir, minGen,maxGen);
end

% if ~exist([nuevosubdir filesep 'comLZfigureLongitud' num2str(maxGen, '%03g') '.mat'], 'file')
%   fprintf('>>>>ESTOY EN complexityProto<<<<\n');
%   complexityProto(nuevosubdir,minGen,maxGen);
% end
% 
if ~exist([basedir filesep 'comEspeciesSig.mat'], 'file')
  fprintf('>>>>ESTOY EN relacionarSpecies<<<<\n');
  [especies ] = relacionarSpeciesNOCOMPLEXITY({poph, basedir},maxGen); %#ok<NASGU>
  save([basedir filesep 'especies.mat'],        'especies');
  clear especies ;
end
% 
% clear poph;
% 
% if ~exist([nuevosubdir filesep 'zdata.mat'], 'file')
%   fprintf('>>>>ESTOY EN loadSpeciation<<<<\n');
%   zdata = loadSpeciation(nuevosubdir, maxGen, 'especies'); %#ok<NASGU>
%   save([nuevosubdir filesep 'zdata.mat'], 'zdata');
%   clear zdata;
% end

