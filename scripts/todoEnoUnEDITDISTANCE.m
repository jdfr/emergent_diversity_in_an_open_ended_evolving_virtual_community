% % %ESTO ES LO QUE HAY QUE EJECUTAR EN CADA NODO
% % 
% % %COMUN
% % -entrar en el ssh con la cuenta "josedavid"
% % 
% % %NODO 1 (T1)
% % -ssh n0001
% % -screen -d -r
% % (el matlab estará arrancado)
% % -todoEnUno('/datadisk/arboles/yahecho/G0.5F1/2009_Jul_29_12_16_39/1=0.0001_3', '/datadisk/arboles/porhacer/G0.5F1/2009_Jul_29_12_16_39/3', 80);
% % 
% % %NODO 2 (T1)
% % -ssh n0002
% % -screen -d -r
% % (el matlab estará arrancado)
% % -todoEnUno('/datadisk/arboles/yahecho/G0.5F1/2009_Jul_29_12_16_39/1=0.0001_7', '/datadisk/arboles/porhacer/G0.5F1/2009_Jul_29_12_16_39/7', 70);
% % 
% % %NODO 3 (T3) AL 11
% % -ssh n0003
% % -screen -d -r
% % (el matlab estará arrancado)
% % -todoEnUno('/datadisk/arboles/yahecho/G0.7F1/2009_Jul_29_12_20_37/1=0.0001_4', '/datadisk/arboles/porhacer/G0.7F1/2009_Jul_29_12_20_37/4', 120);
% % 
% % %NODO 4 (T3)
% % -ssh n0004
% % -screen -d -r
% % (el matlab estará arrancado)
% % -todoEnUno('/datadisk/arboles/yahecho/G0.7F1/2009_Jul_29_12_20_37/1=0.0001_2', '/datadisk/arboles/porhacer/G0.7F1/2009_Jul_29_12_20_37/2', 160);
% % 
% % %NODO 5 (T2) AL 12
% % -ssh n0005
% % -screen -d -r
% % (el matlab estará arrancado)
% % -todoEnUno('/datadisk/arboles/yahecho/G0.8F1/2009_Jul_29_12_21_07/1=0.0001_1', '/datadisk/arboles/porhacer/G0.8F1/2009_Jul_29_12_21_07/1', 150);
% % 
% % %NODO 9 (T2)
% % -ssh n0009
% % -screen -d -r
% % (el matlab estará arrancado)
% % -todoEnUno('/datadisk/arboles/yahecho/G0.8F1/2009_Jul_29_12_21_07/1=0.0001_2', '/datadisk/arboles/porhacer/G0.8F1/2009_Jul_29_12_21_07/2', 340);
% % 
% % %NODO 11 (T4)
% % -ssh n0011
% % -screen -d -r
% % (el matlab estará arrancado)
% % -todoEnUno('/datadisk/arboles/yahecho/G1F1/2009_Jul_29_14_22_12/1=0.0001_7', '/datadisk/arboles/porhacer/G1F1/2009_Jul_29_14_22_12/7', inf);
% % 
% % %NODO 12 (T4)
% % -ssh n0012
% % -screen -d -r
% % (el matlab estará arrancado)
% % -todoEnUno('/datadisk/arboles/yahecho/G1F1/2009_Jul_29_14_22_12/1=0.0001_5', '/datadisk/arboles/porhacer/G1F1/2009_Jul_29_14_22_12/5', inf);
% % 


function todoEnoUnEDITDISTANCE(basedir)
% todoEnoUnEDITDISTANCE('newdists\editdist\2009_Jul_29_12_16_39\1=0.001_1');
% todoEnoUnEDITDISTANCE('newdists\editdist\2009_Jul_29_14_22_12\1=0.001_1');
% todoEnoUnEDITDISTANCE('newdists\editdist\2009_Jul_29_12_21_07\1=0.001_1');
%
%HierarchicalClusteringEDITDISTANCE('rene/src/2009_Jul_29_12_21_07/2/1=0.001_1',[],423,442,true);
%todoEnoUnEDITDISTANCE('rene/src/2009_Jul_29_12_21_07/2/1=0.001_1');
%todoEnoUnEDITDISTANCE('rene/src/2009_Jul_29_14_22_12/5/1=0.001_1');
%HierarchicalClusteringNUMMUTS('rene/src/mutdist/2009_Jul_29_12_21_07/2/1=0.001_1', [], 135, 557, true); 
%HierarchicalClusteringNUMMUTS('rene/src/mutdist/2009_Jul_29_14_22_12/5/1=0.001_1', [], 76, 500, true); 
%todoEnoUnoNUMMUTSE('rene/src/mutdist/2009_Jul_29_12_21_07/2/1=0.001_1');
%todoEnoUnoNUMMUTSE('rene/src/mutdist/2009_Jul_29_14_22_12/5/1=0.001_1');

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
  HierarchicalClusteringEDITDISTANCE(basedir, poph,minGen,maxGen);
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

