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


function todoEnUno(basedir, nuevodir, numGenHastaAhora)

nuevoNumGenerations = 550-numGenHastaAhora;

nuevo_tooLongEvTime = 18000;%7200;
nuevo_stopIfTooTall = 20000;
nuevo_ensayos       = 1;

numero_noseusa      =  0.001;

nuevosubdir = [nuevodir filesep '1=' mat2str(numero_noseusa) '_1'];

if (nuevoNumGenerations>0) && (~exist([nuevosubdir filesep 'graph_fitnessHistogramByGen.png'], 'file'))
  fprintf('>>>>ESTOY EN MAIN<<<<\n');
  main('embedded', {'', 'continue'}, basedir, nuevodir, numero_noseusa, nuevoNumGenerations, 'tooLongEvTime', {nuevo_tooLongEvTime}, 'stopIfTooTall', {nuevo_stopIfTooTall}, 'ensayos', {nuevo_ensayos});
end

if ~exist([nuevosubdir filesep 'poph.mat'], 'file')
  if (nuevoNumGenerations>0)
    fprintf('>>>>ESTOY EN LOADSIMULATION 1<<<<\n');
    poph          = loadSimulation({basedir, nuevosubdir});
  else
    fprintf('>>>>ESTOY EN LOADSIMULATION 2<<<<\n');
    poph          = loadSimulation(basedir);
  end
  poph.minGenome  = cellfun(@sanitizeString, poph.genome, 'uniformoutput', false);
  save([nuevosubdir filesep 'poph.mat'], 'poph');
else
  fprintf('>>>>ESTOY EN LOADSIMULATION 3<<<<\n');
  load([nuevosubdir filesep 'poph.mat'], 'poph');
end
  poph            = rmfield(poph, 'minGenome');

maxGen = max(poph.generation);

if ~exist([nuevosubdir filesep 'P' num2str(maxGen, '%03g') '.mat'], 'file')
  fprintf('>>>>ESTOY EN CALCULATEPMAT<<<<\n');
  calculatePmat('embedded', poph, [0 maxGen], nuevosubdir);
end

minGen = 1;

if ~exist([nuevosubdir filesep 'indicesClasificados' num2str(maxGen, '%03g') '.mat'], 'file')
  fprintf('>>>>ESTOY EN HierarchicalClustering<<<<\n');
  HierarchicalClustering(nuevosubdir, poph,minGen,maxGen);
end

if ~exist([nuevosubdir filesep 'numGroups' num2str(maxGen, '%03g') '.mat'], 'file')
  fprintf('>>>>ESTOY EN numberOfSpecies<<<<\n');
  numberOfSpecies(nuevosubdir, minGen,maxGen);
end

if ~exist([nuevosubdir filesep 'comLZfigureLongitud' num2str(maxGen, '%03g') '.mat'], 'file')
  fprintf('>>>>ESTOY EN complexityProto<<<<\n');
  complexityProto(nuevosubdir,minGen,maxGen);
end

if ~exist([nuevosubdir filesep 'comEspeciesSig.mat'], 'file')
  fprintf('>>>>ESTOY EN relacionarSpecies<<<<\n');
  [especies comEspeciesAnte comEspeciesSig] = relacionarSpecies({poph, nuevosubdir},maxGen); %#ok<NASGU>
  save([nuevosubdir filesep 'especies.mat'],        'especies');
  save([nuevosubdir filesep 'comEspeciesAnte.mat'], 'comEspeciesAnte');
  save([nuevosubdir filesep 'comEspeciesSig.mat'],  'comEspeciesSig');
  clear especies comEspeciesAnte comEspeciesSig;
end

clear poph;

if ~exist([nuevosubdir filesep 'zdata.mat'], 'file')
  fprintf('>>>>ESTOY EN loadSpeciation<<<<\n');
  zdata = loadSpeciation(nuevosubdir, maxGen, 'especies'); %#ok<NASGU>
  save([nuevosubdir filesep 'zdata.mat'], 'zdata');
  clear zdata;
end




% 1-5 y 9-11
% 
% T1 
% 
% srctree/resultados/systematicAlphaF1/G0.5F1/2009_Jul_29_12_16_39/1=0.0001_3
% 80
% 
% srctree/resultados/systematicAlphaF1/G0.5F1/2009_Jul_29_12_16_39/1=0.0001_7
% 70
% 
% T3
% 
% srctree/resultados/systematicAlphaF1/G0.7F1/2009_Jul_29_12_20_37/1=0.0001_4
% 120
% 
% srctree/resultados/systematicAlphaF1/G0.7F1/2009_Jul_29_12_20_37/1=0.0001_2
% 160
% 
% T2
% 
% srctree/resultados/systematicAlphaF1/G0.8F1/2009_Jul_29_12_21_07/1=0.0001_1
% 150
% 
% srctree/resultados/systematicAlphaF1/G0.8F1/2009_Jul_29_12_21_07/1=0.0001_2
% 340
% 
% T4
% 
% srctree/resultados/systematicAlphaF1/G1F1/2009_Jul_29_14_22_12/1=0.0001_7
% 500
% 
% srctree/resultados/systematicAlphaF1/G1F1/2009_Jul_29_14_22_12/1=0.0001_5
% 500

