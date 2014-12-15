%para generar figuras del l-paper

base = 'lsystemdani\dataset';
subs = {'\G1F1\uno\1=0.001_3', 'G0.75F1\uno\1=0.001_1', 'G0.5F1\dos\1=0.001_1'};
keys = {'very harsh', 'harsh', 'mild'};

basedirs = cellfunc(@(x)[base filesep x], subs);

gnmstudy = cellfunc(@(x)load([x filesep 'reducedGenomeStudy.mat']), basedirs);
gnmstudy = cellfunc(@(x)x.res, gnmstudy);
pophs    = cellfunc(@(x)load([x filesep 'poph.mat']), basedirs);
pophs    = cellfunc(@(x)x.poph, pophs);
complexs = cellfunc(@(x)load([x filesep 'compression_9.mat']), basedirs);
complexs = cellfunc(@(x)x.compression, complexs);

lastGens = cellfun(@(x)max(x.generation), pophs);
lastPNames = arrayfunc(@(x)sprintf('P%03g', x), lastGens);
lastPs = cellfunc(@(x, y)load([x filesep y '.mat']), basedirs, lastPNames);
lastPs    = cellfunc(@(x, y)x.(y), lastPs, lastPNames);

showPopScript(1.2);

data1 = load('lsystemdani\dataset\T2dataset.mat');
data2 = load('lsystemdani\dataset\T2_G1F1_3.mat');
data3 = load('lsystemdani\dataset\T2_G0.5F1_2.mat');
data = data1.data;
data(1) = data2.data;
data(3) = data3.data;
clear data1 data2 data3;
{data.name}'
scriptShowRobustez(data);

figurasC = load('lsystemdani\dataset\figurasC_new_G1F1=3.mat');
figurasC = figurasC.figurasC;
figurasR = load('lsystemdani\dataset\figurasR_0.5_dos.mat');
figurasR = figurasR.figurasR;

makeDiversityFigures(figurasC, figurasR);

anglesmall = 'trees\il2\bestia\variable_pop_retune\randomSearch\2009_Jun_20_13_18_45_A=04.6_T=2_O=462.2_G=0.768_F=0.473\1=0.0001_1\entornoRaster380.png';
anglelarge='trees\il2\bestia\variable_pop_retune\randomSearch\2009_Jun_21_00_21_24_A=41.2_T=2_O=115.6_G=0.821_F=0.099\1=0.0001_1\entornoRaster400.png';
makeComparisonGenomes(gnmstudy{[3 2 1]});%, filename);

vals = [pophs([3 2 1]); complexs([3 2 1])];
scripShowCompression(vals{:});%, fname)

scripShowGL(pophs{[3 2 1]}, 1); scripShowGL(pophs{[3 2 1]}, 2);

scripShowGLGL(pophs{[3 2 1]});

% scripShowCompressionVsGnmLen(vals{:}, 1);%, fname);
% scripShowCompressionVsGnmLen(vals{:}, 2);%, fname);

clear vals;

showTreesAndComplexity(pophs{2}, lastGens(2), lastPs{2}, complexs{2});
