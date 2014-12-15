function scriptEvolutionOfComplexity

basedir = 'gemabestia\G0.7F0.5\2009_Jul_05_17_08_32\1=0.0001_4';

maxgen = 243;

state = rand('twister');
rand('twister', 0);


comEspeciesAnte243 = load([basedir filesep 'comEspeciesAnte243.mat']);
comEspeciesAnte243 = comEspeciesAnte243.comEspeciesAnte243;
comEspeciesSig243  = load([basedir filesep 'comEspeciesSig243.mat']);
comEspeciesSig243  = comEspeciesSig243.comEspeciesSig243;

h = showEvolutionOfComplexity(basedir, maxgen, comEspeciesAnte243, comEspeciesSig243);

rand('twister', state);

a = get(h, 'CurrentAxes');

set(a, 'YLim', [0 20000]);%, 'YScale', 'log');

