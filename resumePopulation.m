function resumePopulation(basedir)
%takes an old simulation dir and compresses the text file poblacion.txt by
%replacing genomes identical to the parent's with 'I'

%call this recursively for all simulations under some directory
%doForSeveralSimulationsParallel('torque', makedummy(@resumePopulation), {}, getAllDirs('rene/dataset/test110NB'), 'rene/src', false, @(x)randperm(numel(x)), struct('tag', 'touchPob'))  

%31.8Gb

tic; t1 = cputime;

fprintf('loading untouched poblacion.txt\n');
poph = loadSimulation(basedir);
fprintf('done!!!! %s, %s...\n', mat2str(toc), mat2str(cputime-t1));

dr = dir([basedir filesep 'poblacion.txt']);

ds = datestr(dr(1).datenum, 'yyyymmddHHMM.SS');

fil = fopen([basedir filesep 'poblacion.txtt'], 'w');

GNM = poph.genome;
MNG = poph.minGenome;
g = poph.generation;
pa = poph.tree.parents;
r = poph.rangeid;
i = poph.individual;
f = poph.fitness;
x = poph.xpos;
h = poph.height;
w = poph.width;
n = poph.nbranches;
l = poph.nleafs;
s = poph.nPixelsInSoil;
t = poph.nleafsOnTop;

last = -1;
for k=1:numel(g)
  gen = g(k);
  if gen~=last
    last = gen;
    fprintf('processing gen %d, %s, %s...\n', gen, mat2str(toc), mat2str(cputime-t1));
  end
  gnm = GNM{k};
  mng = MNG{k};
  p = pa(k);
  if p>0
    if strcmp(gnm, GNM{p})
      gnm = 'I';
    end
    if strcmp(mng, MNG{p})
      mng = 'I';
    end
  end
  fprintf(fil,'%d %d %d %s %d %d %d %d %d %d %d %s %s\n',g(k),r(k),i(k),mat2str(f(k)),x(k), h(k), w(k), n(k), l(k), s(k), t(k), mng, gnm);
end

fclose(fil);

delete([basedir filesep 'poblacion.txt']);

movefile([basedir filesep 'poblacion.txtt'], [basedir filesep 'poblacion.txt']);

system(['touch -t ' ds ' ' basedir filesep 'poblacion.txt']);

fprintf('finished!!! %s, %s...\n', mat2str(toc), mat2str(cputime-t1));
