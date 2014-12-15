function res = estimateComputationalCostPerGeneration(basedir, poph)
%to show the curves of computational cost per generation

if not(exist('poph', 'var')) || isempty(poph)
  poph = loadSimulation(basedir);
end

ming = min(poph.generation);
maxg = max(poph.generation);

info = reapTimeStats(basedir,maxg,1,false);

incrementalTimes = [info.totalClusterReaped]';

gens = (1:numel(incrementalTimes));
compTimes = [incrementalTimes(1); diff(incrementalTimes)];

[gens2 difftimes2] = getTimesWithPFiles(basedir, ming, maxg);

if not(isfield(poph, 'pixel'))
  poph = getPixels(basedir, poph);
end

gens0 = (ming:maxg)';
totalsp   = zeros(size(gens0));
parcialsp = zeros(size(gens0));
totalsb   = zeros(size(gens0));
parcialsb = zeros(size(gens0));
for k=1:numel(gens0)
  g = gens0(k);
  thisg = find(poph.generation==g);
  totalsp(k) = sum(poph.pixel(thisg));
  totalsp(k) = sum(poph.nbranches(thisg));
  changed = not(strcmp('=', poph.tree.change(thisg)));
  parcialsp(k) = sum(poph.pixel(thisg(changed)));
  parcialsp(k) = sum(poph.nbranches(thisg(changed)));
end

figure;
plot(gens, compTimes, 'b', gens2, difftimes2, 'r', gens0, totalsp*4, 'g', gens0, parcialsp*4, 'k');
legend({'CPU times for tree renderization (not including fitness calculation)', 'total time for generation calculation', 'total bytes', 'parcial bytes'}); grid on;

res = struct('time1', struct('x', gens, 'y', compTimes), ...
             'time2', struct('x', gens2, 'y', difftimes2), ...
             'pixeltot', struct('x', gens0, 'y', totalsp), ...
             'pixelpar', struct('x', gens0, 'y', parcialsp), ...
             'nbranchtot', struct('x', gens0, 'y', totalsb), ...
             'nbranchpar', struct('x', gens0, 'y', parcialsb));



function [gens difftimes] = getTimesWithPFiles(basedir, ming, maxg)

datenums = zeros(maxg-ming+1,1);
k=1;
for g=ming:maxg
  pn = sprintf('P%03g', g);
  d = dir([basedir filesep pn '.mat']);
  datenums(k) = d.datenum;
  k=k+1;
end
datenums = datenums*86400;
gens = ((ming+1):maxg)';
difftimes = diff(datenums);



function [x y] = functionFittedToPixelTotX4(n)

a =  2.287e-005;%  (3.958e-006, 4.177e-005)
b =       5.097;%  (4.931, 5.264)
c =  1.933e+004;%  (3483, 3.518e+004)

x = (1:n)';
y = a*x.^b+c;
figure; plot(x,y); grid on; title('extrapolation for total number of bytes per generation');


function [x, y] = functionFittedToTime1(n)

 a =  1.481e-029;%  (-1.05e-028, 1.346e-028)
 b =       15.11;%  (13.49, 16.73)
 c =           0;%  (-116.8, 116.8)

x = (1:n)';
y = a*x.^b+c;
figure; plot(x,y); grid on; title('extrapolation for total computational time');



