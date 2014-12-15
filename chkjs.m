jm = findResource('scheduler', 'type', 'torque');
js = findJob(jm);
nt = zeros(size(js)); for k=1:numel(js); nt(k) = numel(js(k).Tasks); end; joba = js(nt>1); js = js(nt==1); joba, js
zs=cell(size(js)); ze = zs; for k=1:numel(js); zs{k}=js(k).State; ze{k} = js(k).Tasks.ErrorMessage; if isempty(ze{k}); ze{k} = ''; end; end; [zu za zb] = unique(zs); zf = [zu, arrayfunc(@(x)sum(zb==x), (1:numel(zu))')], [qa qb qc] = unique(ze); qf=[qa, arrayfunc(@(x)sum(qc==x), (1:numel(qa))')]
xa=find(strcmp('queued', zs)); for k=1:numel(xa); fprintf([mat2str(xa(k)) ' ' any2str(js(xa(k)).Tasks.InputArguments) '\n']);end
find(cellfun('prodofsize', ze)>0)
%z=getAllDirs('rene/dataset/test110NB'); generateDispersionStats('torque', z(randperm(numel(z))), ['rene/src'], true);

%to be executed in the nodes which are now down
% cd /usr/lib
% ln -s /opt/matlabr2008a/sys/opengl/lib/glnxa64/libGL.so.1.5.060001 libGL.so.1
% ln -s /opt/matlabr2008a/sys/opengl/lib/glnxa64/libGLU.so.1.3.060001 libGLU.so.1
% exit

%PARECE QUE TORQUE NO HACE BALANCEO DE CARGA EN IAMUS. PREGUNTAR A JULIO

%doForSeveralSimulationsParallel('torque', makedummy(@resumePopulation), {}, getAllDirs('rene/dataset/testLONGOLD2'), 'rene/src', false, @(x)randperm(numel(x)), struct('tag', 'touchPob'))
