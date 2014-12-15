function paintFitness0(basedir,generation,numIndividuals,N)

pop = reapPopulationHistory(basedir);


for i = 1 : numIndividuals
    i
    fitness = pop.fitness(((generation)*numIndividuals)+i);
    if fitness ~= 0

        axiom = pop.individual{((generation)*numIndividuals)+i};
        [ramas,nramas,alto, ancho, leaf,maxiy] = ls2(axiom, 'canonical'); %#ok<NASGU,ASGLU>
        randN = randint(1,1,[1,N-length(maxiy)]);

%         randN = pop.posicion(((generation)*numIndividuals)+i);
        if ~isempty(ramas)
            r = ramas(:,1:4);
            r(:,1:2) = r(:,1:2)+randN;
            t = r;
            t(:,5) = ramas(:,5);
            drawtree(t,0)
            hold on
        end
    end
end
plot([0 N], [0 0], 'k');
drawnow




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = reapPopulationHistory(basedir)

fpob  = fopen([basedir filesep 'poblacion.txt'], 'r');
datos = textscan(fpob, '%d %d %d %d %s');
fclose(fpob);

p = struct('generation', [], ...
           'rangeid', [], ...
           'index', [], ...
           'fitness', [], ...
           'individual', [] ...
           );
p.generation   = double(datos{1}); datos{1} = [];
p.rangeid      = double(datos{2}); datos{2} = [];
p.index        = double(datos{3}); datos{3} = [];
p.fitness      = double(datos{4}); datos{4} = [];
p.individual   = datos{5};         datos{5} = [];

% p = struct('generation', [], ...
%            'rangeid', [], ...
%            'index', [], ...
%            'fitness', [], ...
%            'posicion', [], ...
%            'individual', [] ...
%            );
% p.generation   = double(datos{1}); datos{1} = [];
% p.rangeid      = double(datos{2}); datos{2} = [];
% p.index        = double(datos{3}); datos{3} = [];
% p.fitness      = double(datos{4}); datos{4} = [];
% p.posicion     = double(datos{5}); datos{5} = [];
% p.individual   = datos{6};         datos{6} = [];

