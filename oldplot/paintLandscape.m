function paintLandscape(basedir, generation, N)
% Plot using the regular plot command the landscape (forest) of a simulation. 
% Inputs:
%   -basedir: the base directory of the simulation (i.e., the directory
%             where poblacion.txt and arbol.txt are placed)
%   -generation: the generation number to plot
%   -N: landscape width

pop = reapPopulationHistory(basedir);
numIndividuals = sum(pop.generation == 0);

for i = 1 : numIndividuals
    axiom = pop.individual{((generation)*numIndividuals)+i};
    ramas = ls2(axiom, 'canonical');
    randN = pop.posicion(((generation)*numIndividuals)+i);
    if ~isempty(ramas)
        r = ramas(:,1:4);
        r(:,1:2) = r(:,1:2)+randN;
        t = r;
        t(:,5) = ramas(:,5);
        drawtree(t,0);
        hold on
    end
end
plot([0 N], [0 0], 'k');
drawnow
print('-dpng', [basedir,'entorno.png']);



    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = reapPopulationHistory(basedir)

fpob  = fopen([basedir filesep 'poblacion.txt'], 'r');
datos = textscan(fpob, '%d %d %d %d %d %s');
fclose(fpob);

p = struct('generation', [], ...
           'rangeid', [], ...
           'index', [], ...
           'fitness', [], ...
           'posicion', [], ...
           'individual', [] ...
           );
p.generation   = double(datos{1}); datos{1} = [];
p.rangeid      = double(datos{2}); datos{2} = [];
p.index        = double(datos{3}); datos{3} = [];
p.fitness      = double(datos{4}); datos{4} = [];
p.posicion     = double(datos{5}); datos{5} = [];
p.individual   = datos{6};         datos{6} = [];
