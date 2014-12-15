function ramas(basedir,finaldir,nG,numIndividuals)

% nG es el numero de generaciones
% numIndividuals es el numero de individuos en la poblacion

pop = reapPopulationHistory(basedir);

for Generation = 1 : nG
    for i = 1: numIndividuals
        axiom = pop.individual{((Generation)*numIndividuals)+i};
        ramasgema = ls2(axiom, 'canonical');
        if ~isempty(ramasgema)
          ramasgema = [round(ramasgema(:,1:4)) ramasgema(:,5)]; %#ok<NASGU>
        end
        tosave=['ramas' num2str(Generation,'%02g') num2str(i,'%03g')];
        eval([tosave ' = ramasgema;']);
        save([finaldir filesep tosave],tosave);
        clear(tosave, 'ramasgema');
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = reapPopulationHistory(basedir)

fpob  = fopen([basedir filesep 'poblacion.txt'], 'r');
datos = textscan(fpob, '%d %d %d %d %d %s');
fclose(fpob);

p = struct('generation', [], ...
           'rangeid', [], ...
           'index', [], ...
           'fitness', [], ...
           'position', [], ...
           'individual', [] ...
           );
p.generation   = double(datos{1}); datos{1} = [];
p.rangeid      = double(datos{2}); datos{2} = [];
p.index        = double(datos{3}); datos{3} = [];
p.fitness      = double(datos{4}); datos{4} = [];
p.position     = double(datos{5}); datos{5} = [];
p.individual   = datos{6};         datos{6} = [];