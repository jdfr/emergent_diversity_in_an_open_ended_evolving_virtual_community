function drawAllIndividuals(basedir,finaldir,Generation,numIndividuals,ejes)

% basedir        = directorio donde esta el fichero poblacion.txt.
% finaldir       = directorio en el que quiero situar los resultados.
% Generation     = generacion de la que quiero obtener los resultados.
% numIndividuals = numero de individuos en la poblacion.
% ejes           = [xmin xmax ymin ymax], fijo los axis de todas las figuras que se van a pintar, para que
%                  todos tengan el mismo escalado.

% pinta a cada individuo individualmente 

pop = reapPopulationHistory(basedir);
f = fopen([finaldir,'IndGen.txt'],'a');


for i = 1 : numIndividuals
    drawing.figsize = [2000 2000];
    drawing.FontSize = 8;
    fig = figure('Visible', 'off', 'Units', 'pixels', 'PaperUnits', 'points');
    axis(ejes)
    axis equal

    i
    axiom = pop.individual{((Generation)*numIndividuals)+i};
    fit = pop.fitness(((Generation)*numIndividuals)+i);
    fprintf(f,'%3d  %7.4f  %s\n', i, fit, axiom);
    ramas = ls2(axiom, 'canonical');

    title(sprintf('%g',i),'FontSize',drawing.FontSize);
    drawtree(ramas, 0, [0 0 0]);
    axis off

    drawnow;
    saveas(fig,[finaldir,'Individuo' num2str(i, '%03g')],'png');
    close(fig);
end
 fclose all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = reapPopulationHistory(basedir)

fpob  = fopen([basedir filesep 'poblacion.txt'], 'r');
datos = textscan(fpob, '%d %d %d %f %d %s');
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

