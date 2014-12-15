function [raster] = paintRasterLandscape(basedir, generation, N)
% Plot using a raster the landscape (forest) of a simulation. 
% Inputs:
%   -basedir: the base directory of the simulation (i.e., the directory
%             where poblacion.txt and arbol.txt are placed)
%   -generation: the generation number to plot
%   -N: landscape width

pop = reapPopulationHistory(basedir);
numIndividuals = sum(pop.generation == 0);

backgroundColor = [1 1 1];
floorColor = [0 0 0];
colors = rand(numIndividuals, 3);
colors(colors>0.8) = 0.8;
colors = [backgroundColor; floorColor; colors];

yOffset = 1;
raster = zeros(2, N); % both index 0 and 1 refers to the first color in the color map
raster(1, :) = 2; % floor

for i = 1 : numIndividuals
    axiom = pop.individual{((generation)*numIndividuals)+i};
    [ramas, nramas] = ls2(axiom, 'canonical');
    
    if ~isempty(ramas)
        minY = min(min(ramas(:, [3 4])));
        if(minY + yOffset < 1)
            raster = [zeros(1-minY-yOffset, N); raster]; % This is not called too much
            yOffset = 1-minY;
        end
        
        maxY = max(max(ramas(:, [3 4])));
        if(maxY + yOffset > size(raster,1))
            raster(end+1:maxY+yOffset, :) = 0; % This is not called too much
        end
        
        randN = pop.posicion(((generation)*numIndividuals)+i);
        ramas(:,[1 2]) = ramas(:,[1 2]) + randN;
        ramas(:,[3 4]) = ramas(:,[3 4]) + yOffset;
        for(j=1:nramas)
            raster = rasterLine(raster, ramas(j,1), ramas(j,2), ramas(j,3), ramas(j,4), i+2);
        end
    end
    
    %rasterTemp = [zeros(5, N); raster(end:-1:1,:); zeros(5, N)];
    %imwrite(full(rasterTemp), colors, sprintf('%s\\entornoRaster%03d.png',basedir, generation));
end

raster = [zeros(5, N); raster(end:-1:1,:); zeros(5, N)];

image(raster);
colormap(colors);
imwrite(full(raster), colors, sprintf('%s\\entornoRaster%03d.png',basedir, generation));


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

