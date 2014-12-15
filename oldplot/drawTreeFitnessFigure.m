function [] = drawTreeFitnessFigure()
% Función que dibuja la figura de las semillas, árboles y rayos de luz en la cuadrícula.

imageFileName{1} = 'arbol_381561.png';
maxLeavesHeightImage{1} = [40 -1 -1 -1 -1 -1 38 -1 -1 -1 -1 -1 46 -1 -1 -1 31 -1 44 -1 -1 39 -1 -1 42 34 -1 37 -1 -1 -1 42 -1 -1 21 27 -1 40 -1 -1 29 35 -1 -1 17 -1 -1 -1 -1 -1 25];
imageFileName{2} = 'arbol_388382.png';
maxLeavesHeightImage{2} =  [24 -1 -1 -1 29 -1 -1 34 -1 31 -1 -1 -1 36 -1 -1 -1 -1 -1 38 -1 -1 -1 -1 -1 40 -1 -1 -1 -1 -1 38 -1 -1 -1 -1 -1 36 -1 -1 -1 31 -1 34 -1 -1 29 -1 -1 -1 24];
imageFileName{3} = 'arbol_381946.png';
maxLeavesHeightImage{3} = [19 6 -1 -1 24 -1 -1 -1 29 -1 16 34 -1 21 -1 -1 -1 36 -1 -1 -1 -1 -1 28 -1 -1 -1 -1 -1 30 -1 -1 -1 -1 -1 28];

treeLocations = [1 51 101];

grillColor = [0.7 0.7 0.7];
lightColor = [150 208 255]/255;

topMargin = 5;

%dateNow = datestr(

for(k=1:3)
    m{k} = imread(imageFileName{k});
    dims(k,:) = size(m{k});
end

dim(1) = max(dims(:, 1)) + topMargin;
dim(2) = treeLocations(3) + dims(3,2) + 1;

fullImage = zeros(dim, 'uint8');

maxLeavesHeight = -ones(1, dim(2));
for(k=1:3)
    for(i=1:dims(k,1))
        for(j=1:dims(k,2))
            ci = i+dim(1)-dims(k,1);
            cj = j+treeLocations(k);
            if(~fullImage(ci, cj))
                fullImage(ci, cj) = m{k}(i,j);
                maxLeavesHeight(cj) = max(maxLeavesHeight(cj), maxLeavesHeightImage{k}(j));
            end
        end
    end
end

figure(1);
clf;
%subplot(3,1,1)
roots = fullImage;
roots(1:end-3,:) = 0;
hold on;
plotGrills(dim, topMargin, grillColor);
alphaChannel = roots;
alphaChannel(roots ~= 0) = 255;
image(roots(end:-1:1, :), 'AlphaData', alphaChannel(end:-1:1, :));
colormap([1 1 1; 0 0 0; 1 0 0]);
hold off;
axis off;
axis equal;
fileName = 'treeFitnessCalcFigure.A.png';
print(fileName, '-dpng', '-r600');
system(['convert ' fileName ' -trim ' fileName]);


clf;
%subplot(3,1,2)
hold on;
plotGrills(dim, topMargin, grillColor);
alphaChannel = fullImage;
alphaChannel(fullImage ~= 0) = 255;
image(fullImage(end:-1:1, :), 'AlphaData', alphaChannel(end:-1:1, :));
colormap([1 1 1; 0 0 0; 1 0 0]);
hold off;
axis off;
axis equal;
fileName = 'treeFitnessCalcFigure.B.png';
print(fileName, '-dpng', '-r600');
system(['convert ' fileName ' -trim ' fileName]);


clf;
%subplot(3,1,3)
hold on;
plotGrills(dim, topMargin, grillColor);

% Light
for(i=1:dim(2))  
    arrow('Start', [i dim(1)+0.5], 'Stop', [i maxLeavesHeight(i)+1.5], 'LineWidth', 0.4, 'Length', 7, 'TipAngle', 10, 'EdgeColor', lightColor, 'FaceColor', lightColor); 
end

alphaChannel = fullImage;
alphaChannel(fullImage ~= 0) = 255;
image(fullImage(end:-1:1, :), 'AlphaData', alphaChannel(end:-1:1, :));
colormap([1 1 1; 0 0 0; 1 0 0]);

hold off;
axis off;
axis equal;
fileName = 'treeFitnessCalcFigure.C.png';
print(fileName, '-dpng', '-r600');
system(['convert ' fileName ' -trim ' fileName]);



treeSeedsLocations = [5 63 120];
fullImage = zeros(dim, 'uint8');
for(i=1:length(treeSeedsLocations))
    for(j=dim(1)-2:dim(1))
        fullImage(j, treeSeedsLocations(i)) = 1;
    end
end

clf;
%subplot(3,1,2)
hold on;
plotGrills(dim, topMargin, grillColor);
alphaChannel = fullImage;
alphaChannel(fullImage ~= 0) = 255;
image(fullImage(end:-1:1, :), 'AlphaData', alphaChannel(end:-1:1, :));
colormap([1 1 1; 0 0 0; 1 0 0]);
hold off;
axis off;
axis equal;
fileName = 'treeFitnessCalcFigure.D.png';
print(fileName, '-dpng', '-r600');
system(['convert ' fileName ' -trim ' fileName]);




function [] = plotGrills(dim, topMargin, grillColor)

% Vert grill
for(i=linspace(1, dim(2)+1.02, dim(2)+1))
    plot([i i]-0.5, [0.5 dim(1)+0.5], 'LineWidth', 0.1, 'Color', grillColor); 
end

% Horiz grill
for(i=1:dim(1)+1) 
    plot([0.5 dim(2)+0.5], [i i]-0.5, 'LineWidth', 0.1, 'Color', grillColor); 
end
