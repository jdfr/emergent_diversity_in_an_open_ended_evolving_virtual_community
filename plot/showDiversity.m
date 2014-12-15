function showDiversity

figs = load('lsystemdani\figurasrene\phylofiguras_0.5dos.mat');

figurasC = figs.figurasC;

figurassNuevo(figurasC, 1:3, 'morpho', true); set(get(gca, 'ylabel'), 'string', 'diversity (Jaccard distance)');
img = imcapture(gcf, 'all', 600);
close all hidden;
imwrite(img, 'diversityMorpho.png');

figurassNuevo(figurasC, 1:3, 'ed', true); set(get(gca, 'ylabel'), 'string', 'diversity (edit distance)');
img = imcapture(gcf, 'all', 600);
close all hidden;
imwrite(img, 'diversityEdit.png');

figurassNuevo(figurasC, 1:3, 'edred', true); set(get(gca, 'ylabel'), 'string', 'diversity (reduced edit distance)');
img = imcapture(gcf, 'all', 600);
close all hidden;
imwrite(img, 'diversityEditRed.png');
