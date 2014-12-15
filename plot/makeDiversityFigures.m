function makeDiversityFigures(figurasC, figurasR)

fgs = {figurasC, figurasR};
% names = {{'zclusters.png', 'zphasecluter.png', 'zclusters2.png', 'zpopsize.png', 'zbiomass.png'}, {'zuniqueinds.png', 'zphaseuniqueinds.png', 'zuniqueinndsmorpho.png'}};
% nms = {'groups', 'unique individuals'};
names = {{'zclusters.png', 'zphasephylo1.png', 'zclusters2.png', 'zpopsize.png', 'zbiomass.png'}, {'zuniqueinds.png', 'zphasephylo10.png', 'zuniqueinndsmorpho.png'}};
nms = {'phylogroups (U=1)', 'phylogroups (U=10)'};

for k=1:numel(fgs)
%   imgs={};
%   h=figurassNuevo(fgs{k}, 1:3, 'ed', true, [0 0 1 0 0]);
%   set(get(gca, 'YLabel'), 'string', [nms{k} ' (edit dist.)']);
%   imgs{end+1,1} = imcapture(h, 'all', 600); close(h);
%   h=figurassNuevo(fgs{k}, 1:3, 'edred', true, [0 0 1 0 0]);
%   set(get(gca, 'YLabel'), 'string', [nms{k} ' (red. edit dist.)']);
%   imgs{end+1,1} = imcapture(h, 'all', 600); close(h);
%   h=figurassNuevo(fgs{k}, 1:3, 'morpho', true, [0 0 1 0 0]);
%   set(get(gca, 'YLabel'), 'string', [nms{k} ' (Jaccard distance)']);
%   imgs{end+1,1} = imcapture(h, 'all', 600); close(h);
%   img = vertcat(imgs{:});
%   imwrite(img, names{k}{1});
  
  imgs={};
  h=figurassNuevo(fgs{k}, 1:3, 'ed', true, [1 0 0 0 0]);
  set(get(gca, 'YLabel'), 'string', [nms{k} ' (edit dist.)']);
  imgs{end+1,1} = imcapture(h, 'all', 600); close(h);
  h=figurassNuevo(fgs{k}, 1:3, 'edred', true, [1 0 0 0 0]);
  set(get(gca, 'XLabel'), 'string', [nms{k} ' (red. edit dist.)']);
  imgs{end+1,1} = imcapture(h, 'all', 600); close(h);
  h=figurassNuevo(fgs{k}, 1:3, 'morpho', true, [1 0 0 0 0]);
  set(get(gca, 'XLabel'), 'string', [nms{k} ' (Jaccard distance)']);
  imgs{end+1,1} = imcapture(h, 'all', 600); close(h);
  img = vertcat(imgs{:});
  imwrite(img, names{k}{2});

%   if k==2
%     h=figurassNuevo(fgs{k}, 1:3, 'morpho', true, [0 0 1 0 0]);
%     set(get(gca, 'YLabel'), 'string', [nms{k} ' (Jaccard distance)']);
%     img = imcapture(h, 'all', 600); close(h);
%     imwrite(img, names{k}{3});
%   end
%   
%   
%   if k==1
%   imgs={};
%     h=figurassNuevo(fgs{k}, 1:3, 'ed', true, [0 0 1 0 0]);
%     set(get(gca, 'YLabel'), 'string', [nms{k} ' (edit dist.)']);
%     imgs{end+1,1} = imcapture(h, 'all', 600); close(h);
%     h=figurassNuevo(fgs{k}, 1:3, 'edred', true, [0 0 1 0 0]);
%     set(get(gca, 'YLabel'), 'string', [nms{k} ' (red. edit dist.)']);
%     imgs{end+1,1} = imcapture(h, 'all', 600); close(h);
%     img = vertcat(imgs{:});
%     imwrite(img, names{k}{3});
%   end
%   
%   if k==1
%     h=figurassNuevo(fgs{k}, 1:3, 'edred', true, [0 0 0 1 0]);
%     img = imcapture(h, 'all', 600); close(h);
%     imwrite(img, names{k}{4});
% 
%     h=figurassNuevo(fgs{k}, 1:3, 'edred', true, [0 0 0 0 1]); set(gca, 'YLim', [0 120000]); legend('Location', 'north');
%     img = imcapture(h, 'all', 600); close(h);
%     imwrite(img, names{k}{5});
%   end
end
