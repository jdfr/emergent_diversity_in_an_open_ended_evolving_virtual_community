function showIndividualHistory(poph, generation, rangeid, individual, range)
%print individual statistics (fitness, etc.) throughout evolutionary
%history

gens          = max(min(poph.generation), range(1)):min(max(poph.generation), range(2));
absolute.fitness       = repmat(nan, size(gens));
absolute.height        = repmat(nan, size(gens));
absolute.width         = repmat(nan, size(gens));
absolute.nbranches     = repmat(nan, size(gens));
absolute.nleafs        = repmat(nan, size(gens));
absolute.nPixelsInSoil = repmat(nan, size(gens));
absolute.nleafsOnTop   = repmat(nan, size(gens));
relative.fitness       = repmat(nan, size(gens));
relative.height        = repmat(nan, size(gens));
relative.width         = repmat(nan, size(gens));
relative.nbranches     = repmat(nan, size(gens));
relative.nleafs        = repmat(nan, size(gens));
relative.nPixelsInSoil = repmat(nan, size(gens));
relative.nleafsOnTop   = repmat(nan, size(gens));

names = fieldnames(absolute);

indiv = poph.triplet2idx(generation, rangeid, individual);

nchar = 0;

while true
  idx = find(poph.idx==indiv);
  if isempty(idx)
    break;
  else
    gen                = find(gens==poph.generation(idx));
    
    for k=1:numel(names)
      name = names{k};
      absolute.(name)(gen)       = poph.(name)(idx);
      thisGen                    = poph.generation==gens(gen);
      mn                         = min(poph.(name)(thisGen));
      range                      = max(poph.(name)(thisGen))-mn;
      if range==0
        relative.(name)(gen)     = 0.5;
      else
        relative.(name)(gen)     = (absolute.(name)(gen)-mn)/range;
      end
    end
    indiv                        = poph.tree.idxA(poph.tree.idxD==indiv);
  end
  if nchar>0
    fprintf(repmat('\b', 1, nchar));
  end
  str   = sprintf('Processed generation %04d', gens(gen));
  fprintf(str);
  nchar = numel(str);
end
if nchar>0
  fprintf(repmat('\b', 1, nchar));
end

figure; hold on; stairs(gens, absolute.fitness, 'b');
xlabel('generations'); ylabel('fitness'); grid on;
% ax1 = gca;
% ax2 = axes('Position',get(ax1,'Position'),...
%            ...%'XAxisLocation','top',...
%            'YAxisLocation','right',...
%            'Color','none',...
%            'XColor','k','YColor','k');
% stairs(gens, relative.fitness, 'c', 'Parent', ax2);
         
figure; hold on; stairs(gens, absolute.height, 'b'); stairs(gens, absolute.width, 'r');
xlabel('generations'); legend({'height', 'width'}); grid on;
% ax1 = gca;
% ax2 = axes('Position',get(ax1,'Position'),...
%            ...%'XAxisLocation','top',...
%            'YAxisLocation','right',...
%            'Color','none',...
%            'XColor','k','YColor','k');
% stairs(gens, relative.height, 'c'); stairs(gens, relative.width, 'm', 'Parent', ax2);

figure; hold on; stairs(gens, absolute.nbranches, 'b'); stairs(gens, absolute.nleafs, 'r'); stairs(gens, absolute.nPixelsInSoil, 'g'); stairs(gens, absolute.nleafsOnTop, 'm');
xlabel('generations'); legend({'branches', 'leafs', 'pixels below soil', 'leafs on top'}); title('absolute values'); grid on;

% figure; hold on; stairs(gens, relative.nbranches, 'b'); stairs(gens, relative.nleafs, 'r'); stairs(gens, relative.nPixelsInSoil, 'g'); stairs(gens, absolute.nleafsOnTop, 'm');
% xlabel('generations'); legend({'branches', 'leafs', 'pixels below soil', 'leafs on top'}); title('relative values'); grid on;

