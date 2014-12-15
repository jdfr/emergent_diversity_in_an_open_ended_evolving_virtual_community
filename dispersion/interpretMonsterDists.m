function pophs = interpretMonsterDists(res, space, pophs)
%show graphic results from multidimensional scaling of giant distance
%matrices, simulation-wise.

n=200;
popindexes = res.popindexes;
gens       = res.gens;
ygrid      = linspace(min(space), max(space), n);

switch size(space,2)
  case 1
    for k=1:numel(popindexes)
      gs = cellfunc(@(x, y) repmat(y, size(x)), popindexes{k}, array2cell(gens{k}));
      gs = vertcat(gs{:});
      inds = vertcat(popindexes{k}{:});
      if nargout>0
        if not(all((inds>0)==(pophs{k}.nPixelsInSoil==0)))
          error('inconsistent!!!');
        end
        metric = zeros(size(inds));
        metric(inds>0) = space(inds(inds>0));
        pophs{k}.mdscale_sstress = metric;
      end
      toDelete = inds==0;
      inds(toDelete) = [];
      gs(toDelete) = [];
      metric = space(inds);
      %ygrid = linspace(min(metric), max(metric), n);
      h = showLogHist3(gs, 'generations', min(gens{k}):max(gens{k}), metric, 'dissimilarity approximation', ygrid, 'evolution of dissimilarity', [1 1 1; jet(255)], false);
    end
  otherwise
    error('not implemented!!!!!');
end

