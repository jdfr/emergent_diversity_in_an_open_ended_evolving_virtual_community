function addDispersions(namein, nameout)

allres = load(namein);
allres = allres.allres;
changed = false;

for k=1:numel(allres)
  d = dir(allres(k).dir);
  n = {d.name}';
  i = find(strcmp('dispersionsss.mat', n));
  if isempty(i)
    i = find(strcmp('dispersionss.mat', n));
    if isempty(i)
      i = find(strcmp('dispersions.mat', n));
    end
  end
  if isempty(i)
    continue;
  end
  data = load([allres(k).dir filesep n{i}]);
  data = data.data;
  allres(k).disps = data;
  changed = true;
end

if changed
  save(nameout, 'allres');
end