function pophs = getPixels(dirs, pophs)
%from a simulation with PXXX.mat files, get number of pixels for each tree

for k=1:numel(dirs)
  d = dir([dirs{k} filesep 'P*.mat']);
  d = cellfun(@(x)str2double(x(2:4)), {d.name}');
  ming = min(d);
  maxg = max(d);
  pophs{k}.pixel = zeros(size(pophs{k}.generation));
  for g=ming:maxg
    fprintf('GEN %03g, %s\n', g, dirs{k});
    pn = sprintf('P%03g', g);
    P = load([dirs{k} filesep pn '.mat']);
    P = P.(pn);
    pixels = cellfun(@(x)numel(x{1}), P.raster);
    pophs{k}.pixel(pophs{k}.generation==g) = pixels;
  end
end



