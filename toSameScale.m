function [sizesO scalesI] = toSameScale(files, dests, scale)

%file{n}='_G1.00_N009_1000_sc22100X213.png';

regul = '([0-9]+)X([0-9]+)';%.png';

if ischar(files)
  basedir = files;
  d = dir(files);
  d = {d.name}';
  ok = cellfun(@(x)isempty(regexp(x, regul, 'once')), d);
  names = d(not(ok));
  files = cellfun(@(x)[basedir filesep x], names, 'uniformoutput', false);
else
  [nevermind names] = cellfun(@fileparts, files, 'uniformoutput', false);
end

if ischar(dests)
  basedir = dests;
  dests = cellfun(@(x)[basedir filesep x], names, 'uniformoutput', false);
end

sizesO = zeros(numel(files), 2);
%scales = zeros(numel(files), 1);
scalesI = zeros(numel(files), 2);

for k=1:numel(files)
  tokens = regexp(names{k}, regul, 'tokens');
  [img map] = imread(files{k});
  orig = [str2double(tokens{1}{1}) str2double(tokens{1}{2})];
  sizesO(k,:) = orig;
  actual = [size(img,2) size(img,1)];
  scaleA = actual./orig;
  scalesI(k,:) = scaleA;
  newscale = scale/mean(scaleA);
  %scales(k) = newscale;
  [img map] = imresize(img, map, newscale);
  imwrite(img, map, dests{k});
end

