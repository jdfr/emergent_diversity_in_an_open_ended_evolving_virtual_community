function getAllFamilyPhotos(basedir, otro)

d = dir(basedir);
for k=1:numel(d)
  if (d(k).isdir) && (d(k).name(1)~='.')
    doit([basedir filesep d(k).name], otro);
  end
end

function doit(basedir, otro)

d=dir(basedir);

nam = '';
num = -inf;

for k=1:numel(d)
  if (d(k).isdir)
    if (d(k).name(1)~='.')
      doit([basedir filesep d(k).name], otro);
    end
  elseif (numel(d(k).name)>13) && ...
         (all(lower(d(k).name(1:13))=='entornoraster')) && ...
         (all(lower(d(k).name(end-3:end))=='.png'))
    newnum = str2num(d(k).name(14:end-3));
    if newnum>num
      nam = d(k).name;
    end
  end
end

if ~isempty(nam)
  newnam = [otro filesep basedir(~ismember(basedir, [':' filesep])) nam];
  copyfile([basedir filesep nam], newnam);
end

fprintf('finishing %s\n', basedir);