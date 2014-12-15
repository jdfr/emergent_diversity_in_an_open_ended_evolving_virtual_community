function [files datenums dates datestrs isdir] = getDates(basedir)
%the format of datestr is suitable for the TOUCH command

%hist(dates, 200700000000:1000000:201108000000);

d                = dir(basedir);
dt               = [d.datenum]';
id               = [d.isdir]';
nm               = {d.name}';
nmsz             = cellfun('prodofsize', nm);
discard          = nmsz<=2;
discard(discard) = cellfun(@(x)all(x=='.'), nm(discard));
nd               = not(discard);
nm               = nm(nd);
dt               = dt(nd);
id               = id(nd);
datenums         = dt;

files    = cellfunc(@(x)[basedir filesep x], nm);
datestrs = arrayfunc(@(x)datestr(x, 'yyyymmddHHMM.SS'), dt);
dates    = cellfun(@str2double, datestrs);
dirs     = files(id);
isdir    = id;

a1 = cell(numel(dirs)+1,1);
a2 = cell(numel(dirs)+1,1);
a3 = cell(numel(dirs)+1,1);
a4 = cell(numel(dirs)+1,1);
a5 = cell(numel(dirs)+1,1);

for k=1:numel(dirs)
  [a1{k} a2{k} a3{k} a4{k} a5{k}] = getDates(dirs{k});
end
a1{end} = files;
a2{end} = datenums;
a3{end} = dates;
a4{end} = datestrs;
a5{end} = isdir;

files    = vertcat(a1{:});
datenums = vertcat(a2{:});
dates    = vertcat(a3{:});
datestrs = vertcat(a4{:});
isdir    = vertcat(a5{:});


