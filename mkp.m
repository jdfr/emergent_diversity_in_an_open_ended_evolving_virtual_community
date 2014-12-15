function zz=mkp(rmdirs)
%make a recursive addpath, with some subfolders excluded if the argument is
%present
%example: unixpth({'newsim\remove', 'newsim\rm1\rm2'})

p    = mfilename('fullpath');
last = find(p==filesep, 1, 'last');
pth  = genpath(p(1:(last-1)));
if nargin>0
  if ischar(rmdirs)
    rmdirs = {rmdirs};
  end
  for k=1:numel(rmdirs)
    matchs = strfind(pth, rmdirs{k});
    pth = pth(cellfun('prodofsize', matchs)==0);
  end
end
pth  = [pth, repmat({pathsep}, size(pth))];
pth  = reshape(pth', 1, []);
pth  = horzcat(pth{:});
addpath(pth);
if nargout>0
  zz=pth;
end
