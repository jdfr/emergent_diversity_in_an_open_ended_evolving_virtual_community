function [mn mx logvalues,logyticks, yticks]= calculateLogScaleTicks(mn, mx, values, puthalfticks)

mn = log10(mn);
mx = log10(mx);
logvalues = log10(values);
yt = ceil(mn):floor(mx);
if puthalfticks
  yt = [log10(realpow(10, yt')/2), yt']';
  yt = yt(:);
  yt = yt(2:end); %skip the first half
  ytceil = ceil(mx);
  ytceil2 = log10(realpow(10, ytceil)/2);
  if ytceil2<=mx
    yt(end+1) = ytceil2;
  end
end
yt10 = realpow(10, yt);
logyticks=yt;
yticks = yt10;

