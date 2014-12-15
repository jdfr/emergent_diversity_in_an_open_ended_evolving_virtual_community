function h = showLogHist3(xaxis, xlab, xdivs, yaxis, ylab, yndivs, titname, colmap, useLogY)

puthalfticks=true;

minc = min(yaxis);
maxc = max(yaxis);
ming = min(xaxis);
maxg = max(xaxis);
if useLogY
  minc = log10(minc);
  maxc = log10(maxc);
  yaxis = log10(yaxis);
  cs = linspace(minc, maxc, yndivs);
else
  if isscalar(yndivs)
    cs = linspace(minc, maxc, yndivs);
  else
    cs = yndivs;
  end
end
if isempty(xdivs)
  gs = ming:maxg;
else
  if isscalar(xdivs)
    gs = linspace(ming, maxg, xdivs);
  else
    gs = xdivs;
  end
end
h = myhist3(xaxis, xlab, gs, [], yaxis, ylab, cs, realpow(10, cs), true, colmap);
if useLogY
  a = get(h, 'currentaxes');
  %yt = get(a, 'YTick');
  yt = ceil(minc):floor(maxc);
  if puthalfticks
    yt = [log10(realpow(10, yt')/2), yt']';
    yt = yt(:);
    yt = yt(2:end); %skip the first half
    ytceil = ceil(maxc);
    if ytceil~=maxc %don't do this if maxc is an integer
      ytceil2 = log10(realpow(10, ytceil)/2);
      if ytceil2<=maxc
        yt(end+1) = ytceil2;
      end
    end
  end
  set(a, 'YTick', yt, 'YTickLabel', realpow(10, yt));
end
if not(isempty(titname))
  title(titname);
end
% set(a, 'YScale', 'log');

