function h=myhist3(val1, nom1, ran1, valran1, val2, nom2, ran2, valran2, dofig, colmap)
if isscalar(ran1)
  ran1 = linspace(min(val1), max(val1), ran1);
end
if isscalar(ran2)
  ran2 = linspace(min(val2), max(val2), ran2);
end
if isempty(valran1)
  valran1 = ran1;
end
if isempty(valran2)
  valran2 = ran2;
end
nlabs = 10;
[n c] = hist3([val1, val2], {ran1, ran2});
nlog = log10(n+1);
dopretty = true;
puthalfticks=true;
if dopretty
  minlog = min(min(nlog(nlog>0)));
  maxlog = max(nlog(:));
  minlogverdad = log10(realpow(10, minlog)-1);
  maxlogverdad = log10(realpow(10, maxlog)-1);
  valsverdad = ceil(minlogverdad):floor(maxlogverdad);
  if puthalfticks
    valsverdad = [log10(realpow(10, valsverdad')/2), valsverdad']';
    valsverdad = valsverdad(:);
    valsverdad = valsverdad(2:end); %skip the first half
    maxlogverdadceil = ceil(maxlogverdad);
    maxlogverdadceil2 = log10(realpow(10, maxlogverdadceil)/2);
    if maxlogverdadceil2<=maxlogverdad
      valsverdad(end+1) = maxlogverdadceil2;
    end
  end
  vals = log10(realpow(10, valsverdad)+1);
  ticks = vals;%ceil(min(nlog(:))+1):floor(max(nlog(:))+1);
  ticknums = realpow(10, valsverdad);%(10.^ticks);
else
  ticks = linspace(0, max(nlog(:)), nlabs);
  ticknums = (10.^ticks)-1;
end
if (nargin>7) && not(dofig)
  h = gcf;
else
  h = figure;
end  
if exist('colmap', 'var') && not(isempty(colmap))
  set(h, 'Colormap', colmap);
end
imagesc(c{1}, c{2}, nlog');
colorbar('YTick', ticks, 'YTickLabel', ticknums);
xlabel(nom1);
ylabel(nom2);
axis xy;
dcm = datacursormode(h);
set(dcm, 'Updatefcn', dtcbFUN(n', c{1}, valran1, c{2}, valran2));

function fun = dtcbFUN(values, ranx, valranx, rany, valrany)

fun = @(obj,event_obj) dtcb(values, ranx, rany, valranx,valrany,obj, event_obj);

function output_txt = dtcb(values, ranx, rany, valx,valy,obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
X = find(abs(ranx-pos(1))<1e-10);
Y = find(abs(rany-pos(2))<1e-10);
val=values(Y,X);
output_txt = {['X: ',num2str(valx(X),4)],...
    ['Y: ',num2str(valy(Y),4)], ...
    ['amount: ' num2str(val)]};
% output_txt = {['X: ',num2str(pos(1),4)],...
%     ['Y: ',num2str(pos(2),4)], ...
%     ['amount: ' num2str(val)]};

% % If there is a Z-coordinate in the position, display it as well
% if length(pos) > 2
%     output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
% end
