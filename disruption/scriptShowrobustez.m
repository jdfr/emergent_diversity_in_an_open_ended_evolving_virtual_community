function scriptShowrobustez(data)

nz = 1000;
%figure;
cols = repmat('b', nz, 1);
ix = 1;
%points = repmat({zeros(nz,2)}, 5, 1);
for k=1:nz
  d  = data;
  r  = randperm(numel(d(ix).pop))';
  rr = randperm(numel(d(ix).pop))';
  n  = 64;
  r  = r(1:n);
  rr = rr(1:n);
  d(ix).disrupMean  = mean(d(ix).disrup, 3);
  d(ix).disrupMeanR = mean(d(ix).disrupR, 3);
  d(ix).disrupMean  = mean(d(ix).disrupMean(:,r), 2);
  d(ix).disrupMeanR = mean(d(ix).disrupMeanR(:,rr), 2);
  for z=1:5;
    line(d(ix).disrupMean(z), d(ix).disrupMeanR(z), 'Marker', '.', 'MarkerFaceColor', cols(k,:), 'MarkerEdgeColor', cols(k,:));
    %points{z}(k,:) = [d(ix).disrupMean(z), d(ix).disrupMeanR(z)];
  end;
end;
% gr = 0:0.0005:0.4;
% [X Y] = meshgrid(gr, gr);
% points = vertcat(points{:});
% h = 0.003;
% mask = ((X<0.05)&(Y<0.05)) | ((X>0.025)&(X<0.25)&(Y>0.175)&(Y<0.3));
% x = X(mask);
% y = Y(mask);
% f = KDE2(points(:,1),points(:,2),h,x,y);
% F = zeros(size(X));
% F(mask) = f;
% q=F/max(F(:))/0.2; q(q>1) = 1;
% image('Clipping', 'off', 'AlphaData', q, 'CData', cat(3, zeros(size(q)), zeros(size(q)), ones(size(q))), 'XData', gr([1 end]), 'YData', gr([1 end]));
% axis(gr([1 end 1 end]));%[0 0.4 0 0.4]);
% axis square;
% axis xy;
gr = [0 0.3];
axis(gr([1 end 1 end]));
axis square;

if nargin<1 || isempty(data)
  d1 = load('T2disp_caguento.mat');
  d2 = load('T3disp.mat');

  data = [d1.data; d2.data];

  clear d1 d2;
end

whattoshow = struct(...
  'idxs', 1:6, ...
  'points', {{1:5, 1:5, 1:5}},...%{cell(1,6)}, ...
  'marker', {{'o' 'd' 's' '.' '.' '.'}}, ...
  'mkcolf', {{'none' 'none' 'none' 'none' 'none' 'none'}}, ...
  'mkcole', {repmat({'k'},1, 6)}, ...%{{'r' 'b' 'g' 'c' 'm' 'k'}}, ...
  'mksiz', [8 8 8], ...%[10 10 10 15 15 15], ...
  'limit', [0.31 0.31 0.31], ...%+0.1, ...
  'ticks', {repmat({0:0.1:0.4}, 3, 1)}, ...
  ...'legend', {{'very harsh T2', 'harsh T2', 'mild T2', 'very harsh T3', 'harsh T3', 'mild T3'}} ...
  'legend', {{'\alpha = 1     (very harsh)', 'harsh T2', '\alpha = 0.5  (mild)', 'very harsh T3', 'harsh T3', 'mild T3'}}, ...
  'legendpos', {repmat({'SouthEast'}, 1, 3)} ...
  );

% fn = fieldnames(whattoshow);
% for k=1:numel(fn)
%   whattoshow.(fn{k}) = whattoshow.(fn{k})(toshow);
% end
% 
% showrobustez(data, whattoshow, general);

%comps = {'veryharsh' [1 4]; 'harsh' [2 5]; 'mild' [3 6]};
comps = {'very harsh' [1 3]; 'harsh' [2 ]; 'mild' [3 ]; };
toskip = [false true true]; %[false true false];
noformatting = [false true true];
fn = fieldnames(whattoshow);
nms = cell(size(comps,1),1);
for z=1:size(comps,1)
  if toskip(z)
    continue;
  end
  w = whattoshow;
  for k=1:numel(fn)
    w.(fn{k}) = w.(fn{k})(comps{z,2});
  end

  nms{z} = showrobustez(data, w, noformatting);
%   title(sprintf('disruption analysis (%s)', comps{z,1}));
end
%title(sprintf('disruption analysis'));
nms = [nms{:}];
