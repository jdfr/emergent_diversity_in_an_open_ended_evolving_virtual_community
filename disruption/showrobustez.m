function nms = showrobustez(data, whattoshow, noformatting)
% t2= [0.0 0.145812678684893, 0.188208509851167, 0.123415714463229, 0.0214214082253076, 0.0161128404336144, 0.0246690310408911, 0.0197555782257453, 0.00796051518195875;...
%      0.0 0.179127990456264, 0.234369164399853, 0.22864023397412, 0.0177764453469618, 0.01240866825042, 0.0178576475214396, 0.0128990204957026, 0.00565520686686067];    
% t3= [0.0 0.113801684343849,  0.159332634229786,  0.314925804465657, 0.115337244118034, 0.0400577186740039, 0.12503023159731, 0.0743241521376008, 0.0552444329866978;...
%      0.0 0.383986615901648, 0.541727642105433, 0.537207844662425, 0.0435899903944593, 0.0410427345678538, 0.047914661813722, 0.0432504164305447, 0.0290759094700387];
% % opFuncs = {@op_alteracionG, ...
% %            @op_eliminacionG, ...
% %            @op_insercionG, ...
% %            @op_dupaleatoriaG, ...
% %            @op_dupnivelG, ...
% %            @op_duplAleatoriaEdit, ...
% %            @op_duplNivelEdit, ...
% %            @op_duplSecuenciaEdit};

if not(exist('noformatting', 'var'))
  noformatting = false;
end

%figure;
strargs= {'fontname', 'times', 'fontsize', 18};
txtargs=[strargs, {'BackgroundColor', 'w','Margin', eps}];%,'EdgeColor','k'}];

if not(noformatting)
  a = axes; %#ok<NASGU>
  hold on;
end
%plot([0 0.4], [0 0.4], 'k');
f  =12;
shiftx = -0.006;
shifty = 0.017;
idxs = whattoshow.idxs;
mx = -inf;
nms = {}; nk=1;
hs = zeros(size(idxs));
for k=1:numel(idxs)
  names = regexp(data(idxs(k)).plotLabels, '([^|]+)', 'tokens');
  points = whattoshow.points{k};
  if isempty(points)
    points = 1:numel(names);
  end
  p1 = data(idxs(k)).disrupMean(points);
  p2 = data(idxs(k)).disrupMeanR(points);
  mx = max(mx, max([p1(:); p2(:)]));
  hs(k)=line(p1,p2,'LineStyle', 'none', 'Marker', whattoshow.marker{k}, 'MarkerSize', whattoshow.mksiz(k), 'MarkerEdgeColor', whattoshow.mkcole{k}, 'MarkerFaceColor', whattoshow.mkcole{k});%whattoshow.mkcole{k});
  %'MarkerFaceColor', whattoshow.mkcolf{k}, 'MarkerEdgeColor', whattoshow.mkcole{k});
  for kk=1:numel(points)
    text(data(idxs(k)).disrupMean(points(kk))+shiftx,data(idxs(k)).disrupMeanR(points(kk))+shifty,names{points(kk)},txtargs{:});
    nms{nk} = names{points(kk)}; nk = nk+1;
  end
end
if noformatting
  return
end
lm = whattoshow.limit;
if isempty(lm)
  lm = ceil(mx*10)/10;
end
lm = max(lm);
useLog = false;
if useLog
  inff = 0.009;
  logargs={'XScale', 'log', 'YScale', 'log'};
else
  logargs={};
  inff = 0;
end
if isempty(whattoshow.ticks)
  targs = {};
else
  targs = {'XTick', whattoshow.ticks{k}, 'YTick', whattoshow.ticks{k}};
end
line([inff lm], [inff lm], 'Color', 'k', 'LineWidth', 1);
if numel(idxs)>1
  legend(hs, whattoshow.legend, 'location', whattoshow.legendpos{k});
end
xlabel('mean disruption on evolved genomes', strargs{:});
ylabel('mean disruption on random genomes', strargs{:});
axis equal;
set(gca, 'XLim', [inff lm], 'YLim', [inff lm], targs{:}, logargs{:}, strargs{:});
%set(gca, 'XLim', [0 lm], 'YLim', [0 lm]);


