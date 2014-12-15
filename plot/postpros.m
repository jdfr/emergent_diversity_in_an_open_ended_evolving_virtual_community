function postpros(h)
%a hack to make a figure prettier
figargs =  {'Color', 'w'};
axesargs = {'FontName','times','FontSize',14};
labargs = {'FontName','times','FontSize',14};
a = get(h, 'currentaxes');
set(h, figargs{:});
set(a, axesargs{:});
set(get(a, 'XLabel'), labargs{:});
set(get(a, 'YLabel'), labargs{:});
set(get(a, 'Title'),  labargs{:});
