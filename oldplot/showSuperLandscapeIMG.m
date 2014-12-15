function varargout = showSuperLandscapeIMG(img, width)
%show an image that is too long to be shown properly on windows

colors = [1 1 1; 0 0 0];% rand(255, 3)*0.8];
fig    = figure('Toolbar', 'figure');
shw    = imshow(img(:,1:min(width, size(img,2))), colors, 'InitialMagnification', 'fit');
if width<size(img,2)
  step   = width/(size(img,2)-width+1);
  if step<0; step=0.1; end
  slid   = uicontrol('Style', 'slider', 'Parent', fig, 'Units', 'normalized', 'Position', [0 0 1 0.03], 'SliderStep', [step/10 step]);
  set(slid, 'Callback', {@repaint, img, fig, slid, shw, width});
end
if nargout>0
  varargout{1} = img;
end

function repaint(obj, event, img, fig, slid, shw, width) %#ok<INUSL>
value = get(slid', 'Value');
value = ceil(value*(size(img,2)-width+1));
if value <= 0; value = 1; end
mx = (value+width-1);
set(slid, 'ToolTipString', sprintf('xinit: %010d, yinit: %010d', value, mx));
set(shw, 'CData', img(:,value:mx));
