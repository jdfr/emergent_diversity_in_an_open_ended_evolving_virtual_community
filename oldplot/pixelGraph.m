function A = pixelGraph(ramas,ran)

xmax = ran(1);
xmin = ran(2);
ymax = ran(3);
ymin = ran(4);
A={};
% xmin,xmax,ymin,ymax de entre todas las matrices ramas de todos los
% individuos de la poblacion
sR = size(ramas,2);
for i = 1 : sR
    h = figure;
% %     set(h,'color','k')
%     axis([xmin xmax ymin ymax])% axis([xmin xmax ymin ymax zmin zmax])
    % Units = 'pixels';
%     set(h,'paperUnits','points')
%     set(h,'paperPosition',[0 0 100 100]) % rect = [left, bottom, width, height]
%     set(h,'Position',[0 0 200 200])

% set(h,'Position', [0 0 200 200], 'paperUnits', 'points', 'paperPosition',[0 0 200 200]);
% set(h,'paperUnits', 'points', 'paperPosition',[0 0 200 200], 'Position', [0 0 200 200]);

drawtree(ramas{i},ran,1)
frame = imcapture(h, 'img', 0);
frame = frame(:,:,1);
A{i} = ~logical(frame);

end


function frame = extractFrameFromFigure(fig)

% The renderer requires DPI units.  Get the necessary conversion factor.
pixelsperinch = get(0,'screenpixelsperinch');

% The easiest way to obtain the conversion is to switch to pixel units.
% This provides robustness to whatever units the user was using previously.
oldUnits = get(fig,'units');
set(fig,'units','pixels');
set(fig, 'paperposition', pos./pixelsperinch);

% Now we obtain the figure size
pos = get(fig,'position');
oldPaperPosition = get(fig,'paperposition');
set(fig, 'paperposition', pos./pixelsperinch);

% Upgrade no renderer and painters renderer to OpenGL for the off-screen
% rendering, because that's what AVIFILE/ADDFRAME does.
renderer = get(fig,'renderer');
if strcmp(renderer,'painters') || strcmp(renderer,'None')
  renderer = 'opengl';
end

% Temporarily turn off warning in case opengl is not supported and
% hardcopy needs to use zbuffer
warnstate = warning('off', 'MATLAB:addframe:warningsTurnedOff');
noanimate('save',fig);
frame = hardcopy(fig, ['-d' renderer],['-r' num2str(round(pixelsperinch))]);
noanimate('restore',fig);
warning(warnstate);

% Restore figure state.  We do it in the opposite order that it was changed
% so dependent state elements work correctly (esp. paperposition's
% interplay with units).
set(fig, 'paperposition', oldPaperPosition);
set(fig, 'units',oldUnits); 

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% steep = abs(y1 - y0) > abs(x1 - x0);
% a = x0;
% b = y0;
% c = x1;
% d = y1;
% if steep
%     x0=b;
%     y0=a;
%     x1=d;
%     y1=c;
% end
% a = x0;
% b = y0;
% c = x1;
% d = y1;
% if x0 > x1
%     x0=c;
%     x1=a;
%     y0=d;
%     y1=b;
% end
% deltax = x1 - x0;
% deltay = abs(y1 - y0);
% error = deltax / 2;
% y = y0;
% if y0 < y1
%     ystep = 1;
% else
%     ystep = -1;
% end
% k = 1;
% for x = x0 : x1
%     if steep
%         axisx(k)= y;
%         axisy(k) = x;
% %         A(y,x)=1;
%         
%     else
%         axisx(k)= y;
%         axisy(k) = x;
% %         A(x,y)=1;
%         
%     end
%     k = k + 1;
%     error = error - deltay;
%     if error < 0
%         y = y + ystep;
%         error = error + deltax;
%     end
% end
% % A = A(end:-1:1,:)';


