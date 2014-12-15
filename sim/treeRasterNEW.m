function [raster, offset, counts, dimension, isLeaf, maxiy, leafs]=treeRasterNEW(axiom, mode)
% axiom      cadena a representar
% pintar     representar graficamente si pintar=true
% tiempo     tiempo (s) para visualizar la grafica
%
% alto       alto del bounding box
% ancho      ancho del bounding box
% ramas      numero de ramrasterLeaveValas

%chanceFun = @()true;
useTransparency  = false;
if isstruct(mode)
  paintRandom = isfield(mode, 'randomize') && logical(mode.randomize);
  if isfield(mode, 'angle')
    alpha     = mode.angle;
  else
    alpha     = -5;
  end
  notjustPaint = not(isfield(mode, 'justPaint') && mode.justPaint);
  mode        = mode.shadowing;
else
  notjustPaint = true;
  paintRandom = false;
  % angle: '+' = rotate anticlockwise
  %        '-' = rotate clockwise
  alpha       = -5; % degrees (angle)
end
switch mode
    case 'canonical'
    case 'probabilistic'
        error('Mode probabilistic not implemented in treeRaster');
        %probability      = 0.8;
        %chanceFun        = @()rand<=probability;
    case 'transparentBranches'
        useTransparency  = true;
    otherwise
        error('Mode %s not understood!!!', any2str(mode));
end

% lengths of the lines F and G
% length_F = 1;
length_G = 10;

% initialization
xT = 0;
yT = 0;
aT = pi/2;
da = alpha/180*pi; % convert deg to rad
stkPtr = 1;
stackInc=1000;
stack = zeros(stackInc,3);

nramas = 0;

rasterInc = 1000;
% if(useTransparency)
    %ALWAYS USE THIS REPRESENTATION, SINCE WE ALWAYS WANT TO TELL LEAVES
    %APART OF BRANCHES
    rasterEmptyFunc = @(a,b)zeros(a,b,'int8');
    rasterEmptyVal          = int8(0);
    rasterFillVal           = int8(1);
    rasterLeaveVal          = int8(2); % Leaves can not be replaced by fill val, except those in the same level of inmediatly upper
% else
%     rasterEmptyFunc = @false;
%     rasterEmptyVal = false;
%     rasterFillVal = true;
% end
raster = rasterEmptyFunc(rasterInc, rasterInc);
npx = numel(raster);
yOffset = 1;
xOffset = round(rasterInc/2);

for i = 1:length(axiom)
    cmdT = axiom(i);

    switch cmdT
        % It is possible to add multiple cases here
        % in order to expand the program.
        case 'G'
            if paintRandom
                rndL  = length_G*(1.25-0.5*rand);
                newxT = xT + rndL*cos(aT);
                newyT = yT + rndL*sin(aT);
            else
                newxT = xT + length_G*cos(aT);
                newyT = yT + length_G*sin(aT);
            end

            minY = min(yT, newyT);
            if(minY + yOffset < 1)
                incSteps = ceil((1-minY-yOffset)/rasterInc);
                raster = [rasterEmptyFunc(incSteps*rasterInc, size(raster, 2)); raster]; %#ok<AGROW>
                npx = numel(raster);
                yOffset = yOffset + incSteps*rasterInc;
            end
            maxY = max(yT, newyT);
            if(maxY + yOffset > size(raster, 1))
                incSteps = ceil((maxY + yOffset - size(raster, 1))/rasterInc);
                %raster(end+1:end+(incSteps*rasterInc), :) = rasterEmptyVal;
                raster(end+(incSteps*rasterInc), 1) = rasterEmptyVal;
                npx = numel(raster);
            end

            minX = min(xT, newxT);
            if(minX + xOffset < 1)
                incSteps = ceil((1-minX-xOffset)/rasterInc);
                raster = [rasterEmptyFunc(size(raster, 1), incSteps*rasterInc), raster]; %#ok<AGROW>
                npx = numel(raster);
                xOffset = xOffset + incSteps*rasterInc;
            end
            maxX = max(xT, newxT);
            if(maxX + xOffset > size(raster, 2))
                incSteps = ceil((maxX + xOffset - size(raster, 2))/rasterInc);
                %raster(:, end+1:end+(incSteps*rasterInc)) = rasterEmptyVal;
                raster(1, end+(incSteps*rasterInc)) = rasterEmptyVal;
                npx = numel(raster);
            end

            %raster = rasterLine(raster, round(xT+xOffset), round(newxT+xOffset), round(yT+yOffset), round(newyT+yOffset), true);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Matlab seems to make a copy of the raster when calling to
            % rasterLine: so let's embed the code here
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            x1 = round(xT+xOffset);
            x2 = round(newxT+xOffset);
            y1 = round(yT+yOffset);
            y2 = round(newyT+yOffset);

            dx = x2 - x1;
            dy = y2 - y1;

            if(dx || dy)
                if(abs(dx) >= abs(dy))
                    if(dx>0)
                        %indexes = sub2ind(size(raster), floor(linspace(y1, y2, dx+1)), x1:x2);
                        indexes = floor(linspace(y1, y2, dx+1)) + ((x1-1):(x2-1))*size(raster,1);
                    else
                        %indexes = sub2ind(size(raster), floor(linspace(y1, y2, 1-dx)), x1:-1:x2);
                        indexes = floor(linspace(y1, y2, 1-dx)) + ((x1-1):-1:(x2-1))*size(raster,1);
                    end
                else
                    if(dy>0)
                        %indexes = sub2ind(size(raster), y1:y2, floor(linspace(x1, x2, dy+1)));
                        indexes = (y1:y2) + (floor(linspace(x1-1, x2-1, dy+1)))*size(raster,1);
                    else
                        %indexes = sub2ind(size(raster), y1:-1:y2, floor(linspace(x1, x2, 1-dy)));
                        indexes = (y1:-1:y2) + (floor(linspace(x1-1, x2-1, 1-dy)))*size(raster,1);
                    end
                end
            else
                indexes = sub2ind(size(raster),y1,x1);
            end
            
            if ~all(raster(indexes))
                nramas = nramas+1;
                %raster(y2,x2) = rasterLeaveVal;
            end
            raster(indexes(raster(indexes)~=rasterLeaveVal)) = rasterFillVal;
            raster(y2,x2) = rasterLeaveVal;
            
            %[ys xs] = find(raster); minX = min(xs); maxX = max(xs); minY = min(ys); maxY = max(ys); raster2 = raster(minY:maxY, minX:maxX); imshow(raster2+1, [1 1 1; 1 0 0; 0 1 0], 'InitialMagnification', 'fit'); set(gca, 'YDir', 'normal');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            xT = newxT;
            yT = newyT;
        case '+' % rotate anticlockwise
            if paintRandom
                aT = aT + da + (rand-0.5)*0.04;
            else
                aT = aT + da;
            end
        case '-' % rotate clockwise
            if paintRandom
                aT = aT - da + (rand-0.5)*0.04;
            else
                aT = aT - da;
            end
        case '[' % save current position
            stack(stkPtr,:) = [xT yT aT];
            stkPtr = stkPtr +1 ;
            if stkPtr>size(stack,1)
                stack(end+stackInc,1) = 0;
            end
        case ']' % return to former position (last save)
            stkPtr = stkPtr -1 ;
            xT = stack(stkPtr,1);
            yT = stack(stkPtr,2);
            aT = stack(stkPtr,3);
        otherwise
            error('Illegal character in the genome!!!!!');
    end

end

[ys xs values] = find(raster);

if isempty(values)
  raster = {zeros(0,1, 'int32'), zeros(0,1, 'int32')};
  offset = [1 1];
  counts = [0 0 0];
  dimension = [0 0];
  isLeaf = false(1,0);
  maxiy = zeros(1,0, 'int32');
  return
end

minX = min(xs);
maxX = max(xs);
minY = min(ys);
maxY = max(ys);

if any(size(raster)~=[maxY-minY+1, maxX-minX+1])
  raster = raster(minY:maxY, minX:maxX);
  [ys xs values] = find(raster);
  %this happens if the tree is just an horizontal branch
  if size(raster,1)==1
    %to return always column vectors
    ys = ys';
    xs = xs';
    values=values';
  end
end

[alto ancho] = size(raster);
clear raster;
dimension = [alto ancho];

xOffset = xOffset - minX + 1;
yOffset = yOffset - minY + 1;

offset = [yOffset xOffset];

if not(notjustPaint)
  raster = {ys xs};% values};
  return
end

% [y,x] = ind2sub(size(raster), find(raster==rasterLeaveVal));
% leaf = [y,x];

onlyLeaves      = values==rasterLeaveVal;

if nargout>6
  %also return leafs
  leafs = {ys(onlyLeaves), xs(onlyLeaves)};
end

nleafs          = sum(onlyLeaves);
nPixelsNeg      = sum(ys<yOffset);

if useTransparency
  tempRaster      = sparse(ys(onlyLeaves), xs(onlyLeaves), double(ys(onlyLeaves)), dimension(1), dimension(2));
else
  nleafs          = sum(values==rasterLeaveVal);
  tempRaster      = sparse(ys, xs, ys, dimension(1), dimension(2));
end

maxiy = int32(full(max(tempRaster,[],1)));

% minval = int32(-Inf);
% 
% tempRaster = repmat((int32(1):alto)', 1, ancho);
% if(useTransparency)
%     tempRaster(raster~=rasterLeaveVal) = minval;
% else
%     tempRaster(~raster) = minval;
% end
% maxiy = max(tempRaster);
% maxiy(maxiy == minval) = 0;
% nleafs = sum(values==rasterLeaveVal);

counts  = [nramas nleafs nPixelsNeg];

if useTransparency
  isLeaf = true(size(maxiy));
else
  raster = sparse(ys(onlyLeaves),xs(onlyLeaves),2, dimension(1), dimension(2));
  isLeaf = int8(full(raster(sub2ind(size(raster), maxiy, int32(1):ancho))))==rasterLeaveVal;
  clear raster;
end
clear onlyLeaves;
ys=int32(ys);
xs=int32(xs);
raster = {ys xs};% values};
maxiy = maxiy-yOffset;