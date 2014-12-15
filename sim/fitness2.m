function [P colors params]=fitness2(P,params,Generation,parentdata,drawToScreen,randN,addparams)
%drawToScreen and randN are params to be used when redoing the landscape
if ~exist('drawToScreen', 'var')
  drawToScreen = false;
end
if ~exist('addparams', 'var')
  addparams = struct('overridePlotting', false, 'basedir', '', 'namefile', '');
end


maxiy       = P.maxiy;
isLeaf      = P.isLeaf;
offset      = P.offset;
rasters     = P.raster;
dimensions  = P.dimensions;
% todos los arboles estan pintados de 1 en adelante
sR = numel(maxiy);
N = params.N;

%this step is to ensure that a too wide tree will not crash the algorithm
N = max(N, max(dimensions(:,2))+2);

backgroundColor = [1 1 1];
floorColor = [0 0 0];
colors = rand(sR, 3)*0.8;
%colors(colors>0.8) = 0.8;

rasterYUpside   = dimensions(:,1)-offset(:,1)+1;
rasterYDownside = offset(:,1)-1;

if ~isempty(params.heightT)
  heights          = cellfun(@safeMax, maxiy);
  areTooTall       = heights>params.heightT;
end
tooTallArePlaced   = isempty(params.heightT) || params.tooTallArePlaced;

tooManyNeg         = params.negYThreshold;
negY               = P.counts(:,3);
haveNotTooManyNeg  = negY<=tooManyNeg;

plotting           = params.plotting && (addparams.overridePlotting || (mod(Generation, params.plottingFactor)==0));

if params.randomPlacementForInitialGeneration && (Generation==params.initialGeneration)
   treePlacement          = 'random';
else
   treePlacement          = params.treePlacement;
end

if ~exist('randN', 'var')
   if ~isempty(parentdata)
     parentwidth            = parentdata.width;
   end
   width                    = dimensions(:,2);
   offsetAdd                = params.treePlacementOffset;
   switch treePlacement
     case 'random' %quasi purely random placement (there are edge effects)
       if isinf(N)
         error('purely random placement is not possible for unbounded environments!!!');
       end
       randN                = ceil(rand(sR,1).*(N-width));
     case 'uniformDist' %historical placement with uniform distribution (overflow at the edges is reflected into the landscape)
       anchoPlacement       = params.treePlacementParam*parentwidth+offsetAdd;
       uniformRandom        = anchoPlacement.*(rand(size(P.randN))-1/2);
       randN                = round(P.randN+uniformRandom);
       randN                = rebotarEnBordes(randN, N, width);
       clear anchoPlacement uniformRandom;
     case 'uniformDistByRoot' %same as before, but using the root as reference instead of the bounding box
       anchoPlacement       = params.treePlacementParam*parentwidth+offsetAdd;
       uniformRandom        = anchoPlacement.*(rand(size(P.randN))-1/2);
       %now, randN does not specify the leftmost position of the tree, but
       %the root's position
       randN                = round(P.randN+uniformRandom);
       %we must use a customized, inlined version of the rebotarEnBordes
       %function
              goon                 = ~isinf(N); %do not handle this if the environment is unbounded
              while goon %this might be done several times, if the offset is extremely big
               goon               = false;
               tooLeft            = (randN-offset(:,2)+1)<1;
               if any(tooLeft)
                 goon             = true;
                 newLeftmost      = 2-(randN(tooLeft)-offset(tooLeft,2)+1);
                 randN(tooLeft)   = newLeftmost+offset(tooLeft,2)-1;
               end
               tooRight           = (randN+width-offset(:,2))>N;
               if any(tooRight)
                 goon             = true;
                 newLeftmost      = (N-width(tooRight)+1)-(randN(tooRight)-offset(tooRight,2)+1-(N-width(tooRight)+1));
                 randN(tooRight)  = newLeftmost+offset(tooRight,2)-1;
               end
              end
       clear anchoPlacement uniformRandom goon tooLeft newLeftmost tooRight;
     case 'gaussianDist' %historical placement with gaussian distribution (overflow at the edges is reflected into the landscape)
       sigma                = params.treePlacementParam*parentwidth+offsetAdd;
       normalRandom         = -realsqrt(2)*erfcinv(2*rand(size(P.randN))).*sigma;
       randN                = round(P.randN + normalRandom);
       randN                = rebotarEnBordes(randN, N, width);
       clear normalRandom sigma;
     case 'covering' %historical placement strictly inside the plant width
       randN                = P.randN+ceil(rand(size(P.randN)).*parentwidth)-offset(:,2)+1;
       randN(randN<1)       = 1;
       tooRight             = (randN+width-1)>N;
       randN(tooRight)      = N-width(tooRight)+1;
       clear tooRight;
     otherwise      
       error('params.treePlacement not understood: %s!!!!!', any2str(treePlacement));
   end
   clear width offsetAdd parentwidth;
end

leftmost = calculateLeftmost(treePlacement, randN, offset);

if isinf(N)
  %unbounded environment: the environment must be given the width to
  %encompass all trees. First, the leftmost leftmost must be adjusted to 2
  leftmost         = leftmost-min(leftmost)+2;
  %Second, calculate the rightmost position, that will mark the
  %environment's width
  environmentWidth = max(leftmost+dimensions(:,2));
else
  environmentWidth = N;
end

count = ones(1,environmentWidth); % cuando dos hojas solapan decidir quien se queda con la luz aleatoriamente, teniendo en cuenta cuantas hojas han pasado por ese x hasta ese momento 

label = zeros(1,environmentWidth);
labelIsLeaf = false(size(label));
% labelHeight = repmat(int32(-Inf), size(label));
labelHeight = zeros(size(label), 'int32'); %no tiene en cuenta las ramas que estan subterraneas
if params.correc.zeroHeightOK
  labelHeight = labelHeight-1;
end

if(plotting)
%     h = figure('Visible', 'off');
%     axis equal;
%     hold on;

    colors = [backgroundColor; floorColor; colors]; 
    
    toUse = haveNotTooManyNeg;
    if ~(isempty(params.heightT) || tooTallArePlaced)
      toUse = toUse & (~areTooTall);
    end
    maxY = max(rasterYUpside(toUse));
    minY = max(rasterYDownside(toUse));
    if isempty(maxY)
      maxY = 0;
      minY = 0;
    end
    clear toUse;
    
    height = maxY + minY + 11; % 10 points of margin (5 top, 5 bottom)
    yOffset = minY + 6;
    %BEWARE: indexed images can be stored in PNG format only with a bit
    %depth of 8, not 16. So, if we have more than 250 trees in the
    %population, we will run into trouble. The solution might be to use a
    %truecolor image, declaring raster as a 'uint8' cubic matrix, and
    %translating [r,g,b] in the range 0..1 to the range 0..255
    raster = zeros(height, environmentWidth, 'uint8');
    raster(yOffset, :) = 1; % floor
end

for i = 1 : sR
    if (~isempty(maxiy{i})) && (tooTallArePlaced || (~areTooTall(i))) && haveNotTooManyNeg(i)
        if(plotting)
          ys     = rasters{i}{1};
          xs     = rasters{i}{2};
          displacement = yOffset-offset(i,1)+1;
          if (displacement>(yOffset+1)) || (displacement<1)
            error('This is plainly absurd: %d %d!!!!', yOffset, offset(i,1));
          end
          %add 1 to i instead of 2, since the raster is uint16, and the
          %base color is 0, not 1
          if (i+1)<=250
            raster(sub2ind(size(raster), ys+displacement, xs+(leftmost(i)-1))) = i+1;
          else
            %This is not very elegant, but the alternative is a cubic
            %matrix, not a good idea for huge matrices
            %assign other colors to other trees
            raster(sub2ind(size(raster), ys+displacement, xs+(leftmost(i)-1))) = mod(i, 250)+3;
          end
        end
        
        segment                     = leftmost(i):(leftmost(i)+dimensions(i,2)-1);
        %clever: we are not counting underground leaves, since they have
        %negative maxiy values, and labelHeight values are 0 at least
        putLabel                    =         maxiy{i} >  labelHeight(segment);
        putLabelS                   = segment(putLabel);
        %conflicting positions: there already is another tree with the same
        %height at that position
        coincidentLabel             = find((maxiy{i} == labelHeight(segment)) & (maxiy{i}>=0));
%         if ~isempty(coincidentLabel) && length(maxiy{i}) > 1
%             disp('hola')
%         end
%         if ~isempty(find(isLeaf{i}==0))
%             disp('jo')
%         end
        coincidentLabelS            = segment(coincidentLabel);%#ok<FNDSB> % & (maxiy{i}>-1));
        
        %update positions where I am the heightest
        labelHeight(putLabelS)      = maxiy{i}(putLabel);
        label(putLabelS)            = i;
        labelIsLeaf(putLabelS)      = isLeaf{i}(putLabel);
        
        %update positions where my height is equal to other trees, if I
        %win the tournament
        count(coincidentLabelS)     = count(coincidentLabelS)+1;
        seizeProbability            = (count(coincidentLabelS)-1)./count(coincidentLabelS);
        randomChoice                = rand(size(coincidentLabelS)) > seizeProbability;
        alsoPutLabelS               = coincidentLabelS(randomChoice);
        %unnecessary: %labelHeight(alsoPutLabelS)  = maxiy{i}(alsoPutLabel);
        label(alsoPutLabelS)        = i;
        labelIsLeaf(alsoPutLabelS)  = isLeaf{i}(coincidentLabel(randomChoice));
        
        clear segment putLabel putLabelS coincidentLabel coincidentLabelS randomChoice alsoPutLabelS
    end

end

clear leftmost;

if(plotting) && (~isempty(raster))
    raster(yOffset,label>-1) = label(label>-1)+1; raster(yOffset,labelIsLeaf) = 0;
    if drawToScreen
      figure; imshow(full(raster(end:-1:1,:)), colors);%, 'InitialMagnification', 'fit');
    else
      if isempty(addparams.basedir)
        basedir = params.nomdir;
      else
        basedir = addparams.basedir;
        if basedir(end)~=filesep
          basedir(end+1) = filesep;
        end
      end
      if isempty(addparams.namefile)
        namefile = ['entornoRaster' num2str(Generation, '%03g') '.png'];
      else
        namefile = addparams.namefile;
      end
      imwrite(full(raster(end:-1:1,:)), colors(1:min(size(colors,1),252),:), [basedir namefile], 'png');
    end
    clear raster;
    colors = colors(3:end,:);
end

valid             = label>0;
label             = label(valid);
labelIsLeaf       = labelIsLeaf(valid);
contributions     = sparse(1:numel(label), label, labelIsLeaf, numel(label), numel(P.fitness));
P.fitness(1:end)  = sum(contributions);
clear valid label labelIsLeaf labelHeight contributions;
switch params.divideFitness
  case 'byLeafs'
    % dividiendo por el numero de hojas del arbol
    nleafs   = P.counts(:,2);
    do       = nleafs>0; % |G|= 10
    P.fitness(do)  = (P.fitness(do)-negY(do))./((nleafs(do)*params.factorG).^params.alphaG);
    P.fitness(~do) = 0;
  case 'byBranches'
    % dividiendo por el numero de ramas distinguibles del arbol *|G|
    nramas         = P.counts(:,1);
    do             = nramas>0; % |G|= 10
    P.fitness(do)  = (P.fitness(do)-negY(do))./((nramas(do)*params.factorG).^params.alphaG);
    P.fitness(~do) = 0;
end
clear negY;

if ~isempty(params.heightT)
  P.fitness(areTooTall) = 0;%min(P.fitness);
end

P.randN(1:end) = randN(1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = safeMax(array)
if isempty(array);  m=int32(0);  else  m=max(array); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function randN = rebotarEnBordes(randN, N, width)
%if any value in randN is outside the range 1..N, place it inside the
%range, by mirroring it in the limits 1 and N
goon                = ~isinf(N); %do not handle this if the environment is unbounded
while goon %this might be done several times, if the offset is extremely big
 goon               = false; 
 tooLeft            = randN<1;
 if any(tooLeft)
   goon             = true;
   randN(tooLeft)   = 2-randN(tooLeft);
 end
 tooRight           = (randN+width-1)>N;
 if any(tooRight)
   goon             = true;
   randN(tooRight)  = (N-width(tooRight)+1)-(randN(tooRight)-(N-width(tooRight)+1));
 end
end

function leftmost = calculateLeftmost(treePlacement, randN, offset)
switch treePlacement
  case 'uniformDistByRoot'
    leftmost  = randN-offset(:,2)+1;
  otherwise
    leftmost  = randN;
end

