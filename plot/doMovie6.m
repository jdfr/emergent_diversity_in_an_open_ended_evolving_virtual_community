function doMovie6(option)
%ffmpeg -y -i zfmildb1/rastermild_b1%05d.png -vcodec libx264 -level 41 -crf 30 -g 20 -r 20 -coder 1 zmildb1.mp4

closeExaminationTransitionFrames = 0;
switch option
  case 1
generations      = 0:10:500;
generationdetail = 70;
images           = arrayfun(@(x)sprintf('..\\..\\4tipos\\veryharsh\\entornoRaster%03d.png', x), generations, 'uniformoutput', false);
moviename        = 'fveryharsh\\rastervharsh';
stillFrames      = 20;
maxpositions     = 15;
stillFrames2     = 4;
transitionFrames = 10;
blend2Frames     = 10;
sizesupersheet   = [600 800];
sizescreen       = [sizesupersheet(1)*2/3,           sizesupersheet(2)];
sizesheet        = [sizesupersheet(1)-sizescreen(1), sizesupersheet(2)];
coverfile        = '..\cover2.png';
  case 0
generations      = 00:10:550;
generationdetail = 70;
images           = arrayfun(@(x)sprintf('..\\..\\4tipos\\harshred\\entornoRaster%03d.png', x), generations, 'uniformoutput', false);
moviename        = 'fharshred2\\rasterharsh';
stillFrames      = 20;
maxpositions     = 15;
stillFrames2     = 4;
transitionFrames = 10;
blend2Frames     = 10;
sizesupersheet   = [600 800];
sizescreen       = [sizesupersheet(1)*2/3,           sizesupersheet(2)];
sizesheet        = [sizesupersheet(1)-sizescreen(1), sizesupersheet(2)];
coverfile        = '..\cover2.png';
  case 2
generations      = 0:10:500;
generationdetail = 70;
images           = arrayfun(@(x)sprintf('harsh2\\entornoRaster%03d.png', x), generations, 'uniformoutput', false);
moviename        = 'fharsh2\\rasterharsh';
stillFrames      = 20;
maxpositions     = 15;
stillFrames2     = 4;
transitionFrames = 10;
blend2Frames     = 10;
sizesupersheet   = [600 800];
sizescreen       = [sizesupersheet(1)*2/3,           sizesupersheet(2)];
sizesheet        = [sizesupersheet(1)-sizescreen(1), sizesupersheet(2)];
coverfile        = '..\cover2.png';
  case 3
generations      = 0:10:550;
generationdetail = 7000;
images           = arrayfun(@(x)sprintf('..\\..\\4tipos\\harshred\\entornoRaster%03d.png', x), generations, 'uniformoutput', false);
moviename        = 'zfharshred\\rasterharsh';
stillFrames      = 20;
maxpositions     = 0;
stillFrames2     = 4;
transitionFrames = 10;
blend2Frames     = 10;
sizesupersheet   = [256 1024];
sizescreen       = [0 0];[sizesupersheet(1)*2/3,           sizesupersheet(2)];
sizesheet        = sizesupersheet;[sizesupersheet(1)-sizescreen(1), sizesupersheet(2)];
coverfile        = '..\cover3.png';
  case 4
generations      = 0:10:500;
generationdetail = 7000;
images           = arrayfun(@(x)sprintf('..\\..\\4tipos\\veryharsh\\entornoRaster%03d.png', x), generations, 'uniformoutput', false);
moviename        = 'zfveryharsh\\rastervharsh';
stillFrames      = 20;
maxpositions     = 0;
stillFrames2     = 4;
transitionFrames = 10;
blend2Frames     = 10;
sizesupersheet   = [256 1024];
sizescreen       = [0 0];[sizesupersheet(1)*2/3,           sizesupersheet(2)];
sizesheet        = sizesupersheet;[sizesupersheet(1)-sizescreen(1), sizesupersheet(2)];
coverfile        = '..\cover3.png';
  case 5
generations      = 0:10:230;
generationdetail = 7000;
images           = arrayfun(@(x)sprintf('..\\..\\4tipos\\mild\\entornoRaster%03d.png', x), generations, 'uniformoutput', false);
moviename        = 'zfmild\\rastermild';
stillFrames      = 40;
maxpositions     = 0;
stillFrames2     = 4;
transitionFrames = 20;
blend2Frames     = 20;
sizesupersheet   = [600 800];
sizescreen       = [0 0];[sizesupersheet(1)*2/3,           sizesupersheet(2)];
sizesheet        = sizesupersheet;[sizesupersheet(1)-sizescreen(1), sizesupersheet(2)];
coverfile        = '..\cover2.png';
  case 6
generations      = 0:10:550;
generationdetail = 0;
images           = arrayfun(@(x)sprintf('..\\..\\4tipos\\harshred\\entornoRaster%03d.png', x), generations, 'uniformoutput', false);
moviename        = 'zfharshred\\rasterharsh';
maxpositions     = {...
  {140 [1 3050]}; ...
  {150 [3270]}; ...
  {160 [170 3300]}; ...
  {170 [2800 3580]}; ...
  {180 [3715]}; ...
  {190 [1800 3005]}; ...
  {200 [785 2000 4490]}; ...
  {210 [980 1972 5180]}; ...
  {220 [624]}; ...
  {230 [2520 4920]}; ...
  {240 [4810]};
  {250 [5105]}; ...
  {260 [2465]}; ...
  {270 [2920 6055]}; ...
  {280 [1 2890 5670]}; ...
  {290 [2980 7130]}; ...
  {300 [1 2860 5575]}; ...
  {310 [6015 8435]}; ...
  {320 [4910 8485]}; ...
  {330 [7345 8580]}; ...
  {340 [4685 7425 8680]}; ...
  {350 [1110 2200 9010]}; ...
  {360 [1150 4865]}; ...
  {370 [1375]}; ...
  {380 [1 1215 2545]}; ...
  {390 [1 4920]}; ...
  {400 [1 4550]}; ...
  {410 [1 3395 8045]}; ...
  {420 [4445 5945]}; ...
  {430 [4545 6330]}; ...
  {440 [5275 ]}; ...
  {450 [455 6325]}; ...
  {460 [655 4810 6090]}; ...
  {470 [5425]}; ...
  {480 [5645]}; ...
  {490 [7105]}; ...
  {500 [5095 11745]}; ...
  {510 [5225 7415 12980]}; ...
  {520 [1 5505 7545]}; ...
  {530 [220 2390 5775]}; ...
  {540 [6210 7755 11790]}; ...
  {550 [1 3195 9145]}; ...
  };
stillFrames      = 20;
stillFrames2     = 15;4;
transitionFrames = 10;
blend2Frames     = 10;
closeExaminationTransitionFrames = 5;
sizesupersheet   = [600 800];
sizescreen       = [sizesupersheet(1)*2/3,           sizesupersheet(2)];
sizesheet        = [sizesupersheet(1)-sizescreen(1), sizesupersheet(2)];
coverfile        = '..\cover2.png';
  case 7
generations      = 0:10:550;
generationdetail = 0;
images           = arrayfun(@(x)sprintf('..\\..\\4tipos\\harshred\\entornoRaster%03d.png', x), generations, 'uniformoutput', false);
moviename        = 'zfharshred2\\rasterharsh';
maxpositions     = {...
  {170 [3580]}; ...
  {250 [5105]}; ...
  {330 [8580]}; ...
  {540 [6210]}; ...
  };
stillFrames      = 20;
stillFrames2     = 20*15;15;4;
transitionFrames = 10;
blend2Frames     = 10;
closeExaminationTransitionFrames = 10;
sizesupersheet   = [600 800];
sizescreen       = [sizesupersheet(1)*2/3,           sizesupersheet(2)];
sizesheet        = [sizesupersheet(1)-sizescreen(1), sizesupersheet(2)];
coverfile        = '..\cover2.png';
  case 8
generations      = 0:10:90;
generationdetail = 7000;
images           = arrayfun(@(x)sprintf('..\\..\\4tipos\\mild_b1\\entornoRaster%03d.png', x), generations, 'uniformoutput', false);
moviename        = 'zfmildb1\\rastermild_b1';
stillFrames      = 40;
maxpositions     = 0;
stillFrames2     = 4;
transitionFrames = 20;
blend2Frames     = 20;
sizesupersheet   = [600 800];
sizescreen       = [0 0];[sizesupersheet(1)*2/3,           sizesupersheet(2)];
sizesheet        = sizesupersheet;[sizesupersheet(1)-sizescreen(1), sizesupersheet(2)];
coverfile        = '..\cover2.png';
end
labelScale       = 0.1;
labelShift       = [5 5];
doMovie3_doit(...
  coverfile, ...
  moviename, ...
  sizesupersheet, ...
  sizesheet, ...
  sizescreen, ...
  images, ...
  generations, ...
  generationdetail, ...
  labelScale, ...
  labelShift, ...
  stillFrames, ...
  maxpositions, ...
  stillFrames2, ...
  transitionFrames, ...
  blend2Frames, ...
  closeExaminationTransitionFrames);
end

function doMovie3_doit(coverfile, moviename, sizesupersheet, sizesheet, sizescreen, images, ...
                       generations, generationdetail, labelScale, labelShift, ...
                       stillFrames, maxpositions, stillFrames2, transitionFrames, blend2Frames, closeExaminationTransitionFrames)


coverParams.imfile          = coverfile;
coverParams.scaleImage      = 1.8;
coverParams.backgroundColor = [255 255 255];
coverParams.numFramesStill  = 20;
coverParams.numFramesBlend  = 10;
coverParams.numFramesStill2 = 10;
                     
                     
resizeParams      = {'Antialiasing', true, 'Method', 'lanczos2'};

labels            = generateLabels(generations);
if ~isempty(labelScale)
  labels          = cellfun(@(x)imresize(x, 'Scale', labelScale, resizeParams{:}), labels, 'uniformoutput', false);
end

supersheet        = repmat(uint8(255), [sizesupersheet(1), sizesupersheet(2), 3]);

screen            = repmat(uint8(255), [sizescreen(1), sizescreen(2), 3]);
ratioscreen       = sizescreen(1)/sizescreen(2);
screenShift       = [1, 1];

sheet             = repmat(uint8(255), [sizesheet(1), sizesheet(2), 3]);
ratiosheet        = sizesheet(1)/sizesheet(2);
sheetShift        = [sizescreen(1)+1, 1];

% sizeall     = areaall(3:4)    - areaall(1:2)   + 1;
% sizedetail  = areadetail(3:4) - areadetail(1:2) + 1;
% ratioall    = sizeall(1)      / sizeall(2);
% ratiodetail = sizedetail(1)   / sizedetail(2);

labelAlphaColor = uint8([255 255 255]);
imageSkyColor   = uint8([250 255 255]);

[thisim map]      = imread(images{1});
if ~isempty(map)
  thisim          = myind2rgb(thisim, map);
end
thisim            = thisim(1:end-6,:,:);
thisim            = replaceSky(thisim, imageSkyColor);
sizethisim        = size(thisim);
maxx              = sizethisim(1);
maxy              = sizethisim(2);
if ratiosheet<maxx/maxy
  scalethis       = sizesheet(1)/maxx;
else
  scalethis       = sizesheet(2)/maxy;
end

colorthis = [0 0 0];
colorotro = [0 0 0]+0.7;
markthickness = 1;

[frame pos]       = putOnSheet(thisim, scalethis, resizeParams, sheet, sizesheet);
img               = frame((pos(1):(pos(1)+pos(3)-1)), (pos(2):(pos(2)+pos(4)-1)), :);
img               = putmark(img, markthickness, colorthis);
frame             = copyOn(frame, pos(1:2), img);

numceros = '5';

%numpixelstomove = 500;

numb = 0;
f = 0;

if iscell(maxpositions)
  gensCloser        = cellfun(@(x)x{1}, maxpositions);
  originalpositions = cellfun(@(x)x{2}, maxpositions, 'uniformoutput', false);
end


% try
  for k=2:numel(images)
    %copy initial image for this generation
%     img           = frame((pos(1):(pos(1)+pos(3)-1)), (pos(2):(pos(2)+pos(4)-1)), :);
%     img           = putmark(img, markthickness, colorthis);
%     frame         = copyOn(frame, pos(1:2), img);
    
    framel        = copyOn(frame, labelShift, labels{k-1}, labelAlphaColor);
    superframe  = copyOn(supersheet, sheetShift, framel);
    for s=1:stillFrames
      fprintf(repmat('\b', 1, numb)); numb = fprintf('Procesing image %03d, still frame %03d...', k-1, s);
      f         = makeframe(f, moviename, numceros, superframe, coverParams);
    end
    closerExamination(k-1);
    clear thisim;
    
    %get the next generation
    [nextim map]  = imread(images{k});
    if ~isempty(map)
      nextim      = myind2rgb(nextim, map);
    end
    nextim        = nextim(1:end-6,:,:);
    nextim        = replaceSky(nextim, imageSkyColor);
    sizenextim    = size(nextim);
    maxx          = max(sizethisim(1), sizenextim(1));
    maxy          = max(sizethisim(2), sizenextim(2));
    if ratiosheet<maxx/maxy
      scalenext   = sizesheet(1)/maxx;
    else
      scalenext   = sizesheet(2)/maxy;
    end
    [frame2 pos2] = putOnSheet(nextim, scalenext, resizeParams, sheet, sizesheet);
    
    %do a transition from the current to the next scale
    scales        = linspace(1, scalenext/scalethis, transitionFrames);
    scales        = scales(1:end);
    for s=scales
      fprintf(repmat('\b', 1, numb)); numb = fprintf('Procesing image %03d, scaling frame %s...', k-1, mat2str(s));
      [frame pos] = putOnSheet(img, s, resizeParams, sheet, sizesheet);
      framel      = copyOn(frame, labelShift, labels{k-1}, labelAlphaColor);
      superframe  = copyOn(supersheet, sheetShift, framel);
      f           = makeframe(f, moviename, numceros, superframe, coverParams);
    end

    %blend from generation k to generation k+1
    fram          = frame;
    img           = frame2((pos2(1):(pos2(1)+pos2(3)-1)), (pos2(2):(pos2(2)+pos2(4)-1)), :);
    img           = putmark(img,  markthickness, colorthis);
    frame2        = copyOn(frame2, pos2(1:2), img);
    param         = linspace(0, 1, blend2Frames);
    for s=1:blend2Frames
      fprintf(repmat('\b', 1, numb)); numb = fprintf('Procesing image %03d, blend2 frame %03d...', k-1, s);
      frame       = uint8(double(frame2)*param(s)+double(fram)*(1-param(s)));
      framel      = copyOn(frame, labelShift, labels{k-1}, labelAlphaColor);
      superframe  = copyOn(supersheet, sheetShift, framel);
      f           = makeframe(f, moviename, numceros, superframe, coverParams);
    end
    
    thisim        = nextim;
    pos           = pos2;
    scalethis     = scalenext;
    sizethisim    = sizenextim;
    clear nextim;
    clear scalenext frame2 pos2 fram framel img2 imgw param scales
  end
  
  %put the last generation
  if numel(images)>1
    framel        = copyOn(frame, labelShift, labels{end}, labelAlphaColor);
    superframe    = copyOn(supersheet, sheetShift, framel);
    for kx=1:stillFrames
      fprintf(repmat('\b', 1, numb)); numb = fprintf('Procesing image %03d, still frame %03d...', numel(images), kx);
      f           = makeframe(f, moviename, numceros, superframe, coverParams);
    end
    closerExamination(k);
  end
  fprintf('\n');
% catch ME
%   rethrow(ME);
% end


  function closerExamination(qw)
    if generations(qw)<generationdetail
    elseif all(size(screen)>0) && (numel(maxpositions)>0) && ((numel(maxpositions)>1) || (maxpositions>0))
      
      stretched   = imresize(thisim, 'OutputSize', [sizescreen(1) nan], resizeParams{:});
      npxtretched = sizescreen(2);%numpixelstomove*sizescreen(1)/sizethisim(1);
      thispos     = [pos(1:2) size(img)];
      if numel(maxpositions)>0
        ty = find(gensCloser==generations(qw));
        if isempty(ty)
          return
        end
        originalpos = originalpositions{ty};
        if isempty(originalpos)
          return
        end
        positions   = floor( (originalpos-1)*size(stretched,2)/sizethisim(2) ) +1;
      else
        positions   = round(linspace(1, size(stretched, 2) - sizescreen(2) + 1, min(maxpositions, round(size(stretched, 2)/npxtretched))));
      end
      minipositions = floor( (positions-1)*thispos(4)/size(stretched,2) ) +1;
      minilen     = floor(size(screen,2)*thispos(4)/size(stretched,2));
      for p=1:(numel(positions))
        %compose framel
        a1 = minipositions(p);
        frame2 = frame;
        a2 = min(size(frame2, 2), a1+minilen-1);
        for cc=1:3
          frame2(thispos(1):end,a1,   cc) = colorotro(cc);
          frame2(thispos(1):end,a2,   cc) = colorotro(cc);
          frame2(thispos(1),    a1:a2,cc) = colorotro(cc);
          frame2(end,           a1:a2,cc) = colorotro(cc);
        end
        frame2(thispos(1):end, a1:a2, :) = frame2(thispos(1):end, a1:a2, :)*0.75;
        framel2   = copyOn(frame2, labelShift, labels{qw}, labelAlphaColor);
        a1 = positions(p);
        a2 = (a1+sizescreen(2)-1);
        if size(stretched, 2) >= a2
          ministret = stretched(:,a1:a2,:);
        else
          ministret = [stretched(:,a1:end,:), repmat(reshape(imageSkyColor, [1 1 3]), [size(stretched, 1), size(stretched,2)-a2])]; %#ok<AGROW>
        end
        ministret = putmark(ministret, markthickness, colorotro);
        oldsuperframe = double(superframe);
        superframe = copyOn(supersheet, screenShift, ministret);
        superframe = copyOn(superframe, sheetShift,  framel2);
        newsuperframe = double(superframe);
        vl            = linspace(0, 1, closeExaminationTransitionFrames);
        for ss=1:closeExaminationTransitionFrames
          fprintf(repmat('\b', 1, numb)); numb = fprintf('Procesing image %03d, transition 1 of gliding frame %03d in position %03d of %03d...', qw, ss, p, numel(positions));
          f         = makeframe(f, moviename, numceros, uint8(newsuperframe*vl(ss)+oldsuperframe*(1-vl(ss))), coverParams);
        end
        for ss=1:stillFrames2
          fprintf(repmat('\b', 1, numb)); numb = fprintf('Procesing image %03d, gliding frame %03d in position %03d of %03d...', qw, ss, p, numel(positions));
          f         = makeframe(f, moviename, numceros, superframe, coverParams);
        end
        if p==numel(positions)
          oldsuperframe = double(superframe);
          superframe    = copyOn(supersheet, sheetShift,  frame);
          newsuperframe = double(superframe);
          for ss=1:closeExaminationTransitionFrames
            fprintf(repmat('\b', 1, numb)); numb = fprintf('Procesing image %03d, transition 2 of gliding frame %03d in position %03d of %03d...', qw, ss, p, numel(positions));
            f         = makeframe(f, moviename, numceros, uint8(newsuperframe*vl(ss)+oldsuperframe*(1-vl(ss))), coverParams);
          end
        end
      end
    end
  end

end

function [frame position] = putOnSheet(image, scale, resizeParams, sheet, sizesheet)
  stretched         = imresize(image, 'Scale', scale, resizeParams{:});
  if size(stretched,1)>sizesheet(1)
    stretched       = imresize(image, 'OutputSize', [sizesheet(1), nan], resizeParams{:});
  elseif size(stretched,2)>sizesheet(2)
    stretched       = imresize(image, 'OutputSize', [nan, sizesheet(2)], resizeParams{:});
  end
  pxthis            = sizesheet(1)-size(stretched,1)+1;
  pythis            = floor((sizesheet(2)-size(stretched,2))/2+1);
  frame             = copyOn(sheet, [pxthis pythis], stretched);
  position          = [pxthis pythis size(stretched)];
end

function f        = makeframe(f, moviename, numceros, frame, coverParams)
  if (f==0) && (~isempty(coverParams))
    f = makeCover(f, frame, @(f, image)makeframe(f, moviename, numceros, image, []), coverParams);
  end
  imwrite(frame, sprintf([moviename '%0' numceros 'd.png'], f));
  f = f+1;
end

function mark     = putmark(mark, markthickness, color)
% for z=1:3;
%   mark(1:markthickness,:,z)         = color(z);
%   mark(end-markthickness+1:end,:,z) = color(z);
%   mark(:,1:markthickness,z)         = color(z);
%   mark(:,end-markthickness+1:end,z) = color(z);
% end;
end


function o = myind2rgb(a,cm)
  if ~isa(a, 'double')
      a = a+1;    % Switch to one based indexing
  end

  if isa(cm, 'double')
      cm = uint8(cm*255);
  elseif isa(cm, 'uint8')
  else
    error('unexpected!!!!!!');
  end

  zz = double(max(max(a)));
  zb = size(cm,1);
  if zz>zb
    cm = [cm; cm(10:(9+zz-zb),:)];
  end

  % % Make sure A is in the range from 1 to size(cm,1)
  % a = max(1,min(a,size(cm,1)));


  r = zeros(size(a), 'uint8'); r(:) = cm(a,1);
  g = zeros(size(a), 'uint8'); g(:) = cm(a,2);
  b = zeros(size(a), 'uint8'); b(:) = cm(a,3);

  o = cat(3, r, g, b);
end


function thisim            = replaceSky(thisim, imageSkyColor)
  u255 = uint8(255);
  for k=1:size(thisim,1)
    l = thisim(k,:,:);
    areSky = find(all(l == u255, 3));
    nl = numel(l)/3;
    for z=1:3
      l(areSky) = imageSkyColor(z);
      areSky = areSky+nl;
    end
    thisim(k,:,:) = l;
  end
end
