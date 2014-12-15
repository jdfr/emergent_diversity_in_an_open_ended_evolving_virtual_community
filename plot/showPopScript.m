function showPopScript(option)


%combs = {init1 end1 height1 col1 width1; init2 end2 height2 col2 width2;...}

switch (option)
  case 1
    imgname = 'foto_1_veryharsh.png';
    imgwrite = 'detail_1_veryharsh.png';
    width = 1500;
    len = 10733;
    nparts = 8;
    sep = 6;
    poss = round(linspace(1, len, nparts+1))';
    inits = array2cell(poss(1:end-1));
    ends  = array2cell(poss(2:end));
    inits = array2cell([1 1365 2699 4026 5367 6709 8050 9392]'); %[1 1343 2684 4026 5367 6709 8050 9392]
    ends  = array2cell([1365 2699 4026 5367 6709 8050 9392 10733]'); %[1343 2684 4026 5367 6709 8050 9392 10733]
    heights = array2cell(nan(size(inits)));
    heights = array2cell([55 35 105 90 1 90 105 130]'); %rows as measured in mspaint (or your favourite image editor
    widths = array2cell(nan(size(inits)));
    a = 255; b = 250;
    cols = {[b a a]; [a b a]; [a a b]; [b b a]; [b a b]; [a b b]; [a a a]; [b b b]};
    combs = [inits ends heights cols(1:numel(inits)) widths];
    doReplaceSky = true;
    showUpper = true;
    showDetailedImage(imgname, imgwrite, sep, width, combs, doReplaceSky, showUpper);
  case 1.1
    imgname = 'entorno_500_G1F1.png';
    imgwrite = 'detail_500_G1F1.png';
    width = 1500;
    len = 8560;
    nparts = 8;
    sep = 6;
    poss = round(linspace(1, len, nparts+1))';
    inits = array2cell(poss(1:end-1));
    ends  = array2cell(poss(2:end));
    inits = array2cell([1 1071 2141 3211 4270 5350 6420 7490]'); %[1 1071 2141 3211 4281 5350 6420 7490]
    ends  = array2cell([1071 2141 3211 4270 5350 6420 7490 8560]'); %[1071 2141 3211 4281 5350 6420 7490 8560]
    heights = array2cell(nan(size(inits)));
    heights = array2cell([100 100 100 100 1 100 100 100]'); %rows as measured in mspaint (or your favourite image editor
    widths = array2cell(nan(size(inits)));
    a = 255; b = 250;
    cols = {[b a a]; [a b a]; [a a b]; [b b a]; [b a b]; [a b b]; [a a a]; [b b b]};
    combs = [inits ends heights cols(1:numel(inits)) widths];
    doReplaceSky = true;
    showUpper = true;
    showDetailedImage(imgname, imgwrite, sep, width, combs, doReplaceSky, showUpper);
  case 1.2
    imgname = 'entornoRaster500_G1F1_3.png';
    imgwrite = 'detail_3_veryharsh.png';
    width = 1500;
    len = 10372;
    nparts = 8;
    sep = 6;
    poss = round(linspace(1, len, nparts+1))';
    inits = array2cell(poss(1:end-1));
    ends  = array2cell(poss(2:end));
%     inits = array2cell([1 1365 2699 4026 5367 6709 8050 9392]'); %[1 1343 2684 4026 5367 6709 8050 9392]
%     ends  = array2cell([1365 2699 4026 5367 6709 8050 9392 10733]'); %[1343 2684 4026 5367 6709 8050 9392 10733]
    heights = array2cell(nan(size(inits)));
%     heights = array2cell([55 35 105 90 1 90 105 130]'); %rows as measured in mspaint (or your favourite image editor
    widths = array2cell(nan(size(inits)));
    a = 255; b = 250;
    cols = {[b a a]; [a b a]; [a a b]; [b b a]; [b a b]; [a b b]; [a a a]; [b b b]};
    combs = [inits ends heights cols(1:numel(inits)) widths];
    doReplaceSky = true;
    showUpper = true;
    showDetailedImage(imgname, imgwrite, sep, width, combs, doReplaceSky, showUpper);
  case 2
    imgname = 'foto_2_harsh.png';
    imgwrite = 'detail_2_harsh.png';
    width = 1500;
    len = 15174;
    nparts = 8;
    sep = 6;
    poss = round(linspace(1, len, nparts+1))';
    inits = array2cell(poss(1:end-1));
    ends  = array2cell(poss(2:end));
    inits = array2cell([1 1898 3860 5691 7588 9430 11381 13277]'); %[1 1898 3794 5691 7588 9484 11381 13277]
    ends  = array2cell([1898 3860 5691 7588 9430 11381 13277 15174]'); %[1898 3794 5691 7588 9484 11381 13277 15174]
    heights = array2cell(nan(size(inits)));
    heights = array2cell([390 260 375 210 1 100 315 265]'); %rows as measured in mspaint (or your favourite image editor
    widths = array2cell(nan(size(inits)));
    a = 255; b = 250;
    cols = {[b a a]; [a b a]; [a a b]; [b b a]; [b a b]; [a b b]; [a a a]; [b b b]};
    combs = [inits ends heights cols(1:numel(inits)) widths];
    doReplaceSky = true;
    showUpper = true;
    showDetailedImage(imgname, imgwrite, sep, width, combs, doReplaceSky, showUpper);
  case 2.1
    imgname = 'entorno_500_G0.75F1.png';
    imgwrite = 'detail_500_G0.75F1.png';
    width = 1500;
    len = 18685;
    nparts = 3;
    sep = 6;
    poss = round(linspace(1, len, nparts+1))';
    inits = array2cell(poss(1:end-1));
    ends  = array2cell(poss(2:end));
    inits = array2cell([1 6229 14100]'); %[1 6229 12457]
    ends  = array2cell([6229 14100 18685]'); %[6229 12457 18685]
    heights = array2cell(nan(size(inits)));
    heights = array2cell([1850 1 2850]'); %rows as measured in mspaint (or your favourite image editor
    widths = array2cell(nan(size(inits)));
    a = 255; b = 250;
    cols = {[b a a]; [a b a]; [a a b]; [b b a]; [b a b]; [a b b]; [a a a]; [b b b]};
    combs = [inits ends heights cols(1:numel(inits)) widths];
    doReplaceSky = true;
    showUpper = true;
    showDetailedImage(imgname, imgwrite, sep, width, combs, doReplaceSky, showUpper);
  case 3
    imgname = 'lsystemdani\dataset\G1F1\uno\1=0.001_3\entornoRaster100.png';
    imgwrite = 'lsystemdani\src\entornoRaster100_G1F1_parts.png';
    width = 1500;
    len = 1926;
    nparts = 2;
    sep = 6;
    poss = round(linspace(1, len, nparts+1))';
    inits = array2cell(poss(1:end-1));
    ends  = array2cell(poss(2:end));
%     inits = array2cell([1 1365 2699 4026 5367 6709 8050 9392]'); %[1 1343 2684 4026 5367 6709 8050 9392]
%     ends  = array2cell([1365 2699 4026 5367 6709 8050 9392 10733]'); %[1343 2684 4026 5367 6709 8050 9392 10733]
    heights = array2cell(nan(size(inits)));
%     heights = array2cell([55 35 105 90 1 90 105 130]'); %rows as measured in mspaint (or your favourite image editor
    widths = array2cell(nan(size(inits)));
    a = 255; b = 250;
    cols = {[b a a]; [a b a]; [a a b]; [b b a]; [b a b]; [a b b]; [a a a]; [b b b]};
    combs = [inits ends heights cols(1:numel(inits)) widths];
    doReplaceSky = false;
    showUpper = false;
    showDetailedImage(imgname, imgwrite, sep, width, combs, doReplaceSky, showUpper);
  case 3.1
    imgname = 'lsystemdani\dataset\G1F1\uno\1=0.001_3\entornoRaster340.png';
    imgwrite = 'lsystemdani\src\entornoRaster340_G1F1_parts.png';
    width = 1500;
    len = 7220;
    nparts = 8;
    sep = 6;
    poss = round(linspace(1, len, nparts+1))';
    inits = array2cell(poss(1:end-1));
    ends  = array2cell(poss(2:end));
%     inits = array2cell([1 1365 2699 4026 5367 6709 8050 9392]'); %[1 1343 2684 4026 5367 6709 8050 9392]
%     ends  = array2cell([1365 2699 4026 5367 6709 8050 9392 10733]'); %[1343 2684 4026 5367 6709 8050 9392 10733]
    heights = array2cell(nan(size(inits)));
%     heights = array2cell([55 35 105 90 1 90 105 130]'); %rows as measured in mspaint (or your favourite image editor
    widths = array2cell(nan(size(inits)));
    a = 255; b = 250;
    cols = {[b a a]; [a b a]; [a a b]; [b b a]; [b a b]; [a b b]; [a a a]; [b b b]};
    combs = [inits ends heights cols(1:numel(inits)) widths];
    doReplaceSky = false;
    showUpper = false;
    showDetailedImage(imgname, imgwrite, sep, width, combs, doReplaceSky, showUpper);
  case 3.2
    imgname = 'lsystemdani\dataset\G0.75F1\uno\1=0.001_1\entornoRaster100.png';
    imgwrite = 'lsystemdani\src\entornoRaster100_G0.75F1_parts.png';
    width = 1500;
    len = 1730;
    nparts = 2;
    sep = 6;
    poss = round(linspace(1, len, nparts+1))';
    inits = array2cell(poss(1:end-1));
    ends  = array2cell(poss(2:end));
%     inits = array2cell([1 1365 2699 4026 5367 6709 8050 9392]'); %[1 1343 2684 4026 5367 6709 8050 9392]
%     ends  = array2cell([1365 2699 4026 5367 6709 8050 9392 10733]'); %[1343 2684 4026 5367 6709 8050 9392 10733]
    heights = array2cell(nan(size(inits)));
%     heights = array2cell([55 35 105 90 1 90 105 130]'); %rows as measured in mspaint (or your favourite image editor
    widths = array2cell(nan(size(inits)));
    a = 255; b = 250;
    cols = {[b a a]; [a b a]; [a a b]; [b b a]; [b a b]; [a b b]; [a a a]; [b b b]};
    combs = [inits ends heights cols(1:numel(inits)) widths];
    doReplaceSky = false;
    showUpper = false;
    showDetailedImage(imgname, imgwrite, sep, width, combs, doReplaceSky, showUpper);
  case 3.3
    imgname = 'lsystemdani\dataset\G0.75F1\uno\1=0.001_1\entornoRaster340.png';
    imgwrite = 'lsystemdani\src\entornoRaster340_G0.75F1_parts.png';
    width = 1500;
    len = 12627;
    nparts = 3;
    sep = 6;
    poss = round(linspace(1, len, nparts+1))';
    inits = array2cell(poss(1:end-1));
    ends  = array2cell(poss(2:end));
%     inits = array2cell([1 1365 2699 4026 5367 6709 8050 9392]'); %[1 1343 2684 4026 5367 6709 8050 9392]
%     ends  = array2cell([1365 2699 4026 5367 6709 8050 9392 10733]'); %[1343 2684 4026 5367 6709 8050 9392 10733]
    %heights = array2cell(nan(size(inits)));
    heights = array2cell([905 1 1500]'); %rows as measured in mspaint (or your favourite image editor
    widths = array2cell(nan(size(inits)));
    a = 255; b = 250;
    cols = {[b a a]; [a b a]; [a a b]; [b b a]; [b a b]; [a b b]; [a a a]; [b b b]};
    combs = [inits ends heights cols(1:numel(inits)) widths];
    doReplaceSky = false;
    showUpper = false;
    showDetailedImage(imgname, imgwrite, sep, width, combs, doReplaceSky, showUpper);
  case 3.4
    imgname = 'lsystemdani\dataset\G0.75F1\uno\1=0.001_1\entornoRaster340.png';
    imgwrite = 'lsystemdani\src\entornoRaster340_G0.75F1_parts4.png';
    finalwidth = 1500;
    scale = 0.46;
    %heights = nan(10,1);
    heights = [1180 1 1560 1930]'; %rows as measured in mspaint (or your favourite image editor
    showDetailedImageByLines(imgname, imgwrite, scale, finalwidth, heights);
  otherwise
    error('option not recognized!!!!');
end

