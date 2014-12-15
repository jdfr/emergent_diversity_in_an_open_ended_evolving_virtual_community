function matriz = numgenerationsExhaustive(fileawk, matriz, doShow, foutput)

if isempty(fileawk)
  %basedir = '\\150.214.56.111\l-system\l-system\codigo\codigo\alphaGalphaF';
  %lsystem@orus:~$ find ./l-system/codigo/codigo/alphaGalphaF -type f -name poblacion.txt -exec awk 'END {print FILENAME "," $1 "," $3}' '{}' \;   
  fileawk = 'fileExhaustiveAWK.txt';
end

if not(exist('foutput', 'var')) 
  foutput = [];%'phasemap.png';
end

  gs = 0.1:0.1:1;
  fs = 0.1:0.1:1;
  is = 1:5;

if isempty(matriz)

  matriz = struct('numgens', nan(numel(gs), numel(fs), numel(is)), 'poplastgen', nan(numel(gs), numel(fs), numel(is)), 'gs', gs, 'fs', fs, 'is', is);
  mtr = matriz.numgens;
  lastg = matriz.poplastgen;
  lns = cell(1000,1);
  lq = 1;
  
  fh = fopen(fileawk, 'r');
  while not(feof(fh))
    ln = fgetl(fh);
    lns{lq}=ln;
    lq=lq+1;
    %G0.6F0.6/2009_Jul_03_17_20_56/3/,220,49
    sep = find(ln=='/');
    comma = find(ln==',');
    G = find(ln=='G');
    F = find(ln=='F');
    sg = str2double(ln(G+1:F-1));
    sf = str2double(ln(F+1:sep(1)-1));
    si = str2double(ln(sep(2)+1:sep(3)-1));
    sng = str2double(ln(comma(1)+1:comma(2)-1));
    snp = str2double(ln(comma(2)+1:end));
    g = find(abs(sg-gs)<1e-6);
    f = find(abs(sf-fs)<1e-6);
    i = find(abs(si-is)<1e-6);
    if isempty(g)
      error('mira g: %s', sg);
    end
    if isempty(f)
      error('mira f: %s', sf);
    end
    if isempty(i)
      error('mira i: %s', si);
    end
    if isnan(sng)
      error('mira numgen for g=%d,f=%d,i=%d: <%s>', sg, sf, si, ln(comma(1)+1:comma(2)-1));
    end
    if isnan(snp)
      error('mira lastnumpop for g=%d,f=%d,i=%d: <%s>', sg, sf, si, ln(comma(2)+1:end));
    end
    if not(isnan(mtr(g,f,i)))
      fprintf('Already loaded g=%f,f=%f,i=%d\n', sg, sf, si);
    else
      mtr(g,f,i) = sng;
      lastg(g,f,i) = snp;
    end
  end
  
  fclose(fh);
  
  matriz.numgens = mtr;
  matriz.poplastgen = lastg;

end

if doShow
  numObs  = sum(not(isnan(matriz.numgens)),3);
  [ind] = find(isnan(matriz.numgens));
  for k=1:numel(ind)
    [a b c]= ind2sub(size(matriz.numgens), ind(k));
    %this is until we rerun these experiments
    matriz.numgens(a,b,c) = max(matriz.numgens(a,b,:));
    matriz.poplastgen(a,b,c) = max(matriz.poplastgen(a,b,:));
  end
  medgens = median(matriz.numgens,3);
  mengens =   mean(matriz.numgens,3);
  maxgens =   max(matriz.numgens,[],3);
  medpop = median(matriz.poplastgen,3);
  menpop =   mean(matriz.poplastgen,3);
  maxpop =   max(matriz.poplastgen,[],3);
  
  gmp = gray(512);
  %gmp = gmp(end:-1:end/2,:);
  gmp = gmp(end/2:end,:);
  xlab = '{\alpha} parameter';
  ylab = '{\beta} parameter';
  toimg = [];
  
%   showimg(toimg, gmp, gs, fs, numObs',  xlab, ylab, 'number of observations');
%   showimg(toimg, gmp, gs, fs, maxgens', xlab, ylab, 'max last generation');
  res1 = showimg(toimg, gmp, gs, fs, medgens', xlab, ylab, 'median last generation');
%   showimg(toimg, gmp, gs, fs, mengens', xlab, ylab, 'mean last generation');
%   showimg(toimg, gmp, gs, fs, maxpop',  xlab, ylab, 'max number of individuals in last generation');
  res2 = showimg(toimg, gmp, gs, fs, medpop',  xlab, ylab, 'median number of individuals in last generation');
%   showimg(toimg, gmp, gs, fs, menpop',  xlab, ylab, 'mean number of individuals in last generation');
  if toimg
    recA1 = 400;
    recA2 = 400;
    recB1 = 400;
    recB2 = 400;
    res = [res1(:,recA1:end-recA2,:), res2(:,recB1:end-recB2,:)];
    imwrite(res, foutput);
  end
end


function res = showimg(showres, colmap,gs, fs, img, xlab, ylab, tit)

if isempty(showres)
  h = figure;
else
  h = figure('Visible', 'off');
end
set(h, 'colormap', colmap)
imagesc(gs, fs, img); colorbar;
a = get(h, 'currentaxes');
set(a, 'FontName','times','FontSize',14,...
       'XTick', gs, 'YTick', fs);
axis xy;
axis square;
xlabel(xlab);
ylabel(ylab);
title(tit);

if isempty(showres)
  res = [];
else
  res = imcapture(h, 'all', 600);
  close(h);
end

