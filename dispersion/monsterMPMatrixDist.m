function dists = monsterMPMatrixDist(basedir, filename, path, numTasks, type)
%for jaccard/logbrnaches distances across several simulations

fprintf('LOADING %s\n', filename);
res = load([basedir filesep filename]); %result from monsterMPMatrixDistData
res = res.res;

fprintf('INITIALIZING\n');
ni = numel(res.indivs);

if not(isnan(numTasks))
  [pairs numpairs] = getOrderedPairs(ni);
  idxpairs = randperm(numpairs)';
  pairs = pairs(idxpairs,:);
  [ranges rangesizes] = calculateRanges(numpairs, numTasks);

  pairsd = mat2cell(pairs, rangesizes, 2);
  pairsd = cellfunc(@sortrows, pairsd);

  filenames = repmat({[basedir, filesep, filename]}, size(pairsd));

  jobArgs = {'Tag', 'monsterMatrix', 'PathDependencies',  genpath(path)};
  % dists = jobArgs;
  % return;

  fprintf('GOING PARALLEL\n');
  switch type
    case 'MP'
      fun = @calculatePairsMP;
    case 'LOGNB'
      fun = @calculatePairsLOGNB;
    case 'GNM'
      fun = @calculatePairsGNM;
  end
  ress = my_dfeval(fun, filenames, pairsd, jobArgs{:});

  fprintf('SAVING POR SI LAS MOSCAS\n');
  ressfile = [basedir filesep 'monsterMatrix_ress.mat'];

  save(ressfile, 'ress');

  clear pairsd idxpairs pairs filenames res;

  fprintf('REORDERING DATA\n');
  dists = zeros(ni);
  for k=1:numel(ress)
    res = ress{k};
    i1 = res(:,1);
    i2 = res(:,2);
    ii = [i1+(i2-1)*ni; i2+(i1-1)*ni];
    dists(ii) = [res(:,3); res(:,3)];
  %   for m=1:size(res,1)
  %     a = res(m,1);
  %     b = res(m,2);
  %     r = res(m,3);
  %     dists(a,b)=r;
  %     dists(b,a)=r;
  %   end
  end
  fprintf('OK!!!!\n');
  delete(ressfile);


else

  switch type
    case 'LOGBN'
      nramas = double(res.datos(:,5));
      dists = zeros(ni);
      for a=1:(ni-1)
        if mod(a,100)==0
          fprintf('Going for %d/%d\n', a, ni);
        end
        b1 = a+1;
        b2 = ni;
        dists(a,b1:b2) = abs(log10(nramas(b1:b2))-log10(nramas(a)));
        dists(b1:b2,a) = dists(a,b1:b2)';
      end
    case 'MP'
      indivs = res.indivs;
      datos = res.datos;
      dists = zeros(ni);
      offs = datos(:,[3 4]);
      for a=1:(ni-1)
        fprintf('Going for %d\n', a);
        for b=a+1:ni
          if mod(b,500)==0
            fprintf('  In (%d,%d)\n', a,b);
          end
          dists(a,b) = 1 - similitud2(indivs{a},offs(a,:),indivs{b},offs(b,:));
          dists(b,a) = dists(a,b);
        end
      end
    case 'GNM'
      indivs = res.indivs;
      dists = zeros(ni);
      for a=1:(ni-1)
        fprintf('Going for %d\n', a);
        for b=a+1:ni
%           if mod(b,500)==0
%             fprintf('  In (%d,%d)\n', a,b);
%           end
          dists(a,b) = strdstMX(indivs{a}, indivs{b});
          dists(b,a) = dists(a,b);
        end
      end
  end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [results] = calculatePairsLOGNB(filename, pairs)

res = load(filename); %result from monsterMPMatrixDistData
res = res.res;

%[siz1 siz2 off1 off2 nramas nleafs npixelsneg]
%[01   02   03   04   05     06     07        ]
nramas = res.datos(:,5);
nramas = double(nramas(pairs));
if size(pairs,1)==1
  nramas = nramas(:)';
end
results = [pairs, abs(diff(log10(nramas),1,2))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [results] = calculatePairsMP(filename, pairs)

res = load(filename); %result from monsterMPMatrixDistData
res = res.res;

indivs = res.indivs;
results = [pairs, zeros(size(pairs,1),1)];
offs = res.datos(:,[3 4]);
for k=1:size(pairs,1)
  a = pairs(k,1);
  b = pairs(k,2);
  results(k,3) = 1 - similitud2(indivs{a},offs(a,:),indivs{b},offs(b,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [results] = calculatePairsGNM(filename, pairs)

res = load(filename); %result from monsterMPMatrixDistData
res = res.res;

indivs = res.indivs;
results = [pairs, zeros(size(pairs,1),1)];
for k=1:size(pairs,1)
  a = pairs(k,1);
  b = pairs(k,2);
  results(k,3) = strdstMX(indivs{a}, indivs{b});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pairs numpairs] = getOrderedPairs(ni)
numpairs = ni*(ni-1)/2;
k=1;
pairs = zeros(numpairs,2);
for a=1:(ni-1)
  %fprintf('Going for %d\n', a);
  nb = ni-a-1;
  pairs(k:k+nb,1) = a;
  pairs(k:k+nb,2) = (a+1:ni)';
  k=k+nb+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sim = similitud2(raster1,offset1,raster2,offset2)%, show)
%     pixel1 = horzcat(raster1{:});
%     pixel2 = horzcat(raster2{:});
%     v = isempty(pixel1) + isempty(pixel2);
%     switch v
%       case 2
%         sim = 1;
%       case 1
%         sim = 0;
%       otherwise
%         if size(pixel1,1) == 1 
%             pixel1 = [];
% 
%             pixel1(:,1) = raster1{1};
%             pixel1(:,2) = raster1{2};
%         end
%         if size(pixel2,1) == 1 
%             pixel2 = [];
% 
%             pixel2(:,1) = raster2{1};
%             pixel2(:,2) = raster2{2};
%         end
%
%         pixel1(:,1) = pixel1(:,1) - offset1(1);
%         pixel1(:,2) = pixel1(:,2) - offset1(2);
% 
%         pixel2(:,1) = pixel2(:,1) - offset2(1);
%         pixel2(:,2) = pixel2(:,2) - offset2(2);

        % sim = size(intersect(pixel1,pixel2,'rows'),1)/size(union(pixel1,pixel2,'rows'),1);

        %this is faster than the line above, but only works if the arrays of pixels
        %have no repeated rows (that is to say, if the arrays are row sets)

        sortedrows = [raster1{1} - offset1(1), raster1{2} - offset1(2); raster2{1} - offset2(1), raster2{2} - offset2(2)];
        %sort rows of both trees
%         sortedrows      = [pixel1; pixel2];
        [junk,ind]    = sort(sortedrows(:,  2));
        [junk,ind2]   = sort(sortedrows(ind,1));
        ind             = ind(ind2);
        sortedrows      = sortedrows(ind,:);
        %find matching entries
%         matching        = all(sortedrows(1:end-1,:) == sortedrows(2:end,:), 2);
%         sizeintersec    = sum(matching);
        matching        = sortedrows(1:end-1,:) == sortedrows(2:end,:);
        sizeintersec    = sum(matching(:,1) & matching(:,2));
%         sim             = sizeintersec/(size(pixel1,1)+size(pixel2,1)-sizeintersec);
        sim             = sizeintersec/(size(raster1{1},1)+size(raster2{1},1)-sizeintersec);
%     end



