function array = genome2array(genome, mode)
%give the 0/1 bitmap of a tree (leaves not shown)

% poph1 = loadSimulation('trees\il2\bestia\variable_pop_retune\randomSearch\2009_Jun_20_13_18_45_A=04.6_T=2_O=462.2_G=0.768_F=0.473\1=0.0001_1');    
% thisg = find(poph1.generation==300); thisg(poph1.nPixelsInSoil(thisg)>0) = [];    
% arrays = genome2array(poph1.genome(thisg([1 3 8 26 27])), struct('shadowing', 'canonical', 'angle',poph1.ACIParams.T.angle, 'justPaint', true));   
% allarray = paintTrees(arrays([5 2 3 4 1]), [], false); imshow(allarray);
% imwrite(allarray, 'anglesmall.png');
%
% poph1 = loadSimulation('trees\il2\bestia\variable_pop_retune\randomSearch\2009_Jun_21_00_21_24_A=41.2_T=2_O=115.6_G=0.821_F=0.099\1=0.0001_1');
% thisg = find(poph1.generation==400); thisg(poph1.nPixelsInSoil(thisg)>0) = []; thisg(poph1.height(thisg)<200) = []; 
% arrays = genome2array(poph1.genome(thisg([20 29 47 66 67 81 84 85 87 90])), struct('shadowing', 'canonical', 'angle',poph1.ACIParams.T.angle, 'justPaint', true));   
% allarray = paintTrees(arrays([3 4 6 10]), [], false); imshow(allarray);
% imwrite(allarray, 'anglelarge.png');

if iscell(genome)
  array = cellfunc(@(x)genome2array(x, mode), genome);
else
  raster=ls2(genome, mode);
  dimension = [max(raster{1}) max(raster{2})];
  array = drawraster(raster, dimension);
end



function array = drawraster(raster, dimensions)
ys = raster{1}; if size(ys,2)>1; ys = ys(:); end;
xs = raster{2}; if size(xs,2)>1; xs = xs(:); end;
%array = accumarray([ys,xs],uint16(indexInColors+2), uint16(dimensions));
array = true(dimensions);
array(sub2ind(dimensions, ys, xs)) = false;
array = array(end:-1:1,:);
