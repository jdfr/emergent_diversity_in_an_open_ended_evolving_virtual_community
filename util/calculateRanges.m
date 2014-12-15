function [ranges rangesizes] = calculateRanges(popsize, numRanges)
rangesizes  = diff(floor(linspace(0, popsize, numRanges+1)));
rangesizes  = rangesizes(rangesizes~=0);
ranges      = [0; cumsum(rangesizes)'];
ranges      = [ranges(1:(end-1))+1, ranges(2:end)];
