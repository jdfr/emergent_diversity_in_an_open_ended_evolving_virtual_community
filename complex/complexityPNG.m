function sizePNG = complexityPNG(basedir,numInd)

tama = dir(basedir);
for i = 3 : numInd+2
    sizePNG(i-2) = tama(i).bytes;
end