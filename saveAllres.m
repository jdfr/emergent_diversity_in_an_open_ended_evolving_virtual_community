function saveAllres(basedir)
%b = 'rene/dataset/test'; bd = cellfunc(@(x) [b x], {'110NB_old', '110NB', 'LONG', 'LONGG', 'LONGOLD', 'LONGOLD2', 'OKPRFX'}); job = hmy_dfevalasyncJM(jm, makedummy(@saveAllres), 1, bd, 'Tag', 'saveAllres', 'PathDependencies',  genpath('rene/src'));
allres = showFitnessPerGenRec(false, true, false, basedir); %#ok<NASGU>
save([basedir filesep 'allres.mat'], 'allres');


