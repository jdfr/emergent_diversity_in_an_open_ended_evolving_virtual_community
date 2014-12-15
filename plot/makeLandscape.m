function makeLandscape(basedir, gen, placeTallTrees, namefile)
%makeLandscape('rene/dataset/G0.5F1/uno/1=0.001_1', 79, true, '');
%makeLandscape('rene/dataset/G0.5F1/tres/1=0.001_1', 128, true, '');
if ~exist('placeTallTrees', 'var')
  placeTallTrees = [];
end
if ~exist('namefile', 'var')
  namefile = '';
end

params = load([basedir filesep '..' filesep 'estado.mat']);
params = params.ACIParams;
params = putDefaultParameters(params);
nm = sprintf('P%03d', gen);
P = load([basedir filesep nm '.mat']);
P = P.(nm);
if not(isempty(placeTallTrees))
  params.tooTallArePlaced = placeTallTrees;
end

if isempty(namefile)
  fitness2(P,params,gen,[],true,P.randN,struct('overridePlotting', true, 'basedir', [], 'namefile', []));
else
  fitness2(P,params,gen,[],false,P.randN,struct('overridePlotting', true, 'basedir', basedir, 'namefile', namefile));
end

