%%%%%%%%%%%%%%%%%%%%%% REMOTE FUNCTION (TO EXECUTE IN CLUSTER)
function [timeSpent, individuals,dimensions,offset,maxiy,leaf,counts, mingenomes] = evaluar(P,returnIndividuals,T)

% Evalua una poblaci�n calculando el �rbol y la distancia a la relaci�n de
% aspecto buscada.
%
%  e = evaluar(P,r,ramasojb)
%      P : poblaci�n (lista de cadenas)
%      r : relaci�n de aspecto objetivo
%      ramasobj : numero de ramas requeridas
timeStart   = cputime;
n           = numel(P);
counts      = zeros(n,3);
offset      = zeros(n,2);
maxiy       =  cell(n,1);
leaf        =  cell(n,1);
dimensions  = zeros(n, 2);
individuals =  cell(n,1);
mingenomes  = cell(n,1);
retInds     = returnIndividuals>0;

for i=1:n
    if retInds && (i<=returnIndividuals)
      [individuals{i}, offset(i,:), counts(i,:), dimensions(i,:), leaf{i}, maxiy{i}, mingenomes{i}] = ls2(P{i},T, true, true); %#ok<ASGLU> VERSION CON TREERASTER
    else
      [nevermind,      offset(i,:), counts(i,:), dimensions(i,:), leaf{i}, maxiy{i}, mingenomes{i}] = ls2(P{i},T, true, true); %#ok<ASGLU> VERSION CON TREERASTER
      clear nevermind;
    end
end
timeSpent = cputime-timeStart;
end