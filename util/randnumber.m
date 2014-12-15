function [n] = randnumber(interv)
% genera un vector de números aleatorios en el rango especificado por la
% matriz interv, que contiene el inicio del rango en la primera columna y el
% fin en la segunda
%
% >> randnumber([0 1;3 4;20 30])'
% ans =
%     0.5977    3.1240   21.5666


n = interv(:,1) + rand(size(interv,1),1).*(interv(:,2) - interv(:,1));