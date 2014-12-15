function [esBalanceado] = corchetesBalanceados(regla)
% Chequea que la regla tiene los corchetes balanceados

cor = zeros(1, length(regla));
cor(regla=='[') = 1;
cor(regla==']') = -1;
cor = cumsum(cor);
esBalanceado = (numel(cor)>0) && (~any(cor < 0) && cor(end) == 0);
