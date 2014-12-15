function operators(basedir, generation, rangeid, individual)

%      2   op_alteracion        A simbolo  posicion
%      3   op_dupaleatoria      R segmento posicion
%      4   op_dupnivel          L segmento posicion
%      5   op_dupsecuencia      C segmento posicion
%      6   op_eliminacion       D posicion
%      7   op_insercion         I simbolo  posicion
%      8   op_transferencia     T segmento posicion

g = reapIndividualGenealogy(basedir, generation, rangeid, individual);

generacionesOp = g.gDescendant;
cambios = g.change;

axis([0,individual,0,8]);
tama = size(cambios,1);

vsA = zeros(1,generation);
vsR = zeros(1,generation);
vsL = zeros(1,generation);
vsC = zeros(1,generation);
vsD = zeros(1,generation);
vsI = zeros(1,generation);
vsT = zeros(1,generation);


for i = min(generacionesOp) : generation
    k = 1;
    op1 = 0;
    op2 = 0;
    op3 = 0;
    op4 = 0;
    op5 = 0;
    op6 = 0;
    op7 = 0;
    op8 = 0;

    sA = 0;
    sR = 0;
    sL = 0;
    sC = 0;
    sD = 0;
    sI = 0;
    sT = 0;
    if ~isempty(find(generacionesOp == i))
        p = find(generacionesOp == i);
        total =  ~isempty(find(cambios{p}== 'A')) + ~isempty(find(cambios{p}== 'R')) + ~isempty(find(cambios{p}== 'L')) + ~isempty(find(cambios{p}== 'C')) + ~isempty(find(cambios{p}== 'D')) + ~isempty(find(cambios{p}== 'I')) + ~isempty(find(cambios{p}== 'T'));
        while total > 0

            if ~isempty(find(cambios{p}== 'A'))&& op7 == 0;
                vsA(i+1) = sA + 1;
                ejeY{i+1}(k)= 7;
                k = k + 1;
                total = total - 1;
                op7 = 1;
            elseif ~isempty(find(cambios{p}== 'R'))&& op2 == 0;
                vsR(i+1) = sR + 1;
                ejeY{i+1}(k)= 2;
                k = k + 1;
                total = total - 1;
                op1 = 2;
            elseif ~isempty(find(cambios{p}== 'L'))&& op3 == 0;
                vsL(i+1) = sL + 1;
                ejeY{i+1}(k)= 3;
                k = k + 1;
                total = total - 1;
                op3 = 1;
            elseif ~isempty(find(cambios{p}== 'C'))&& op4 == 0;
                vsC(i+1) = sC + 1;
                ejeY{i+1}(k)= 4;
                k = k + 1;
                total = total - 1;
                op4 = 1;
            elseif ~isempty(find(cambios{p}== 'D'))&& op5 == 0;
                vsD(i+1) = sD + 1;
                ejeY{i+1}(k)= 5;
                k = k + 1;
                total = total - 1;
                op15 = 1;
            elseif ~isempty(find(cambios{p}== 'I'))&& op6 == 0;
                vsI(i+1) = sI + 1;
                ejeY{i+1}(k)= 6;
                k = k + 1;
                total = total - 1;
                op6 = 1;
            elseif ~isempty(find(cambios{p}== 'T'))&& op1 == 0;
                vsT(i+1) = sT + 1;
                ejeY{i+1}(k)= 1;
                k = k + 1;
                total = total - 1;
                op1 = 1;
            end
        end
    else
        ejeY{i+1}(k)= 0;
    end
end



sumvsA = cumsum(vsA);
sumvsR = cumsum(vsR);
sumvsL = cumsum(vsL);
sumvsC = cumsum(vsC);
sumvsD = cumsum(vsD);
sumvsI = cumsum(vsI);
sumvsT = cumsum(vsT);

plot(sumvsA,'b')
hold on
plot(sumvsR,'g')
hold on
plot(sumvsL,'r')
hold on
plot(sumvsC,'c')
hold on
plot(sumvsD,'m')
hold on
plot(sumvsI,'y')
hold on
plot(sumvsT,'k')

xlabel('generations')
ylabel('cumulative sum of the operators')
legend('opalteracion','opdupaleatoria','opdupnivel','opdupsecuencia','opeliminacion','opinsercion','optransferencia')

h2 = figure;
for i =  0 : generation - 1
    sejeY = length(ejeY{i+1});
    for j = 1 : sejeY
        plot(i,ejeY{i+1}(j),'m*') % puntos
        hold on
    end
end

xlabel('generations')
ylabel('operators')
set(gca,'YTickLabel','None|op_transferencia|op_dupaleatoria|op_dupnivel|op_dupsecuencia|op_eliminacion|op_insercion|op_alteracion');

% for i = min(generacionesOp) : generation
%     k = 1;
%     op1 = 0;
%     op2 = 0;
%     op3 = 0;
%     op4 = 0;
%     op5 = 0;
%     op6 = 0;
%     op7 = 0;
%     op8 = 0;
%     if ~isempty(find(generacionesOp == i))
%         p = find(generacionesOp == i);
%         total = ~isempty(find(cambios{p}== '=')) + ~isempty(find(cambios{p}== 'A')) + ~isempty(find(cambios{p}== 'R')) + ~isempty(find(cambios{p}== 'L')) + ~isempty(find(cambios{p}== 'C')) + ~isempty(find(cambios{p}== 'D')) + ~isempty(find(cambios{p}== 'I')) + ~isempty(find(cambios{p}== 'T'));
%         while total > 0
%             if ~isempty(find(cambios{p}== '=')) && op1 == 0;
%                 ejeY{i+1}(k)= 1;
%                 k = k + 1;
%                 total = total - 1;
%                 op1 = 1;
%             elseif ~isempty(find(cambios{p}== 'A'))&& op8 == 0;
%                 ejeY{i+1}(k)= 8;
%                 k = k + 1;
%                 total = total - 1;
%                 op2 = 1;
%             elseif ~isempty(find(cambios{p}== 'R'))&& op3 == 0;
%                 ejeY{i+1}(k)= 3;
%                 k = k + 1;
%                 total = total - 1;
%                 op3 = 1;
%             elseif ~isempty(find(cambios{p}== 'L'))&& op4 == 0;
%                 ejeY{i+1}(k)= 4;
%                 k = k + 1;
%                 total = total - 1;
%                 op4 = 1;
%             elseif ~isempty(find(cambios{p}== 'C'))&& op5 == 0;
%                 ejeY{i+1}(k)= 5;
%                 k = k + 1;
%                 total = total - 1;
%                 op5 = 1;
%             elseif ~isempty(find(cambios{p}== 'D'))&& op6 == 0;
%                 ejeY{i+1}(k)= 6;
%                 k = k + 1;
%                 total = total - 1;
%                 op16 = 1;
%             elseif ~isempty(find(cambios{p}== 'I'))&& op7 == 0;
%                 ejeY{i+1}(k)= 7;
%                 k = k + 1;
%                 total = total - 1;
%                 op7 = 1;
%             elseif ~isempty(find(cambios{p}== 'T'))&& op2 == 0;
%                 ejeY{i+1}(k)= 2;
%                 k = k + 1;
%                 total = total - 1;
%                 op8 = 1;
%             end
%         end
%     else
%         ejeY{i+1}(k)= 0;
%     end
% end
%
% for i =  0 : generation - 1
%     sejeY = length(ejeY{i+1});
%     for j = 1 : sejeY
%         plot(i,ejeY{i+1}(j),'m*') % puntos
%         hold on
%     end
% end
%
% xlabel('generations')
% ylabel('operators')
% set(gca,'YTickLabel','None|op_without_change|op_transferencia|op_dupaleat
% oria|op_dupnivel|op_dupsecuencia|op_eliminacion|op_insercion|op_alteracio
% n');