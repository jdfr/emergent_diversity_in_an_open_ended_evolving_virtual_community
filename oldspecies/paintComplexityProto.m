function paintComplexityProto(basedir,GenI,GenF)

stdt = zeros(GenF,3);
if ispc 
  fsep = '\\';
else
  fsep = '/';
end
for i = GenI : GenF
    if i == 33
        disp('hola')
    end
    load(sprintf(['%s' fsep 'comLZfigureLongitud%03d.mat'],basedir,i));
    comLZfigureLongitud = eval(sprintf('comLZfigureLongitud%03d', i));
    if length(comLZfigureLongitud) == 1
        T{i}=1;

    else
        Y = distancia(comLZfigureLongitud);
        Z = linkage(Y,'single');
        T{i} = cluster(Z,'maxclust',3);
    end
    a = comLZfigureLongitud(find(T{i}==1));
    mediaa = mean(a);
    b = comLZfigureLongitud(find(T{i}==2));
    mediab = mean(b);
    c = comLZfigureLongitud(find(T{i}==3));
    mediac = mean(c);

    medias = [mediaa mediab mediac];
    [ordenar pos] = sort(medias);

    if ~isempty(find(T{i}==pos(1)))

        stdt(i,1) = max(comLZfigureLongitud(find(T{i}==pos(1))));
    else
        stdt(i,1) = 0;
    end

    if ~isempty(find(T{i}==pos(2)))

        stdt(i,2) = max(comLZfigureLongitud(find(T{i}==pos(2))));
    else
        stdt(i,2) = 0;
    end
    if ~isempty(find(T{i}==pos(3)))

        stdt(i,3) = max(comLZfigureLongitud(find(T{i}==pos(3))));
    else
        stdt(i,3) = 0;
    end


    colors = [1 0 0; 0 1 0; 0 0 1];

    colorsD(1,:) = colors(find(pos == 1),:);
    colorsD(2,:) = colors(find(pos == 2),:);
    colorsD(3,:) = colors(find(pos == 3),:);



    if ~isempty(find(T{i}==1))
        plot(i,comLZfigureLongitud(find(T{i}==1)),'*','Color',colorsD(1,:))
        hold on
    end
    if ~isempty(find(T{i}==2))
        plot(i,comLZfigureLongitud(find(T{i}==2)),'*','Color',colorsD(2,:))
        hold on
    end
    if ~isempty(find(T{i}==3))
        plot(i,comLZfigureLongitud(find(T{i}==3)),'*','Color',colorsD(3,:))
        hold on
    end

    %     com(i) = mean(comLZfigureLongitud);
        maxi(i) = max(comLZfigureLongitud);
    %     mini(i) = min(comLZfigureLongitud);
    % %     varianza(i) = var(comLZfigureLongitud);
    %     desviacionTipica(i) = std(comLZfigureLongitud);
end

% plot(stdt(:,1),'m')
% hold on
% plot(stdt(:,2),'c')
% hold on
% 
%h = figure;

fit = polyfit(1:GenF,stdt(:,1)',1);
plot(1:GenF,polyval(fit,1:GenF),'k');
hold on
fit = polyfit(1:GenF,stdt(:,2)',1);
plot(1:GenF,polyval(fit,1:GenF),'k');
hold on
%axis([0 250 0 20000])

% 
% y1 = stdt(:,1)';
% x1 = 1:GenF;
% Y1 = log(y1);
% M1 = [length(x1),sum(x1);sum(x1),sum(x1.^2)];
% u1 = [sum(Y1);sum(x1.*Y1)];
% sol1 = M1\u1;
% a1 = exp(sol1(1));
% b1 = sol1(2);

% pol1 = a1*exp(b1*x1);
% plot(x1,pol1,'k');
% hold on
% 
% y2 = stdt(:,2)';
% noceros = find(y2~=0);
% y2 = y2(noceros);
% x2 = x1(noceros);
% Y2 = log(y2);
% M2 = [length(x2),sum(x2);sum(x2),sum(x2.^2)];
% u2 = [sum(Y2);sum(x2.*Y2)];
% sol2 = M2\u2;
% a2 = exp(sol2(1));
% b2 = sol2(2);
% 
% pol2 = a2*exp(b2*x1);
% plot(x1,pol2,'k');
% 


% hold on
% plot(stdt(:,3),'k')



% mediaDes = mean(desviacionTipica);
% [maxDes pos]= max(desviacionTipica);
% x1 = 0;
% y1 = 0;
% x21 = i;
% x22 = pos;
% y21 = mediaDes;
% y22 = maxDes;
% % r1 = (y21-y1/x2-y1)*(GenF-x1)+y1;
% % r2 = (y22-y1/x2-y1)*(GenF-x1)+y1;
% m2 = y22/x22;
% r2 = GenF*m2;
%
% plot(com)
% hold on
% plot(maxi,'g')
% hold on
% plot(mini,'r')
% hold on
% % plot(varianza,'m')
% % hold on
% plot(desviacionTipica,'m')
% hold on
% plot([x1 x21],[y1 y21],'k')
% hold on
% plot([x1 GenF],[y1 r2],'k')


% xlabel('generation')
% ylabel('complexity of the prototypes of each species')
% legend('mean','maximum','minimum','standar derivation')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dist=distancia(X)

long = length(X);
k = 1;
for i = 1 : long-1
    for j = (i+1): long
        dist(k) = abs(X(j) - X(i));
        k = k + 1;
    end
end
