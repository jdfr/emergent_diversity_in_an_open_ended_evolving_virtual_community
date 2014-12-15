function numeroEspecies(basedir,G)

% G = number of generations

x = zeros(1,G);
for i = 1 : G
    load(sprintf('%s/numGroups%03d.mat',basedir,i));
    x(i)=eval(sprintf('numGroups%03d', i));
end

plot(x);
xlabel('generations')
ylabel('number of groups')



