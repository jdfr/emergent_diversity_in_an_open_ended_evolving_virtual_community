function showDisruptionHist(meanDisrup, meanDisrupR, plotXTickLabels)

Y = [meanDisrup, meanDisrupR];


%Y = [meanDisrupGenomeA meanDisrupGenomeRA; meanDisrupGenomeE meanDisrupGenomeRE; meanDisrupGenomeI meanDisrupGenomeRI; meanDisrupGenomeDA meanDisrupGenomeRDA; meanDisrupGenomeDN meanDisrupGenomeRDN; meanDisrupGenomeDAE meanDisrupGenomeRDAE; meanDisrupGenomeDNE meanDisrupGenomeRDNE; meanDisrupGenomeDSE meanDisrupGenomeRDSE];
fig = figure;
bar(Y,'group');
xlabel('operators');
ylabel('mean disruption');
set(gca,'XTickLabel', plotXTickLabels);
legend('Evolved Individuals','Random Individuals');
saveas(fig,sprintf('%sRobustez_%s',basedir, datestr(now, 30)),'epsc');
close(fig);
