addpath('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem')

goEnr = goAnalysis(targets.geneId,'BG',1:6701,'type',goType);
goEnr.Var4 = categorical(goEnr.Var4);
idx = goEnr.Var4 == 'C';
goEnrNew = goEnr(~idx, :);

figure('Color',[1 1 1], 'Renderer','painters')
scatter(goEnrNew.n, -log10(goEnrNew.pVals), 30, goEnrNew.Var4, 'filled','MarkerEdgeColor',[1 1 1].*.35)
xlabel("Number of Targets")
ylabel("p-value")

labels = goEnrNew.Var5;
text(goEnrNew.n, -log10(goEnrNew.pVals), labels)
