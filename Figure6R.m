%% Load samples
load('20231005_medians_BH_RM.mat','medianSumPromNewAll')
cluster_colors=[0.913725490196078,0.596078431372549,0.478431372549020;0.400000000000000,0.760784313725490,0.647058823529412;0.196078431372549,0.533333333333333,0.741176470588235;0.721568627450980,0.670588235294118,0.827450980392157];
load('targetFresh.mat')
GP=load('group_imp.mat');
%% Strains order for plots and combine medians
currStrainsmod = {'GLN3_orig'; 'GLN3_norm_hap2';  'GLN3_bHAP2_pCglab';
                  'GLN3_bHAP2_d2_63'; 'GLN3_bHAP2_d78_131'; 'GLN3_bHAP2_d222_266'; 
                  'hap2'; 'HAP2_rGLN3'; 'HAP2_pCglab'; 
                  'HAP2_d2_63'; 'HAP2_d78_131'; 'HAP2_d222_266';
                  'GLN3_250N150C_d381_480'; 'GLN3_250N350C_WT'};




%% Combine profiles in one table
for i = 1:length(currStrainsmod)
       medSumPromHap2(:,i) = medianSumPromNewAll.(currStrainsmod{i}); 
end
clearvars i



figure('Color',[1 1 1], 'Renderer','painters')
%% Correlation of DBD+R2 to DBD (Figure 6B)
subplot(1,2,1)
for i=1:4
    scatter(log2(medSumPromHap2(targets.geneId(ismember(targets.cluster,i)),14)+700),log2(medSumPromHap2(targets.geneId(ismember(targets.cluster,i)),13)+700), 80, cluster_colors(targets.cluster(ismember(targets.cluster,i)),:),'filled','MarkerEdgeColor',[1 1 1].*.25)
    hold on
end
xlabel("Gln3 DBD")
ylabel("Gln3 DBD+R2")
plot([10 19], [10 19], '--', 'color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off')
axis tight
idxReal=targets.geneId;
text(log2(medSumPromHap2(targets.geneId,14)+700),log2(medSumPromHap2(targets.geneId,13)+700), GP.gene_infoR64.name(idxReal), 'FontSize', 8)
saveas(gcf, 'Fig6B.svg');


%% Correlation of DBD+R2 to Gln3 FL (Figure 6C)
subplot(1,2,2)
for i=1:4
    scatter(log2(medSumPromHap2(targets.geneId(ismember(targets.cluster,i)),1)+700),log2(medSumPromHap2(targets.geneId(ismember(targets.cluster,i)),13)+700), 70, cluster_colors(targets.cluster(ismember(targets.cluster,i)),:),'filled','MarkerEdgeColor',[1 1 1].*.25)
    hold on
end
xlabel("Gln3 FL")
ylabel("Gln3 DBD+R2")
plot([10 19], [10 19], '--', 'color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off')
axis tight
ylim([10 17])
idxReal=targets.geneId;
text(log2(medSumPromHap2(targets.geneId,1)+700),log2(medSumPromHap2(targets.geneId,13)+700), GP.gene_infoR64.name(idxReal), 'FontSize', 8)
saveas(gcf, 'Fig6C.svg');


%% Hap2 IDRs deletion (Figure 6E)
currStrainsTFs = {'hap2';'GLN3_orig'};
for i = 1:length(currStrainsTFs)
    medianSumPromTF(:,i) = medianSumPromNewAll.(currStrainsTFs{i});
end
clearvars i

sumPromColl = {medSumPromHap2(:, [9 10 11 12]);medSumPromHap2(:, [3 4 5 6])}
currStrainColl={currStrainsmod([9, 10, 11, 12]);currStrainsmod([3, 4, 5, 6])}
nameColl={'hap2','Gln3'}
figure('Color',[1 1 1],'Renderer','painters')
for g=1:2
    medianSumPromTarg = sumPromColl{g};
    currStrainsTarg = currStrainColl{g};


    currStrainsLabels = strrep(currStrainsTarg, '_', ' ');

    for i = 1:4
        cluster = i;
        clusterProm = targets.geneId(targets.cluster==cluster);
        relativeBinding = median(log2(medianSumPromTarg(clusterProm, :)+700) - log2(medianSumPromTarg(clusterProm, 1)+700), 'omitnan');
        chgMat(:, i) = relativeBinding';
        clearvars cluster clusterProm fracBindingCluster fc relativeBinding
    end
    clearvars i

    for i = 1:4
        cluster = i;
        clusterProm = targets.geneId(targets.cluster==cluster);
        dotProdCluster = medianSumPromTF(clusterProm, g)'*medianSumPromTF(clusterProm, g)./ sum(medianSumPromTF(clusterProm, g),  'omitnan');
        dotProdAll(:, i) = dotProdCluster.*ones(numel(currStrainsTarg),1);
        clearvars cluster clusterProm dotProdCluster sumBindingCluster
    end
    clearvars i

    xCoords = repmat([1:4],numel(currStrainsTarg),1);
    yCoords = repmat([1:numel(currStrainsTarg)]',1,4);

    subplot(1,2,g)
    hold off
    scatter(xCoords(:), yCoords(:), rescale(dotProdAll(:),50,300, 'InputMin', 10000, 'InputMax', 100000), chgMat(:), 'filled', 'LineWidth', 1, 'MarkerEdgeColor',[1 1 1].*.4)
    colorbar
    colormap(gca,brighten(flipud(brewermap(1000,'RdBu')),0.3))
    caxis([-2.5 2.5])
    xticks([1:4])
    yticks(1:length(currStrainsTarg))
    yticklabels(currStrainsLabels)
    set(gca,'ydir','reverse')
    set(gcf,'renderer','painters')
    xlim([0.5 4.5])
    ylim([0.5 length(currStrainsTarg)+.5])
    title(nameColl{g})
end
saveas(gcf, 'Fig6E.svg');


%% Correlation Hap2 to Hap2 without Gln3 (Figure 6F, left)
figure('Color',[1 1 1], 'Renderer','painters')
hold off
for i=1:4
    scatter(log2(medSumPromHap2(targets.geneId(ismember(targets.cluster,i)),7)+700),log2(medSumPromHap2(targets.geneId(ismember(targets.cluster,i)),8)+700), 80, cluster_colors(targets.cluster(ismember(targets.cluster,i)),:),'filled','MarkerEdgeColor',[1 1 1].*.25)
    hold on
end
xlabel("Hap2 wt")
ylabel("Hap2 rGln3")
plot([10 19], [10 19], '--', 'color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off')
axis tight
idxReal=targets.geneId;
%text(log2(medSumPromHap2(targets.geneId,7)+700),log2(medSumPromHap2(targets.geneId,8)+700), GP.gene_infoR64.name(idxReal), 'FontSize', 8)
saveas(gcf, 'Fig6F1.svg');


%% Gln3/Hap2 effect (Figure 6F, right)
figure('Color',[1 1 1], 'Renderer','painters')
hold off
gln3_effect = log2(medSumPromHap2(targets.geneId, 1)+700) -  log2(medSumPromHap2(targets.geneId, 2)+700);
hap2_effect = log2(medSumPromHap2(targets.geneId, 7)+700) -  log2(medSumPromHap2(targets.geneId, 8)+700);
for i=1:4
    scatter(gln3_effect(targets.cluster==i), hap2_effect(targets.cluster==i), 80, cluster_colors(targets.cluster(ismember(targets.cluster,i)),:),'filled','MarkerEdgeColor',[1 1 1].*.35)
    hold on
end
xlabel("log2(Gln3 wt) - log2(Gln3 rHap2)")
ylabel("log2(Hap2 wt) - log2(Hap2 rGln3)")
axis tight
saveas(gcf, 'Fig6F2.svg');


%% Figure S5AB
clearvars
load('20231005_medians_BH_RM.mat','medianSumPromNewAll')
load('targetFresh.mat')
cmpStrains={'GLN3_orig','GLN3_bHAP2_pCglab';'hap2','HAP2_pCglab'}
cluster_colors=[0.913725490196078,0.596078431372549,0.478431372549020;0.400000000000000,0.760784313725490,0.647058823529412;0.196078431372549,0.533333333333333,0.741176470588235;0.721568627450980,0.670588235294118,0.827450980392157];

for g=1:2
    subplot(1,2,g)
    for i=1:4
        scatter(log2(medianSumPromNewAll.(cmpStrains{g,1})(targets.geneId(ismember(targets.cluster,i)))+700),log2(medianSumPromNewAll.(cmpStrains{g,2})(targets.geneId(ismember(targets.cluster,i)))+700), 80, cluster_colors(targets.cluster(ismember(targets.cluster,i)),:),'filled','MarkerEdgeColor',[1 1 1].*.25)
        hold on
    end
    xlabel("Gln3 DBD")
    ylabel("Gln3 DBD+R2")
    plot([10 19], [10 19], '--', 'color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off')
    axis tight
end



%% Figure S5C
clear all
load("nBS.mat") %  Hap2 and GLn3 binding sites in every promoter (6701)
load('targetFresh.mat')
load('promoterLengthsORF')
promCluster=full(sparse(targets.geneId,1,targets.newCluster,6701,1));
promCluster(isnan(promoterLengthsORF))=nan;

for i=1:2
    subplot(1,2,i)
    hold off
    selProms=promCluster==0;
    histB=histcounts(nBS(selProms,i),[-.5:1:10.5],'Normalization','pdf')
    plot([0:10],histB)
    
    selProms=promCluster==2;
    histB=histcounts(nBS(selProms,i),[-.5:1:10.5],'Normalization','pdf')
    hold on
    plot([0:10],histB)
    
        selProms=ismember(promCluster,[1,3,4]);
    histB=histcounts(nBS(selProms,i),[-.5:1:10.5],'Normalization','pdf')
    hold on
    plot([0:10],histB)
    xlim([0 7])
    ylim([0 1])
end
saveas(gcf, 'FigS5.svg');