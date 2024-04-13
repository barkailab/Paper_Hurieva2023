%% Load samples
load(['/home/labs/barkailab/bohdana/001_sequencing_data/004_medians/20231005_medians_BH_RM'], 'medianSumPromNewAll')
    
GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');
subtelomereGenes = GP.groups{1, 23}{1, 2}{1, 61};

fldNames = fieldnames(medianSumPromNewAll); 
for i = 1:length(fldNames)
    medianSumPromGln3.(fldNames{i}) = medianSumPromNewAll.(fldNames{i});
    medianSumPromGln3.(fldNames{i}) (subtelomereGenes) = nan;
end
strainNamesNewGln3 = fldNames;
clearvars fldNames i


%% Strains order for plots and combine medians
currStrainsmod = {'GLN3_orig'; 'GLN3_norm_hap2';  'GLN3_bHAP2_pCglab';
                  'GLN3_bHAP2_d2_63'; 'GLN3_bHAP2_d78_131'; 'GLN3_bHAP2_d222_266'; 
                  'hap2'; 'HAP2_rGLN3'; 'HAP2_pCglab'; 
                  'HAP2_d2_63'; 'HAP2_d78_131'; 'HAP2_d222_266';
                  'GLN3_250N150C_d381_480'; 'GLN3_250N350C_WT'};

currStrainsTFs = {'hap2'; 'hap2'; 'hap2'; 'hap2'}; 

currStrainsTFs = {'GLN3_orig'; 'GLN3_orig'; 'GLN3_orig';  'GLN3_orig'}; 


%% Combine profiles in one table
for i = 1:length(currStrainsmod)
       medSumPromHap2(:,i) = medianSumPromGln3.(currStrainsmod{i}); 
end
clearvars i

for i = 1:length(currStrainsTFs)
       medianSumPromTF(:,i) = medianSumPromGln3.(currStrainsTFs{i}); 
end
clearvars i


%% Correlation and labels
currStrainsLabels = strrep(currStrainsmod, '_', ' ');

addpath('/home/labs/barkailab/bohdana/001_sequencing_data/000_Scripts/functions')
color_scheme = brewermap(8, 'Greys');
color_spec = brewermap(11, 'Spectral');
cluster_colors=color_spec([8:11],:);

load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/targetFresh.mat')


%% Correlation of DBD+R2 to DBD (Figure 6B)
figure('Color',[1 1 1], 'Renderer','painters')
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
figure('Color',[1 1 1], 'Renderer','painters')
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
medianSumPromTarg = medSumPromHap2(:, [9 10 11 12]);
medianSumPromTarg = medSumPromHap2(:, [3 4 5 6]);

medianSumPromTarg(isnan(medianSumPromTarg)) = 0;
medianSumPromTF(isnan(medianSumPromTF)) = 0;

currStrainsTarg = currStrainsmod([9, 10, 11, 12]);
currStrainsTarg = currStrainsmod([3, 4, 5, 6]);

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
    dotProdCluster = medianSumPromTF(clusterProm, 1)'*medianSumPromTF(clusterProm, :)./ sum(medianSumPromTF(clusterProm, 1),  'omitnan');
    dotProdAll(:, i) = dotProdCluster';
    clearvars cluster clusterProm dotProdCluster sumBindingCluster
end
clearvars i

xCoords = repmat([1:4],numel(currStrainsTarg),1);
yCoords = repmat([1:numel(currStrainsTarg)]',1,4);

figure('Color',[1 1 1],'Renderer','painters')
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
saveas(gcf, 'Fig6E.svg');


%% Correlation Hap2 to Hap2 without Gln3 (Figure 6F, left)
figure('Color',[1 1 1], 'Renderer','painters')
for i=1:4
    scatter(log2(medSumPromHap2(targets.geneId(ismember(targets.cluster,i)),7)+700),log2(medSumPromHap2(targets.geneId(ismember(targets.cluster,i)),8)+700), 80, cluster_colors(targets.cluster(ismember(targets.cluster,i)),:),'filled','MarkerEdgeColor',[1 1 1].*.25)
    hold on
end
xlabel("Hap2 wt")
ylabel("Hap2 rGln3")
plot([10 19], [10 19], '--', 'color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off')
axis tight
idxReal=targets.geneId;
text(log2(medSumPromHap2(targets.geneId,7)+700),log2(medSumPromHap2(targets.geneId,8)+700), GP.gene_infoR64.name(idxReal), 'FontSize', 8)
saveas(gcf, 'Fig6F1.svg');


%% Gln3/Hap2 effect (Figure 6F, right)
figure('Color',[1 1 1], 'Renderer','painters')
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

