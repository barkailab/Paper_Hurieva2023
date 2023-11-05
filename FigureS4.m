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
currStrainsmod = {'GLN3_orig'; 'GLN3_bHAP2_pCglab';
                  'hap2'; 'HAP2_pCglab'};


%% Combine profiles in one table
for i = 1:length(currStrainsmod)
       medSumPromHap2(:,i) = medianSumPromGln3.(currStrainsmod{i}); 
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
    scatter(log2(medSumPromHap2(targets.geneId(ismember(targets.cluster,i)),1)+700),log2(medSumPromHap2(targets.geneId(ismember(targets.cluster,i)),2)+700), 80, cluster_colors(targets.cluster(ismember(targets.cluster,i)),:),'filled','MarkerEdgeColor',[1 1 1].*.25)
    hold on
end
xlabel("Gln3 wt")
ylabel("Gln3 with Hap2 from Cgla")
plot([10 19], [10 19], '--', 'color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off')
axis tight
idxReal=targets.geneId;
text(log2(medSumPromHap2(targets.geneId,1)+700),log2(medSumPromHap2(targets.geneId,2)+700), GP.gene_infoR64.name(idxReal), 'FontSize', 8)
saveas(gcf, 'FigS4_1.svg');


%% Correlation of DBD+R2 to Gln3 FL (Figure 6C)
figure('Color',[1 1 1], 'Renderer','painters')
for i=1:4
    scatter(log2(medSumPromHap2(targets.geneId(ismember(targets.cluster,i)),3)+700),log2(medSumPromHap2(targets.geneId(ismember(targets.cluster,i)),4)+700), 70, cluster_colors(targets.cluster(ismember(targets.cluster,i)),:),'filled','MarkerEdgeColor',[1 1 1].*.25)
    hold on
end
xlabel("Hap2 wt")
ylabel("Hap2 from Cgla")
plot([10 19], [10 19], '--', 'color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off')
axis tight
ylim([10 17])
idxReal=targets.geneId;
text(log2(medSumPromHap2(targets.geneId,3)+700),log2(medSumPromHap2(targets.geneId,4)+700), GP.gene_infoR64.name(idxReal), 'FontSize', 8)
saveas(gcf, 'FigS4_2.svg');
