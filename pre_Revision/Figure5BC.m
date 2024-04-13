%% Load datasets: all TFs (x3) and Gln3
load(['/home/labs/barkailab/bohdana/001_sequencing_data/006_labTF/20230810_sumprom_allTF'], 'sumPromAllFelix', ...
        'sumPromAllOffir', 'sumPromAllHahn')

load(['/home/labs/barkailab/bohdana/001_sequencing_data/004_medians/20230704_medians_BH_RM'], 'medianSumPromNewAll')
     
GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');
subtelomereGenes = GP.groups{1, 23}{1, 2}{1, 61};

fldName = fieldnames(sumPromAllFelix); %get field names from the struct
shortMatr = rmfield(sumPromAllFelix, fldName(find(contains(fldName, {'del', 'DBD', 'Sfp1_OL'}))));
fldNames = fieldnames(shortMatr); 
[~, firstAppearance, ~]=unique(extractBefore(fldNames,'_')); %for Felix's
for i = 1:length(firstAppearance)
    medianSumPromFelix.(fldNames{firstAppearance(i)}) = sumPromAllFelix.(fldNames{firstAppearance(i)});
    medianSumPromFelix.(fldNames{firstAppearance(i)}) (subtelomereGenes) = nan;
end
for i = 1:length(firstAppearance)
    strainNamesNewFelix(i, 1) = cellstr(fldNames{firstAppearance(i), 1});
end
clearvars fldName fldNames i shortMatr firstAppearance

fldNames = fieldnames(sumPromAllOffir);
for i = 1:length(fldNames)
    medianSumPromOffir.(fldNames{i}) = sumPromAllOffir.(fldNames{i});
    medianSumPromOffir.(fldNames{i}) (subtelomereGenes) = nan;
end
strainNamesNewOffir = fldNames;
clearvars fldNames i

fldName = fieldnames(sumPromAllHahn); 
shortMatr = rmfield(sumPromAllHahn, fldName(find(contains(fldName, '_'))));
fldNames = fieldnames(shortMatr); 
for i = 1:length(fldNames)
    medianSumPromHahn.(fldNames{i}) = sumPromAllHahn.(fldNames{i});
    medianSumPromHahn.(fldNames{i}) (subtelomereGenes) = nan;
end
strainNamesNewHahn = fldNames;
clearvars fldName fldNames i shortMatr

fldNames = fieldnames(medianSumPromNewAll); 
for i = 1:length(fldNames)
    medianSumPromGln3.(fldNames{i}) = medianSumPromNewAll.(fldNames{i});
    medianSumPromGln3.(fldNames{i}) (subtelomereGenes) = nan;
end
strainNamesNewGln3 = fldNames;
clearvars fldNames i


%% Combine all medians in one table
currStrainsFelix = ['GLN3_orig'; 'Gln3DBD_wt'; 'Hap4_BH'; 'Dal82_Hahn'; 'Hap3_Hahn'; 'Hap5_Hahn'; strainNamesNewFelix];

medianSumPromFelix.('GLN3_orig') = medianSumPromGln3.('GLN3_orig');
medianSumPromFelix.('Gln3DBD_wt') = medianSumPromGln3.('GLN3_250N350C_WT');
medianSumPromFelix.('Hap4_BH') = medianSumPromGln3.('hap4');
medianSumPromFelix.('Gcn4_TS') = medianSumPromOffir.('Gcn4');
medianSumPromFelix.('Yap1_DK') = medianSumPromOffir.('Yap1');
medianSumPromFelix.('Yhp1_TG') = medianSumPromOffir.('Yhp1');
medianSumPromFelix.('Dal82_Hahn') = medianSumPromHahn.('DAL82');
medianSumPromFelix.('Hap3_Hahn') = medianSumPromHahn.('HAP3');
medianSumPromFelix.('Hap5_Hahn') = medianSumPromHahn.('HAP5');

for i = 1:length(currStrainsFelix)
       sumPromFelix(:,i) = medianSumPromFelix.(currStrainsFelix{i}); 
end
clearvars i

currStrainsLabelsFelix = extractBefore(currStrainsFelix,'_');


%% Correlation and labels 
corr_toGln3Felix= corr(sumPromFelix, sumPromFelix(:,1), 'rows', 'pairwise');
corr_toGln3dbd= corr(sumPromFelix, sumPromFelix(:,2), 'rows', 'pairwise');

addpath('/home/labs/barkailab/bohdana/001_sequencing_data/000_Scripts/functions')
color_scheme = brewermap(8, 'Greys');
color_spec = brewermap(11, 'Spectral');
cluster_colors=color_spec([8:11],:);

load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/targetFresh.mat')
load('/home/labs/barkailab/divyakr/matlab/checProfileCreationFiles/GP.mat');


%% Correlation to Gln3 / DBD (Figure 5B)
figure('Color',[1 1 1],'Renderer','painters')
scatter(corr_toGln3Felix, corr_toGln3dbd, [], color_scheme(5, :), 'filled', 'MarkerEdgeColor',[1 1 1].*.35)
xlabel('corr to Gln3')
ylabel('corr to Gln3 DBD')
hold on
xlim([0 1])
ylim([0 1])
plot(xlim, ylim, '--', 'color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off')
[~, k] = maxk(corr_toGln3Felix, 7);
[~, n] = maxk(corr_toGln3dbd, 4);
list = unique([k; n]);
tfs = unique([currStrainsLabelsFelix(list); 'Stp1'; 'Stp2'; 'Dal81'; 'Hap2'; 'Hap4'; 'Gat1'; 'Gzf3'; 'Dal80'; 'Gcn4'; 'Rtg1'; 'Dal82'; 'Hap3'; 'Hap5']);
Y = find(contains(currStrainsLabelsFelix, tfs, 'ignorecase', true));
hold on
scatter(corr_toGln3Felix(Y), corr_toGln3dbd(Y), [], cluster_colors(3,:),'filled','MarkerEdgeColor',[1 1 1].*.35)
text(corr_toGln3Felix(Y), corr_toGln3dbd(Y), currStrainsLabelsFelix(Y))
axis square

saveas(gcf, 'fig5_corrGln3DBD3.svg');


%% Dot product heatmap (Figure 5C)
corr_toGln3 = corr_toGln3Felix;
currStrainsLabels = currStrainsLabelsFelix;
currStrains = currStrainsFelix;
sumProm = sumPromFelix;
sumProm(isnan(sumProm)) = 0;

gln3_pos = find(contains(currStrains, 'Gln3', 'IgnoreCase', true));
currStrainsLabelsFilt = currStrainsLabels(find(~contains(currStrains, 'Gln3', 'IgnoreCase', true)));
corr_toGln3Filt = corr_toGln3(find(~contains(currStrains, 'Gln3', 'IgnoreCase', true)));
[val, idx] = sort(corr_toGln3Filt, 'descend');
idx = idx(~isnan(val));

for i = 1:4
    cluster = i;
    clusterProm = targets.geneId(targets.cluster==cluster);
    dotProdCluster = sumProm(clusterProm, 1)'*sumProm(clusterProm, :)./sum(sumProm(clusterProm, 1));
    dotProdAll(i, :) = dotProdCluster;
    clearvars cluster clusterProm dotProdCluster
end
clearvars i
dotProdAllFilt = dotProdAll(:, ~ismember(1:175, gln3_pos));

figure('Color',[1 1 1],'Renderer','painters')
imagesc(dotProdAllFilt(:, idx))
colorbar
colormap(brewermap(1000, 'Blues'))

%[~, k] = maxk(max(dotProdAllFilt, [], 2), 10);
[~, k] = maxk(max(dotProdAllFilt), 10);
tfs =  unique([currStrainsLabelsFilt(k); 'Stp1'; 'Stp2'; 'Dal81'; 'Hap2'; 'Hap4'; 'Gat1'; 'Gzf3'; 'Dal80'; 'Gcn4'; 'Rtg1'; 'Dal82'; 'Hap3'; 'Hap5']);
Y = find(contains(currStrainsLabelsFilt(idx), tfs, 'ignorecase', true));

xticks(Y)
new_label = currStrainsLabelsFilt(idx);
xticklabels(new_label(Y))
xtickangle(90)
yticks([1:4])
caxis([10000 100000])
ylabel("Gln3 target clusters")
clearvars tfs Y k

saveas(gcf, 'fig5_clusterHM.svg');

