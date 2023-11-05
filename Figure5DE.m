%% Load samples
load(['/home/labs/barkailab/bohdana/001_sequencing_data/004_medians/20231005_medians_BH_RM'], 'medianSumPromNewAll')

load(['/home/labs/barkailab/bohdana/001_sequencing_data/004_medians/20230813_sumprom_Gln3OL'], 'sumPromAllOL')

load(['/home/labs/barkailab/bohdana/001_sequencing_data/006_labTF/20230810_sumprom_allTF'], 'sumPromAllFelix', ...
        'sumPromAllOffir', 'sumPromAllHahn')
    
GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');
subtelomereGenes = GP.groups{1, 23}{1, 2}{1, 61};

fldNames = fieldnames(medianSumPromNewAll); 
for i = 1:length(fldNames)
    medianSumPromGln3.(fldNames{i}) = medianSumPromNewAll.(fldNames{i});
    medianSumPromGln3.(fldNames{i}) (subtelomereGenes) = nan;
end
strainNamesNewGln3 = fldNames;
clearvars fldNames i

fldNames = fieldnames(sumPromAllOL); 
for i = 1:length(fldNames)
    medianSumPromOL.(fldNames{i}) = sumPromAllOL.(fldNames{i});
    medianSumPromOL.(fldNames{i}) (subtelomereGenes) = nan;
end
strainNamesNewOL = fldNames;
clearvars fldNames i
 
currStrainsGln3 = [strainNamesNewOL; strainNamesNewGln3];
medianSumPromGln3.('Gln3_3del') = medianSumPromOL.('Gln3_3del');
medianSumPromGln3.('Gln3_7del') = medianSumPromOL.('Gln3_7del');

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

currStrainsFelix = ['GLN3_orig'; 'Dal82_Hahn'; 'Hap3_Hahn'; 'Hap5_Hahn'; strainNamesNewFelix];

medianSumPromFelix.('GLN3_orig') = medianSumPromGln3.('GLN3_orig');
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

clearvars medianSumPromOL sumPromAllOL sumPromAllOffir medianSumPromOffir strainNamesNewOffir strainNamesNewOL


%% Strains order for plots
currStrainsKO = {'GLN3_orig'; 'GLN3_norm_dal81'; 'GLN3_norm_dal82'; 'GLN3_rSTP1'; 'GLN3_rSTP2'; 
                            'GLN3_norm_hap2'; 'GLN3_rHAP3'; 'GLN3_rHAP5'; 'GLN3_rRTG1'; 'GLN3_rRTG3'; 'GLN3_norm_gat1'; 'Gln3_3del'; 
                            'GLN3_norm_gcn4'; 'Gln3_7del'};
                        
currStrainsTFs = {'GLN3_orig'; 'Dal81_IL'; 'Dal82_Hahn'; 'Stp1_TG'; 'Stp2_TG'; 
                  'Hap2_OL'; 'Hap3_Hahn'; 'Hap5_Hahn'; 'Rtg1_DK'; 'Rtg3_DK'; 
                  'Gat1_RM'; 'GLN3_orig'; 'Gcn4_TS'; 'GLN3_orig'};
                        
currStrainsHap2 = {'GLN3_orig'; 'GLN3_150_250C_WT'; 'GLN3_norm_hap2'; 'GLN3_d481_580_rHAP2'; 'GLN3_d503_505_rHAP2'};


%% Combine all medians in one table
for i = 1:length(currStrainsKO)
       medianSumPromKO(:,i) = medianSumPromGln3.(currStrainsKO{i}); 
end
clearvars i

for i = 1:length(currStrainsTFs)
       medianSumPromTF(:,i) = medianSumPromFelix.(currStrainsTFs{i}); 
end
clearvars i

for i = 1:length(currStrainsHap2)
       medianSumPromHap2(:,i) = medianSumPromGln3.(currStrainsHap2{i}); 
end
clearvars i


%% Correlation and labels
corr_toGln3KO = corr(medianSumPromKO, medianSumPromKO(:,1), 'rows', 'pairwise');
corr_toGln3TF = corr(medianSumPromTF, medianSumPromTF(:,1), 'rows', 'pairwise');
corr_toGln3Hap2 = corr(medianSumPromHap2, medianSumPromHap2(:,1), 'rows', 'pairwise');

addpath('/home/labs/barkailab/bohdana/001_sequencing_data/000_Scripts/functions')
color_scheme = brewermap(8, 'Greys');
color_spec = brewermap(11, 'Spectral');
cluster_colors=color_spec([8:11],:);

load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/targetFresh.mat')


%% FC fraction of binding: dot plot; correlation to Gln3 (Figure 5D)
currStrains = currStrainsKO;
medianSumProm = medianSumPromKO;
medianSumProm(isnan(medianSumProm)) = 0;
medianSumPromTF(isnan(medianSumPromTF)) = 0;

currStrainsLabels = strrep(currStrains, 'GLN3', '');
currStrainsLabels = strrep(currStrainsLabels, 'Gln3_', '');
currStrainsLabels = strrep(currStrainsLabels, 'norm_', '');
currStrainsLabels = strrep(currStrainsLabels, '_r', '');
currStrainsLabels = strrep(currStrainsLabels, '_', ' ');

for i = 1:4
    cluster = i;
    clusterProm = targets.geneId(targets.cluster==cluster);
    relativeBinding = median(log2(medianSumProm(clusterProm, :)+700) - log2(medianSumProm(clusterProm, 1)+700));
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

dotProdAll([12 14 15], :) = 50000;

xCoords = repmat([1:4],numel(currStrains),1);
yCoords = repmat([1:numel(currStrains)]',1,4);

figure('Color',[1 1 1],'Renderer','painters')
subplot(2, 3, [1 2 4 5])
scatter(xCoords(:), yCoords(:), rescale(dotProdAll(:),20,200, 'InputMin', 10000, 'InputMax', 100000), chgMat(:), 'filled', 'LineWidth', 1, 'MarkerEdgeColor',[1 1 1].*.4)
colorbar
colormap(gca,brighten(flipud(brewermap(1000,'RdBu')),0.3))
caxis([-2.5 2.5])
xticks([1:4])
yticks(1:length(currStrains))
yticklabels(currStrainsLabels)
set(gca,'ydir','reverse')
set(gcf,'renderer','painters')
xlim([0.5 4.5])
ylim([0.5 14.5])

subplot(2, 3, [3 6])
xAxis = categorical(currStrainsLabels);
xAxis = reordercats(xAxis, currStrainsLabels);

corr_toGln3 = corr(log2(medianSumProm(targets.geneId,:)+700), log2(medianSumProm(targets.geneId,1)+700), 'rows', 'pairwise');
scatter(xAxis, corr_toGln3, 80, 'MarkerFaceColor', cluster_colors(3,:), 'MarkerEdgeColor', [1 1 1].*.35)
xticklabels([1:14])
ylim([0 1])
grid on
view(90,90)
title("corr to Gln3 FL (only targets)")

saveas(gcf, '5d.svg');


%% Correlation matrix (Figure 5E)
currStrains = currStrainsHap2;
medianSumProm = medianSumPromHap2;
medianSumProm(isnan(medianSumProm)) = 0;
currStrainsLabels = strrep(currStrains, '_', ' ');

clusterProm = targets.geneId(targets.cluster==2);
relativeBinding = log2(medianSumProm(clusterProm, :)+700);
for i=1:5
    col(:,:,i) = relativeBinding - relativeBinding(:, i);
    relMat(i, :) = median(col(:, :, i));
    clearvars col
end
relMatR = round(relMat, 2);
clearvars clusterProm relativeBinding matrix

figure('Color',[1 1 1],'Renderer','painters')
imagesc(corr(log2(medianSumProm(targets.geneId,:)+700), log2(medianSumProm(targets.geneId,:)+700), 'rows', 'pairwise'))
caxis([0.4 1]);
yticks(1:length(currStrains))
xticks(1:length(currStrains))
yticklabels(currStrainsLabels)
xticklabels(currStrainsLabels)
xtickangle(45)
colormap(brewermap(1000, 'Blues'))
colorbar
hold on
plot([1:width(medianSumProm)]+[.5;.5], ylim()','-','Color',[.7 .7 .7])
plot(xlim()', [1:width(medianSumProm)]+[.5;.5],'-','Color',[.7 .7 .7])
plot(ylim, xlim,  'color', [1 1 1], 'LineWidth', 1, 'HandleVisibility', 'off')
saveas(gcf, '5e1.svg');


figure('Color',[1 1 1],'Renderer','painters')
imagesc(relMatR)
caxis([-3 0]);
yticks(1:length(currStrains))
xticks(1:length(currStrains))
yticklabels(currStrainsLabels)
xticklabels(currStrainsLabels)
xtickangle(45)
colormap(brewermap(1000, 'Greens'))
colorbar
hold on
plot([1:width(relMat)]+[.5;.5], ylim()','-','Color',[.7 .7 .7])
plot(xlim()', [1:width(relMat)]+[.5;.5],'-','Color',[.7 .7 .7])
plot(ylim, xlim,  'color', [1 1 1], 'LineWidth', 1, 'HandleVisibility', 'off')
saveas(gcf, 'Fig5E.svg');

