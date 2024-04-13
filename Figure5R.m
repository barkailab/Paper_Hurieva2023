%% Figure 5A
clear all
load('allTFs.mat')
GP=load('group_imp.mat');
load('targetFresh.mat')
load('promType')

corr_toGln3Felix= corr(sumPromFelix, sumPromFelix(:,1), 'rows', 'pairwise');
corr_toGln3dbd= corr(sumPromFelix, sumPromFelix(:,2), 'rows', 'pairwise');
currStrainsLabelsFelix = extractBefore(currStrainsFelix,'_');

color_scheme = brewermap(8, 'Greys');
color_spec = brewermap(11, 'Spectral');
cluster_colors=color_spec([8:11],:);



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


%% Dot product heatmap (Figure 5C) + S4
corr_toGln3 = corr_toGln3Felix;
currStrainsLabels = currStrainsLabelsFelix;
currStrains = currStrainsFelix;
sumProm = sumPromFelix;
%sumProm(isnan(sumProm)) = 0;

%gln3_pos = find(contains(currStrains, 'Gln3', 'IgnoreCase', true));
%currStrainsLabelsFilt = currStrainsLabels(find(~contains(currStrains, 'Gln3', 'IgnoreCase', true)));
%corr_toGln3Filt = corr_toGln3(find(~contains(currStrains, 'Gln3', 'IgnoreCase', true)));
%[val, idx] = sort(corr_toGln3Filt, 'descend');
%idx = idx(~isnan(val));

for i = 1:4
    cluster = i;
    clusterProm = targets.geneId(targets.cluster==cluster);
    dotProdCluster = (sumProm(clusterProm, 1)'*sumProm(clusterProm, :))./sum(sumProm(clusterProm, 1));
    dotProdAll(i, :) = dotProdCluster;
    clearvars cluster clusterProm dotProdCluster
end
clearvars i
%dotProdAllFilt = dotProdAll(:, ~ismember(1:width(dotProdAll), gln3_pos));
[val,idx]=sort(corr_toGln3,'descend')
idx(isnan(val))=[];
idx(ismember(idx,find(contains(currStrains, 'Gln3', 'IgnoreCase', true))))=[];
figure('Color',[1 1 1],'Renderer','painters')
imagesc(dotProdAll(:, idx))
colorbar
colormap(brewermap(1000, 'Blues'))
%[~, k] = maxk(max(dotProdAllFilt, [], 2), 10);
[~, k] = maxk(max(dotProdAll), 10);
tfs =  unique([currStrainsLabels(k); 'Stp1'; 'Stp2'; 'Dal81'; 'Hap2'; 'Hap4'; 'Gat1'; 'Gzf3'; 'Dal80'; 'Gcn4'; 'Rtg1'; 'Dal82'; 'Hap3'; 'Hap5']);
Y = find(contains(currStrainsLabels(idx), tfs, 'ignorecase', true));

xticks(Y)
new_label = currStrainsLabels(idx);
xticklabels(new_label(Y))
xtickangle(90)
yticks([1:4])
caxis([10000 100000])
ylabel("Gln3 target clusters")
clearvars tfs Y k

saveas(gcf, 'fig5_clusterHM.svg');
xlim([0 30]+.5)
saveas(gcf, 'fig5_clusterHMzoom.svg');


%% Figure 5D
clearvars
load('allTFs.mat')
load('targetFresh.mat')
load('20231005_medians_BH_RM.mat','medianSumPromNewAll')

currStrainsLabelsFelix = extractBefore(currStrainsFelix,'_');

currStrainsKO = {'GLN3_orig'; 'GLN3_norm_dal81'; 'GLN3_norm_dal82'; 'GLN3_rSTP1'; 'GLN3_rSTP2'; 
                            'GLN3_norm_hap2'; 'GLN3_rHAP3'; 'GLN3_rHAP5'; 'GLN3_rRTG1'; 'GLN3_rRTG3'; 'GLN3_norm_gat1'; 'Gln3_3del'; 
                            'GLN3_norm_gcn4'; 'Gln3_7del'};
                        
currStrainsTFs = {'GLN3_orig'; 'Dal81_IL'; 'Dal82_Hahn'; 'Stp1_TG'; 'Stp2_TG'; 
                  'Hap2_OL'; 'Hap3_Hahn'; 'Hap5_Hahn'; 'Rtg1_DK'; 'Rtg3_DK'; 
                  'Gat1_RM'; 'GLN3_orig'; 'Gcn4_TS'; 'GLN3_orig'};
                       
%% Combine all medians in one table
for i = 1:length(currStrainsKO)
       medianSumPromKO(:,i) = medianSumPromNewAll.(currStrainsKO{i}); 
end
clearvars i

for i = 1:length(currStrainsTFs)
       medianSumPromTF(:,i) = sumPromFelix(:,ismember(currStrainsFelix,currStrainsTFs{i}));
end
clearvars i



%% Correlation and labels
corr_toGln3KO = corr(medianSumPromKO(promType<3,:), medianSumPromKO(promType<3,1), 'rows', 'pairwise');
corr_toGln3TF = corr(medianSumPromTF(promType<3,:), medianSumPromTF(promType<3,1), 'rows', 'pairwise');

color_scheme = brewermap(8, 'Greys');
color_spec = brewermap(11, 'Spectral');
cluster_colors=color_spec([8:11],:);


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

dotProdAll([12 14], :) = 50000;

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
clearvars -except medianSumPromNewAll targets
currStrainsHap2 = {'GLN3_orig'; 'GLN3_150_250C_WT'; 'GLN3_norm_hap2'; 'GLN3_d481_580_rHAP2'; 'GLN3_d503_505_rHAP2'};
for i = 1:length(currStrainsHap2)
       medianSumPromHap2(:,i) = medianSumPromNewAll.(currStrainsHap2{i}); 
end
clearvars i
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
subplot(1,2,1)
imagesc(corr(log2(medianSumProm(targets.geneId,:)+700), log2(medianSumProm(targets.geneId,:)+700), 'rows', 'pairwise'))
caxis([0.4 1]);
yticks(1:length(currStrains))
xticks(1:length(currStrains))
yticklabels(currStrainsLabels)
xticklabels(currStrainsLabels)
xtickangle(45)
colormap(gca,brewermap(1000, 'Blues'))
colorbar
hold on
plot([1:width(medianSumProm)]+[.5;.5], ylim()','-','Color',[.7 .7 .7])
plot(xlim()', [1:width(medianSumProm)]+[.5;.5],'-','Color',[.7 .7 .7])
plot(ylim, xlim,  'color', [1 1 1], 'LineWidth', 1, 'HandleVisibility', 'off')
subplot(1,2,2)
imagesc(relMatR)
caxis([-3 0]);
yticks(1:length(currStrains))
xticks(1:length(currStrains))
yticklabels(currStrainsLabels)
xticklabels(currStrainsLabels)
xtickangle(45)
colormap(gca,brewermap(1000, 'Greens'))
colorbar
hold on
plot([1:width(relMat)]+[.5;.5], ylim()','-','Color',[.7 .7 .7])
plot(xlim()', [1:width(relMat)]+[.5;.5],'-','Color',[.7 .7 .7])
plot(ylim, xlim,  'color', [1 1 1], 'LineWidth', 1, 'HandleVisibility', 'off')

saveas(gcf, 'Fig5E.svg');
