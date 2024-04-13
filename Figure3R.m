%% Colors and targets 
clear all
load('targetFresh.mat')
GP=load('group_imp.mat');
cluster_colors=[0.913725490196078,0.596078431372549,0.478431372549020;0.400000000000000,0.760784313725490,0.647058823529412;0.196078431372549,0.533333333333333,0.741176470588235;0.721568627450980,0.670588235294118,0.827450980392157];


%% Conditions (Figure 3B)
clearvars -except GP targets cluster_colors

repeatData=load('repeatData.mat')
currStrains={'GLN3_norm_ure2_Am','GLN3_norm_ure2_NG','GLN3_norm_ure2_Pr','GLN3_norm_ure2_sd15','GLN3_norm_ure2_sd25',...
    'GLN3_norm_wt_Am','GLN3_norm_ure2_NG','GLN3_norm_ure2_Pr','GLN3_norm_ure2_sd15','GLN3_norm_ure2_sd25'}

[~,repeatData.samples.strainOrder]=ismember(repeatData.samples.strain,currStrains)
repeats=repeatData.samples(repeatData.samples.strainOrder>0,:)
sumProm=repeatData.sumPromRep(:,repeatData.samples.strainOrder>0)
[repeats,idx]=sortrows(repeats,"strainOrder")
sumProm=sumProm(:,idx);
clearvars idx 


figure('Color',[1 1 1],'Renderer','painters','Position',[1905 554 1533 301])
subplot(1,2,1)
bL=find(diff(repeats.strainOrder)>0);
bLmin=[1:height(repeats)]'

imagesc(corr(log2(sumProm(targets.geneId,:)+700)))
set(gca,'Ytick',1:height(repeats),'YTickLabel',repeats.name)
hold on
plot(bL'+0.5.*[1;1],ylim','-','LineWidth',2,'Color',[1 1 1].*0.85)
plot(xlim',bL'+0.5.*[1;1],'-','LineWidth',2,'Color',[1 1 1].*0.85)
plot(bLmin'+0.5.*[1;1],ylim','-','LineWidth',1,'Color',[1 1 1].*0.85)
plot(xlim',bLmin'+0.5.*[1;1],'-','LineWidth',1,'Color',[1 1 1].*0.85)
caxis([0.8 1])
title('good target genes -log2')
ylabel(colorbar(),'correlation across targets')
colormap(gca,brighten((brewermap(128,'blues')),0.1))

relativeBinding = log2(sumProm+700)-log2(mean(sumProm(:,contains(repeats.strain,'GLN3_norm_ure2_sd')),2)+700);
pVal=nan(6701,max(repeats.strainOrder));
mainDir=[-1;-1;1;-1]
for i=1:max(repeats.strainOrder)
    for g=find(~isnan(relativeBinding(:,1)))'
        [~,pVal(g,i)]=ttest2((700+sumProm(g,repmat(find(repeats.strainOrder==i),1,1))),(700+sumProm(g,repmat(find(repeats.strainOrder>3),1,1))));
    end
%     subplot(1,6,i)
%     scatter(mean(relativeBinding(targets.geneId,repeats.strainOrder==i),2),-log10(pVal(:,i)),[],targets.newCluster,'filled')
%     caxis([.5 4.5])
%     colormap(gca,cluster_colors)
end
targets.pVal=max(-log10(pVal(targets.geneId,[1,3])),[],2).*sign(mean(relativeBinding(targets.geneId,[1,2,5,6]),2)).*mainDir(targets.newCluster)
targets=sortrows(targets,{'newCluster','pVal'},{'ascend','descend'})
subplot(1,10,6:8)
imagesc(relativeBinding(targets.geneId,:))
caxis([-1 1].*1.5)
set(gca,'Xtick',1:height(repeats),'XTickLabel',repeats.name)
hold on
plot(bLmin'+0.5.*[1;1],ylim','-','LineWidth',1,'Color',[1 1 1].*0.85)
plot(bL'+0.5.*[1;1],ylim','-','LineWidth',2,'Color',[1 1 1].*0.85)
gbL=find(diff(targets.newCluster));
plot(xlim',gbL'+[.5;.5],'k-','linewidth',1)
ylabel(colorbar(),'log2 change')
colormap(gca,brighten(flipud(brewermap(128,'RdBu')),0.3))

subplot(1,10,9:10)
imagesc(max(-log10(pVal(targets.geneId,[1,3])),[],2))
caxis([0 3])

colormap(gca,brighten(brewermap(128,'YlGn'),0.3))
colorbar()
saveas(gcf,'Figure3-Cond.svg')

%% Figure 3D
clearvars

load('20231005_medians_BH_RM.mat','medianSumPromNewAll')     
GP=load('group_imp.mat');
load('targetFresh.mat')

currStrainsPhospA = {'GLN3_250NA_wt'; 'GLN3_200NA_wt'; 'GLN3_150NA_wt'; 'GLN3_100NA_wt'; 'GLN3_50NA_wt';
                                    'GLN3_orig';
                                     'GLN3_50CA_wt'; 'GLN3_100CA_wt'; 'GLN3_150CA_wt'; 'GLN3_200CA_wt';
                                     'GLN3_300CA_wt'; 'GLN3_350CA_wt'; 'GLN3_375CA_wt';};
                                 
currStrainsPhospD = {'GLN3_250ND_wt'; 'GLN3_200ND_wt'; 'GLN3_150ND_wt'; 'GLN3_100ND_wt'; 'GLN3_50ND_wt'; 
                                    'GLN3_orig';
                                     'GLN3_50CD_wt'; 'GLN3_100CD_wt'; 'GLN3_150CD_wt'; 'GLN3_200CD_wt';
                                     'GLN3_300CD_wt'; 'GLN3_350CD_wt'; 'GLN3_375CD_wt'; };
                              
for i = 1:length(currStrainsPhospA)
       medianSumPromA(:,i) = medianSumPromNewAll.(currStrainsPhospA{i}); 
end

for i = 1:length(currStrainsPhospD)
       medianSumPromD(:,i) = medianSumPromNewAll.(currStrainsPhospD{i}); 
end

strainsLists={currStrainsPhospA,currStrainsPhospD}
sumProm={medianSumPromA,medianSumPromD}
mutNames={'ST to A','ST to D'}
figure('Color', [1 1 1], 'Renderer', 'painters')
for i=1:numel(strainsLists)
    subplot(1,2,i)
    hold off
    relativeBinding = log2(sumProm{i}+700)-log2(mean(sumProm{i}(:,ismember(strainsLists{i},'GLN3_orig')),2)+700);  
    imagesc(relativeBinding(targets.geneId,:),[-1 1].*4.5)
    colormap(gca,brighten(flipud(brewermap(1000,'RdBu')),0.3))
    colorbar()
    hold on
    plot(xlim',find(diff(targets.cluster)>0)'+[.5;.5],'-','linewidth',1,'Color',[1 1 1].*0.25)
    ylabel("Gln3 targets by truncation clusters (n=184)")
    title(mutNames{i}, 'FontWeight', 'normal',  'FontSize', 11)
end
saveas(gcf, 'Fig3d.svg')


%% Evolution Figure 3E,F
clearvars

load('20231005_medians_BH_RM.mat','medianSumPromNewAll')     
GP=load('group_imp.mat');
load('targetFresh.mat')
cluster_colors=[0.913725490196078,0.596078431372549,0.478431372549020;0.400000000000000,0.760784313725490,0.647058823529412;0.196078431372549,0.533333333333333,0.741176470588235;0.721568627450980,0.670588235294118,0.827450980392157];

currStrainsEvol = {'GLN3_orig'; 'GLN3_Smik_wt'; 'GLN3_Skud_wt'; 'GLN3_gla_wt'; 'GLN3_klu_wt'; 'GLN3_lac_wt';};
for i = 1:length(currStrainsEvol)
       medianSumPromEvol(:,i) = medianSumPromNewAll.(currStrainsEvol{i}); 
end
clearvars i

figure('Color',[1 1 1],'Renderer','painters')
imagesc(corr(log2(medianSumPromEvol(targets.geneId,:)+700), log2(medianSumPromEvol(targets.geneId,:)+700), 'rows', 'pairwise'))
hold on
ylabel(colorbar(),'sumProm corr')
set(gca,'Ytick',1:numel(currStrainsEvol),'Yticklabel',currStrainsEvol,'Xtick',1:numel(currStrainsEvol),'Xticklabel',currStrainsEvol)
title('sumprom correlation all targets')
colormap(gca,brighten(brewermap(128,'blues'),0.3))
plot([1.5:numel(currStrainsEvol)].*[1;1],ylim','-','Color',[1 1 1].*0.65)
plot(xlim',[1.5:numel(currStrainsEvol)].*[1;1],'-','Color',[1 1 1].*0.65)
caxis([.75 1])
saveas(gcf, 'fig3_evolCorr.svg');

figure('Color',[1 1 1],'Renderer','painters')
subplot(1,2,1)
scatter(log2(medianSumPromEvol(targets.geneId,1)+700), log2(medianSumPromEvol(targets.geneId,4)+700), 100, targets.cluster,'filled','MarkerEdgeColor',[1 1 1].*.35)
colormap(gca, cluster_colors)
xlabel("S. cerevisiae")
ylabel("C. glabrata")
hold on
plot([10 18.5], [10 18.5], '--', 'color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off')
axis tight
axis square
xlim([10 18.5])
ylim([10 18.5])

subplot(1,2,2)
scatter(log2(medianSumPromEvol(targets.geneId,1)+700), log2(medianSumPromEvol(targets.geneId,6)+700), 100, targets.cluster,'filled','MarkerEdgeColor',[1 1 1].*.35)
colormap(gca, cluster_colors)
xlabel("S. cerevisiae")
ylabel("K. lactis")
hold on
plot([10 18.5], [10 18.5], '--', 'color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off')
axis tight
axis square
xlim([10 18.5])
ylim([10 18.5])

saveas(gcf, 'fig3_evolScatter.svg');


