%% Load samples
load(['/home/labs/barkailab/bohdana/001_sequencing_data/004_medians/20230704_medians_BH_RM'], 'medianSumPromNewAll')
     
GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');
subtelomereGenes = GP.groups{1, 23}{1, 2}{1, 61};
fldNames = fieldnames(medianSumPromNewAll); %get field names from the struct
for i = 1:length(fldNames)
    medianSumProm.(fldNames{i}) = medianSumPromNewAll.(fldNames{i});
    medianSumProm.(fldNames{i}) (subtelomereGenes) = nan;
end
strainNamesNew = fldNames;
clearvars fldNames i


%% Strains order for plots                                  
currStrainsPhospA = {'GLN3_250NA_wt'; 'GLN3_200NA_wt'; 'GLN3_150NA_wt'; 'GLN3_100NA_wt'; 'GLN3_50NA_wt';
                                    'GLN3_orig';
                                     'GLN3_50CA_wt'; 'GLN3_100CA_wt'; 'GLN3_150CA_wt'; 'GLN3_200CA_wt';
                                     'GLN3_300CA_wt'; 'GLN3_350CA_wt'; 'GLN3_375CA_wt';};
                                 
currStrainsPhospD = {'GLN3_250ND_wt'; 'GLN3_200ND_wt'; 'GLN3_150ND_wt'; 'GLN3_100ND_wt'; 'GLN3_50ND_wt'; 
                                    'GLN3_orig';
                                     'GLN3_50CD_wt'; 'GLN3_100CD_wt'; 'GLN3_150CD_wt'; 'GLN3_200CD_wt';
                                     'GLN3_300CD_wt'; 'GLN3_350CD_wt'; 'GLN3_375CD_wt'; };
                                  
currStrainsEvol = {'GLN3_orig'; 'GLN3_Smik_wt'; 'GLN3_Skud_wt'; 'GLN3_gla_wt'; 'GLN3_klu_wt'; 'GLN3_lac_wt';};


%% Selected sumProms
for i = 1:length(currStrainsPhospA)
       medianSumPromA(:,i) = medianSumProm.(currStrainsPhospA{i}); 
end
clearvars i

for i = 1:length(currStrainsPhospD)
       medianSumPromD(:,i) = medianSumProm.(currStrainsPhospD{i}); 
end
clearvars i

for i = 1:length(currStrainsEvol)
       medianSumPromEvol(:,i) = medianSumProm.(currStrainsEvol{i}); 
end
clearvars i


%% Colors and targets 
addpath('/home/labs/barkailab/bohdana/001_sequencing_data/000_Scripts/functions')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/targetFresh.mat')

color_scheme = brewermap(11, 'Spectral');
cluster_colors=color_scheme([8:11],:);


%% Conditions (Figure 3B)
selSmp=find((ismember(sTable.class,[1,2,4,5,6])&sTable.exp==3)|contains(sTable.cond,'od3'))
figure('Color',[1 1 1],'Renderer','painters','Position',[1905 554 1533 301])

bL=[find(diff(sTable.exp(selSmp))),0];
bLmin=find(diff(sTable.class(selSmp))~=0)
bL=bLmin;
bLmin=[1:numel(selSmp)]'

subplot(1,2,1)
imagesc(corr(log2(sumProm(targets.geneId(targets.bad==0),selSmp)+700)))
set(gca,'Ytick',1:numel(selSmp),'YTickLabel',sTable.label(selSmp))
hold on
plot(bL'+0.5.*[1;1],ylim','-','LineWidth',2,'Color',[1 1 1].*0.85)
plot(xlim',bL'+0.5.*[1;1],'-','LineWidth',2,'Color',[1 1 1].*0.85)
plot(bLmin'+0.5.*[1;1],ylim','-','LineWidth',1,'Color',[1 1 1].*0.85)
plot(xlim',bLmin'+0.5.*[1;1],'-','LineWidth',1,'Color',[1 1 1].*0.85)
caxis([0.8 1])
title('good target genes -log2')
ylabel(colorbar(),'correlation across targets')
colormap(gca,brighten((brewermap(128,'blues')),0.1))

for e=3
    subplot(1,2,2)
    hold off
    selSmp=sTable.exp==e & ismember(sTable.class,[1,2,4,5,6]);
    targets.dyn=-std(log2(sumProm(targets.geneId,selSmp)+700),[],2);
    targets2=sortrows(targets,{'cluster','dyn'});
    imagesc(log2(sumProm(targets2.geneId,selSmp)+700)-log2(mean(sumProm(targets2.geneId,selSmp),2)+700),[-1 1].*1.15)
    set(gca,'Xtick',1:sum(selSmp),'XTickLabel',sTable.label(selSmp))
    hold on
    blTemp=find(diff(sTable.class(selSmp)))
    plot(xlim',gbL'+[.5;.5],'k-','LineWidth',2)
    plot(blTemp'+0.5.*[1;1],ylim','k:','LineWidth',2)
    ylabel(colorbar(),'log2 change')
    colormap(gca,brighten(flipud(brewermap(128,'RdBu')),0.3))
end


%% Phosphomimics relative binding (Figure 3D)
medianSumPromPh = medianSumPromA;

relativeBinding = log2(medianSumPromPh(targets.geneId,:)+700)-log2(medianSumPromPh(targets.geneId,6)+700);
relativeBinding(any(isnan(relativeBinding),2),:) = [];

figure('Color', [1 1 1], 'Renderer', 'painters')
bL=find(diff(targets.newCluster));
imagesc(relativeBinding)
colorbar
colormap(gca,brighten(flipud(brewermap(1000,'RdBu')),0.3))
caxis([-4.5 4.5])
xticks(1:length(medianSumPromPh))
hold on
plot(xlim',bL'+[.5;.5],'k-','linewidth',1)
plot([1.5:numel(medianSumPromPh)].*[1;1],ylim','-','linewidth',.5,'Color',[1 1 1].*0.75)
yticks([])
ylabel("Gln3 targets by truncation clusters (n=184)")
title('relative binding to A', 'FontWeight', 'normal',  'FontSize', 11)
saveas(gcf, 'Fig3d.svg')


%% Evolution heatmap (Figure 3E)
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


%% Evolution scatter plot (Figure 3F)
figure('Color',[1 1 1],'Renderer','painters')
scatter(log2(medianSumPromEvol(targets.geneId,1)+700), log2(medianSumPromEvol(targets.geneId,4)+700), 100, targets.cluster,'filled','MarkerEdgeColor',[1 1 1].*.35)
%scatter(log2(medianSumPromEvol(targets.geneId,1)+700), log2(medianSumPromEvol(targets.geneId,6)+700), 100, targets.cluster,'filled','MarkerEdgeColor',[1 1 1].*.35)
colormap(gca, cluster_colors)
xlabel("S. cerevisiae")
ylabel("C. glabrata")
%ylabel("K. lactis")
plot([10 18.5], [10 18.5], '--', 'color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off')
axis tight
axis square
xlim([10 18.5])
ylim([10 18.5])
saveas(gcf, 'fig3_evolScatter.svg');


