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
currStrainsTrunc = {'GLN3_250N_WT'; 'GLN3_200N_WT'; 'GLN3_100N_WT'; 'GLN3_050N_WT';
                                     'GLN3_orig';
                                     'GLN3_050C_WT'; 'GLN3_100C_WT'; 'GLN3_150C_WT'; 'GLN3_200C_WT'; 'GLN3_250C_WT'; 'GLN3_300C_WT'; 'GLN3_350C_WT';
                                      'GLN3_050N350C_WT'; 'GLN3_100N350C_WT'; 'GLN3_200N350C_WT'; 'GLN3_250N350C_WT';
                                      'GLN3_nonDBD_WT'};
                                  

%% Selected sumProms
for i = 1:length(currStrainsTrunc)
       medianSumPromShort(:,i) = medianSumProm.(currStrainsTrunc{i}); 
end
clearvars i


%% Gln3 target definition by clustering
[~,top100Wt]=maxk(medianSumPromShort(:,5),100);
zTh=3.5;
thLevel=mean(log2(medianSumPromShort+1),'omitnan')+zTh*std(log2(medianSumPromShort+1),[],'omitnan');
selTargets=union(find(sum(log2(medianSumPromShort(:,currStrainsTrunc)+1)>thLevel(currStrainsTrunc),2)>2),top100Wt);

clusterGenes=nan(6701,1);
clusterMat=log2(500+medianSumPromShort(selTargets,currStrainsTrunc))-median(log2(500+medianSumPromShort(selTargets,currStrainsTrunc)),2);
clusterGenes(selTargets,1)=kmeans(clusterMat,4,'Distance','correlation',"Replicates",100,'MaxIter',500);
%clusterGenes saved in targetFresh.mat


%% Targets and colors
addpath('/home/labs/barkailab/bohdana/001_sequencing_data/000_Scripts/functions')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/targetFresh.mat')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/promoterIDXvecFimp.mat','promoterIDXvecF')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/promoterLengthsORF.mat','promoterLengthsORF')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/gln3sites.mat')
addpath('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rossi2021')

color_scheme = brewermap(11, 'Spectral');
cluster_colors=color_scheme([8:11],:);


%% Relative binding heatmap and number of GATA motifs in promoters (Figure 2A)
metaPro=metaProfilePromLenDivya(promoterIDXvecF,'promEnd','position','afterTss',0,'promLen',promoterLengthsORF);
promType=mode(metaPro,2);

relativeBinding = log2(medianSumPromShort(targets.geneId,:)+700)-log2(medianSumPromShort(targets.geneId,5)+700);
relativeBinding(any(isnan(relativeBinding),2),:) = [];

onesBS=metaProfilePromLenDivya(full(sparse(round(gln3sites.pos(gln3sites.score>9)),1,1,GP.chrIdx(18),1)),'promEnd','position','afterTss',100,'promLen',promoterLengthsORF);
nBS = sum(onesBS, 2, 'omitnan');

figure('Color', [1 1 1], 'Renderer', 'painters')
bL=find(diff(targets.newCluster));
subplot(4, 4, [1 2 5 6 9 10 13 14])
imagesc(relativeBinding)
colorbar
colormap(gca,brighten(flipud(brewermap(1000,'RdBu')),0.3))
caxis([-4.5 4.5])
xticks(1:length(currStrainsTrunc))
hold on
plot(xlim',bL'+[.5;.5],'k-','linewidth',1)
plot([1.5:numel(currStrainsTrunc)].*[1;1],ylim','-','linewidth',.5,'Color',[1 1 1].*0.75)
yticks([])
ylabel("Gln3 targets by truncation clusters (n=184)")
title('relative binding', 'FontWeight', 'normal',  'FontSize', 11)

subplot(4, 4, [3 7 11 15])
absBinding = log2(medianSumPromShort(targets.geneId,5)+700);
absBinding(any(isnan(absBinding),2),:) = [];
bL=find(diff(targets.newCluster));
imagesc(absBinding)
colorbar
colormap(brewermap(1000, 'Greens'))
xticks(1:length(currStrainsTrunc))
hold on
plot(xlim',bL'+[.5;.5],'k-','linewidth',1)
plot([1.5:numel(currStrainsTrunc)].*[1;1],ylim','-','linewidth',.5,'Color',[1 1 1].*0.75)
yticks([])
title('absolute binding Gln3 FL', 'FontWeight', 'normal',  'FontSize', 11)

subplot(4, 4, [4 8 12 16])
barh(nBS(targets.geneId(~ismember(targets.geneId, [2522 5021]))),'FaceColor',[1 1 1].*0.55)
axis tight
set(gca,'YDir','reverse')
yticks([])
xlabel('#BS')
hold on
plot(xlim',bL'+[.5;.5],'k-','linewidth',1)
title('Motifs', 'FontWeight', 'normal',  'FontSize', 11)
saveas(gcf, 'fig2_HM.svg')


%% Median relative binding by cluster (Figure 2B)
figure('Color', [1 1 1], 'Renderer', 'painters')
yticks([])
for i = 1:4
    subplot(4,4,4*i+[-3 0])
    cluster = i;    
    plot(1:length(currStrainsTrunc), median(relativeBinding(targets.cluster==cluster, :), 'omitnan'), 'LineWidth', 2, 'Color', cluster_colors(i,:))
    ylabel('')
    hold on
    ylabel('log2 vs wt')
    axis tight
    clearvars cluster
end
clearvars i
saveas(gcf, 'fig2_medBS.svg');


%% Gene ontology (Figure 2C)
addpath('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem')
pScale=flipud(gray(180));
pScale=pScale(1:128,:);%brighten(flipud),0.5)%brewermap(128,'bupu');

figure('Color',[1 1 1],'Renderer','painters','Position',[2127 288 560 361])
goType=[];
subplot(2,1,1)
goEnr=goAnalysis(targets.geneId,'BG',1:6701,'type',goType);
goEnr(contains(goEnr{:,1},'C'),:)=[];
goEnr=goEnr(contains(goEnr{:,2},{'transmembrane transporter activity';'oxidoreductase activity';'DNA binding';'cellular respiration';'amino acid transport'}),:);
hold off
for b=1:size(goEnr,1)
    pIdx=round(rescale(-log10(goEnr.pVals(b)),0.5,128.4,'InputMin',-log10(.05),'InputMax',10));
    barh(b,goEnr.n(b),'BarWidth',0.8,'FaceColor',pScale(pIdx,:))
    hold on
end
text(goEnr.n,1:size(goEnr,1),num2str(goEnr.n),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',6)
yticks(1:size(goEnr,1))
if numel(goType)==0
    ytL=strcat(goEnr{:,1},',',regexp(goEnr{:,2},'.{1,20}','match','once'));
else
    ytL=regexp(goEnr{:,2},'.{1,20}','match','once');
end

clusterGenes=full(sparse(targets.geneId,1,targets.newCluster,6701,1));
set(gca,'yticklabel',ytL,'FontSize',6,'ydir','reverse')
xlabel('n genes')
title(sprintf('all Genes n= %d',numel(targets.geneId)))
caxis([-log10(0.05) 10])
ylabel(colorbar(),'-log10 pVal')
colormap(gca,pScale)
ylim([.5 5.5])
xlim([0 50])
goDist=cell2mat(arrayfun(@(b)accumarray(clusterGenes(goEnr.genes{b}),1,[max(clusterGenes),1]),1:size(goEnr,1),'UniformOutput',false));
subplot(2,1,2)
hold off
b=barh(goDist','stacked');
for i=1:numel(b)
    b(i).FaceColor=cluCol(i,:);
end
set(gca,'YDir','reverse')
title('Cluster distribution of genes')
ylim([.5 5.5])
xlim([0 50])


%% By gene groups (Figure 2E)
geneTable = readtable('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/FigureScripts/forFigure2.xlsx');
geneTable.Gene = strrep(geneTable.Gene,'DUR12','DUR1%2C2');

[~,geneTable.geneId] = ismember(geneTable.Gene, GP.gene_infoR64.nameNew);
targetMarks={'','*'};

relativeBinding = log2(medianSumPromShort+700)-log2(medianSumPromShort(:,5)+700);
cluIdx=zeros(6701,1);
cluIdx(targets.geneId)=targets.cluster;

figure('Color',[1 1 1],'Renderer','painters')
g=4;
selGenes=geneTable.geneId(geneTable.Group==g);
selGenes=sortrows(table(selGenes,cluIdx(selGenes)),2);
subplot(1, 5, [1:4])
imagesc(relativeBinding(selGenes.selGenes, :))
colormap(gca,brighten(flipud(brewermap(1000,'RdBu')),0.3))
colorbar()
caxis([-4.5 4.5])
hold on
plot(xlim',[1.5:numel(selGenes)].*[1;1],'k-')
plot([1.5:numel(currStrainsTrunc)].*[1;1],ylim','k-')
set(gca, 'Ytick',1:height(selGenes),'Yticklabel', GP.gene_infoR64.nameNew(selGenes.selGenes))

subplot(1, 5, 5)
imagesc(selGenes.Var2)
colormap(gca,[[1 1 1].*.45; cluster_colors])
caxis([-.5 4.5])
hold on
plot(xlim',[1.5:numel(selGenes)].*[1;1],'k-')
saveas(gcf, 'fig2_groups4.svg');
