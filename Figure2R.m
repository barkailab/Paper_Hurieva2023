%% Figure 2A,B,E
%% Colors and targets 

%% load repeats
clear all
cluster_colors=[0.913725490196078,0.596078431372549,0.478431372549020;0.400000000000000,0.760784313725490,0.647058823529412;0.196078431372549,0.533333333333333,0.741176470588235;0.721568627450980,0.670588235294118,0.827450980392157];
load('targetFresh.mat')
repeatData=load('repeatData.mat')

currStrainsTrunc = {'GLN3_250N_WT'; 'GLN3_200N_WT'; 'GLN3_100N_WT'; 'GLN3_050N_WT';
    'GLN3_orig';
    'GLN3_050C_WT'; 'GLN3_100C_WT'; 'GLN3_150C_WT'; 'GLN3_200C_WT'; 'GLN3_250C_WT'; 'GLN3_300C_WT'; 'GLN3_350C_WT';
    'GLN3_050N350C_WT'; 'GLN3_100N350C_WT'; 'GLN3_200N350C_WT'; 'GLN3_250N350C_WT';
    'GLN3_nonDBD_WT'};


[~,repeatData.samples.strainOrder]=ismember(repeatData.samples.strain,currStrainsTrunc)
samples=repeatData.samples(repeatData.samples.strainOrder>0,:)
sumProm=repeatData.sumPromRep(:,repeatData.samples.strainOrder>0)
[samples,idx]=sortrows(samples,"strainOrder")
sumProm=sumProm(:,idx);


clearvars -except cluster_colors color_scheme currStrainsTrunc promType samples sumProm targets

relativeBinding=log2(sumProm+700)-mean(log2(sumProm(:,ismember(samples.strain,'GLN3_orig'))+700),2)
%2 A
nRows=12;
mFac=nRows./accumarray(samples.strainOrder,1)
sampleOrder=cell2mat(arrayfun(@(x)sort(repmat(find(samples.strainOrder==x),mFac(x),1))',1:max(samples.strainOrder),'UniformOutput',0))
subplot(1,2,1)
hold off
imagesc(relativeBinding(targets.geneId,sampleOrder))
colorbar
colormap(gca,brighten(flipud(brewermap(1000,'RdBu')),0.3))
%colormap(brewermap(1000, 'RdBu'))
caxis([-4.5 4.5])
bL=find(diff(targets.newCluster));
hold on
plot(xlim',bL'+[.5;.5],'k-','linewidth',1)
plot(find(diff(sampleOrder)~=0)+[.5;.5],ylim','-','linewidth',.5,'Color',[1 1 1].*0.75)
plot([1:numel(currStrainsTrunc)].*nRows+[.5;.5],ylim','-','linewidth',1,'Color',[1 1 1].*0.5)

xticks(1:length(currStrainsTrunc))
% 2B
clusterDynamics=cell2mat(arrayfun(@(x)median(relativeBinding(targets.geneId(targets.newCluster==x),:)),[1:max(targets.newCluster)]','UniformOutput',false))
meanClusterDynamics=cell2mat(arrayfun(@(x)mean(clusterDynamics(:,samples.strainOrder==x),2),1:max(samples.strainOrder),'UniformOutput',false))
stdClusterDynamics=cell2mat(arrayfun(@(x)std(clusterDynamics(:,samples.strainOrder==x),[],2),1:max(samples.strainOrder),'UniformOutput',false))

figure('Color', [1 1 1], 'Renderer', 'painters','Position',[2006 101 1101 791])

cluster_colors=[0.913725490196078,0.596078431372549,0.478431372549020;0.400000000000000,0.760784313725490,0.647058823529412;0.196078431372549,0.533333333333333,0.741176470588235;0.721568627450980,0.670588235294118,0.827450980392157];
subplot(1,1,1)
for i = 1:4
    subplot(4,2,2*i)
    hold off
    cluster = i;
    plot(1:length(currStrainsTrunc)-1, meanClusterDynamics(i,1:end-1), 'LineWidth', 2, 'Color', cluster_colors(i,:))
    ylabel('')
    hold on
    scatter(1:length(currStrainsTrunc)-1, meanClusterDynamics(i,1:end-1), [], cluster_colors(i,:),'filled')
    errorbar(1:length(currStrainsTrunc)-1, meanClusterDynamics(i,1:end-1),stdClusterDynamics(i,1:end-1),'LineStyle','none','CapSize',0,'Color',cluster_colors(i,:))
    hold off
    ylabel('log2 vs wt')
    axis tight
    xlim([1 length(currStrainsTrunc)-1])
    clearvars cluster
end
clearvars i


saveas(gcf, 'fig2_medBS.svg');


% 2E By gene groups
geneTable = readtable('forFigure2.xlsx');
geneTable.Gene = strrep(geneTable.Gene,'DUR12','DUR1%2C2');
GP=load('group_imp.mat');
load('targetFresh.mat')

[~,geneTable.geneId] = ismember(geneTable.Gene, GP.gene_infoR64.nameNew);
targetMarks={'','*'};

cluIdx=zeros(6701,1);
cluIdx(targets.geneId)=targets.cluster;
relativeBinding=log2(sumProm+700)-mean(log2(sumProm(:,ismember(samples.strain,'GLN3_orig'))+700),2)

bL=find(diff(samples.strainOrder(sampleOrder)))

figure('Color',[1 1 1],'Renderer','painters')

subplot(1,1,1)
for g=1:max(geneTable.Group)    
    selGenes=geneTable.geneId(geneTable.Group==g);
    selGenes=sortrows(table(selGenes,cluIdx(selGenes)),2);
    subplot(4, 5, [1:4]+(g-1)*5)
    imagesc(relativeBinding(selGenes.selGenes, sampleOrder))
    colormap(gca,brighten(flipud(brewermap(1000,'RdBu')),0.3))
    colorbar()
    caxis([-4.5 4.5])
    hold on
    plot(xlim',[1.5:numel(selGenes)].*[1;1],'k-')
    plot([bL'+.5].*[1;1],ylim','k-')
    set(gca, 'Ytick',1:height(selGenes),'Yticklabel', GP.gene_infoR64.nameNew(selGenes.selGenes))
    
    subplot(4, 5, 5+(g-1)*5)
    imagesc(selGenes.Var2)
    colormap(gca,[[1 1 1].*.45; cluster_colors])
    caxis([-.5 4.5])
    hold on
    plot(xlim',[1.5:numel(selGenes)].*[1;1],'k-')
end
saveas(gcf, 'fig2_groups4.svg');

%% Gene ontology

clearvars -except targets GP cluster_colors
cluster_colors=[0.913725490196078,0.596078431372549,0.478431372549020;0.400000000000000,0.760784313725490,0.647058823529412;0.196078431372549,0.533333333333333,0.741176470588235;0.721568627450980,0.670588235294118,0.827450980392157];
GP=load('group_imp.mat');

figure('Color',[1 1 1],'Renderer','painters','Position',[2127 288 560 361])
subplot(2,1,1)
goEnr=goAnalysis(targets.geneId,'BG',find(ismember(GP.gene_infoR64.status,'Verified ORF')),'type',[])
goEnr(contains(goEnr{:,1},'C'),:)=[];
goEnr=goEnr(contains(goEnr{:,2},{'transmembrane transporter activity';'oxidoreductase activity';'DNA binding';'cellular respiration';'amino acid transport'}),:);
hold off
clusterGenes=full(sparse(targets.geneId,1,targets.newCluster,6701,1));
goDist=cell2mat(arrayfun(@(b)accumarray(clusterGenes(goEnr.genes{b}),1,[max(clusterGenes),1]),1:size(goEnr,1),'UniformOutput',false))
subplot(2,1,2)
hold off
b=barh(goDist','stacked')
for i=1:numel(b)
    b(i).FaceColor=cluster_colors(i,:);
end
set(gca,'YDir','reverse','YTick',1:height(goEnr),'YTickLabel',goEnr.Var5)
title('Cluster distribution of genes')
ylim([.5 5.5])
xlim([0 50])

%% Figure S2 A,D
clear all
load('targetFresh.mat')
GP=load('group_imp.mat');
load('promType.mat')

goodGenes=find(promType<3)
intGOTerms={'GO:0022857','GO:0003677','GO:0016491','GO:0045333','GO:0006865','GO:0006091'};
intGONames={'transmembrane transporter activity','DNA binding'  ,'oxidoreductase activity'        ,'cellular respiration','amino acid transport'      ,'generation of precursor metabolites and energy'     }
figure('Color',[1 1 1], 'Renderer', 'painters','Position',[146 507 2506 413])
subplot(1,1,1)

goEnr=goAnalysis(targets.geneId,'BG',find(goodGenes),'type',[])
subplot(1,1,1)
for type={'C','P','F'}
    selRows=ismember(goEnr.Var4,type)
    subplot(1,2,1)
    scatter(goEnr.n(selRows),-log10(goEnr.pVals(selRows)+1e-15),200,'filled')
    hold on
    subplot(1,2,2)
    scatter(goEnr.log2Enr(selRows),-log10(goEnr.pVals(selRows)+1e-15),200,'filled')
    hold on
end
[~,selRows]=ismember(intGOTerms,goEnr.Var6)
subplot(1,2,1)
text(goEnr.n(selRows),-log10(goEnr.pVals(selRows)+1e-15),goEnr.Var5(selRows))

subplot(1,2,2)
text(goEnr.log2Enr(selRows),-log10(goEnr.pVals(selRows)+1e-15),goEnr.Var5(selRows))

g=0;
for selBG={find(goodGenes),targets.geneId}
    pValMat=nan(numel(intGOTerms),5);
    enrMat=nan(numel(intGOTerms),5);
    xMat=[1:width(enrMat)].*ones(size(enrMat));
    yMat=[1:height(enrMat)]'.*ones(size(enrMat));
    for i=1:4
        goEnr=goAnalysis(targets.geneId(targets.cluster==i),'BG',selBG{1},'type',[]);
        [~,intRows]=ismember(intGOTerms,goEnr.Var6);
        pValMat(intRows>0,i)=-log10(goEnr.pVals(intRows(intRows>0))+1e-15);
        enrMat(intRows>0,i)=goEnr.log2Enr(intRows(intRows>0));
    end
    goEnr=goAnalysis(targets.geneId,'BG',selBG{1},'type',[])
    [~,intRows]=ismember(intGOTerms,goEnr.Var6);
    pValMat(intRows>0,5)=-log10(goEnr.pVals(intRows(intRows>0))+1e-15);
    enrMat(intRows>0,5)=goEnr.log2Enr(intRows(intRows>0));
        
    g=g+1
    subplot(1,3,g+1)
    scatter(xMat(:),yMat(:),rescale(enrMat(:),1,10,'InputMax',5,'InputMin',0).*50,pValMat(:),'filled')
    xlim([0 width(enrMat)]+0.5)
    ylim([0 height(enrMat)]+0.5)
    caxis([0 10])
    if g==1
    set(gca,'YDir','reverse','XAxisLocation','top','XTick',1:width(enrMat),'YTick',1:height(enrMat),'YTickLabel',intGONames);
    else
        set(gca,'YDir','reverse','XAxisLocation','top','XTick',1:width(enrMat),'YTick',1:height(enrMat),'YTickLabel',[]);
    end
    xlabel('Cluster')
    colorbar()
    set(gca,'box','on')
    colormap(gca,brewermap(128,'BuPu'))
end
 

%% Figure S2 B,C

clear all
load('targetFresh.mat')
repeatData=load('repeatData.mat')

currStrainsTrunc = {'GLN3_250N_WT'; 'GLN3_200N_WT'; 'GLN3_100N_WT'; 'GLN3_050N_WT';
'GLN3_orig';
'GLN3_050C_WT'; 'GLN3_100C_WT'; 'GLN3_150C_WT'; 'GLN3_200C_WT'; 'GLN3_250C_WT'; 'GLN3_300C_WT'; 'GLN3_350C_WT';
'GLN3_050N350C_WT'; 'GLN3_100N350C_WT'; 'GLN3_200N350C_WT'; 'GLN3_250N350C_WT'};

[~,repeatData.samples.strainOrder]=ismember(repeatData.samples.strain,currStrainsTrunc)
samples=repeatData.samples(repeatData.samples.strainOrder>0,:)
sumProm=repeatData.sumPromRep(:,repeatData.samples.strainOrder>0)
[samples,idx]=sortrows(samples,"strainOrder")
sumProm=sumProm(:,idx);

clearvars -except cluster_colors color_scheme currStrainsTrunc promType samples sumProm targets

relativeBinding=log2(sumProm+700)-mean(log2(sumProm(:,ismember(samples.strain,'GLN3_orig'))+700),2)
ddSP=relativeBinding-mean(relativeBinding(targets.geneId,:));
[crWT,pWT]=corr(ddSP(targets.geneId,:)');
subplot(2,3,1)
imagesc(crWT,[0 1.2])
colormap(gca,brighten(brewermap(128,'YlGn'),0.2))
bL=find(diff(targets.newCluster));
hold on
plot(xlim',bL'+[.5;.5],'k-','linewidth',1)
plot(bL'+[.5;.5],xlim','k-','linewidth',1)

subplot(2,3,2)
imagesc(-log10(pWT),[0 7])
hold on
plot(xlim',bL'+[.5;.5],'r-','linewidth',1)
plot(bL'+[.5;.5],xlim','r-','linewidth',1)
colormap(gca,brewermap(128,'BuPu'))

subplot(2,3,3)
selPW=((targets.newCluster-targets.newCluster')==0);
hold off
imagesc(-log10(pWT).*selPW,[0 7])
hold on
plot(xlim',bL'+[.5;.5],'r-','linewidth',1)
plot(bL'+[.5;.5],xlim','r-','linewidth',1)
colormap(gca,brewermap(128,'BuPu'))
colorbar()
saveas(gcf,'crMapWT.svg')
%% ure2 
clearvars -except targets cluster_colors
repeatData=load('repeatData.mat')
currStrainsTrunc = {'GLN3_250N_ure2'; 'GLN3_200N_ure2'; 'GLN3_100N_ure2'; 'GLN3_050N_ure2';
    'GLN3_norm_ure2';
    'GLN3_050C_ure2'; 'GLN3_100C_ure2'; 'GLN3_150C_ure2'; 'GLN3_200C_ure2'; 'GLN3_250C_ure2'; 'GLN3_300C_ure2'; 'GLN3_350C_ure2';
    'GLN3_050N350C_ure2'; 'GLN3_100N350C_ure2'; 'GLN3_200N350C_ure2'; 'GLN3_250N350C_ure2';
    'GLN3_nonDBD_ure2'};

[~,repeatData.samples.strainOrder]=ismember(repeatData.samples.strain,currStrainsTrunc)
samples=repeatData.samples(repeatData.samples.strainOrder>0,:)
sumProm=repeatData.sumPromRep(:,repeatData.samples.strainOrder>0)
[samples,idx]=sortrows(samples,"strainOrder")
sumProm=sumProm(:,idx);


clearvars -except cluster_colors currStrainsTrunc promType samples sumProm targets repeatData

relativeBinding=log2(sumProm+700)-mean(log2(sumProm(:,ismember(samples.strain,'GLN3_norm_ure2'))+700),2)
%%
nRows=12;
mFac=nRows./accumarray(samples.strainOrder,1)
sampleOrder=cell2mat(arrayfun(@(x)sort(repmat(find(samples.strainOrder==x),mFac(x),1))',1:max(samples.strainOrder),'UniformOutput',0))
figure('Color', [1 1 1], 'Renderer', 'painters','Position',[2006 101 1101 791])
subplot(1,2,1)
hold off
imagesc(relativeBinding(targets.geneId,sampleOrder))
colorbar
colormap(gca,brighten(flipud(brewermap(1000,'RdBu')),0.3))
%colormap(brewermap(1000, 'RdBu'))
caxis([-4.5 4.5])
bL=find(diff(targets.newCluster));
hold on
plot(xlim',bL'+[.5;.5],'k-','linewidth',1)
plot(find(diff(sampleOrder)~=0)+[.5;.5],ylim','-','linewidth',.5,'Color',[1 1 1].*0.75)
plot([1:numel(currStrainsTrunc)].*nRows+[.5;.5],ylim','-','linewidth',1,'Color',[1 1 1].*0.5)
xT=movmean([0;find(diff(samples.strainOrder(sampleOrder))>0);numel(sampleOrder)],2,'Endpoints','discard')
xticks(xT)

%% correaltion within cluster

ddSP=relativeBinding-mean(relativeBinding(targets.geneId,:));
[crWT,pWT]=corr(ddSP(targets.geneId,:)');
subplot(2,2,2)
imagesc(crWT,[0 1.2])
colormap(gca,brighten(brewermap(128,'YlGn'),0.2))
bL=find(diff(targets.newCluster));
hold on
plot(xlim',bL'+[.5;.5],'k-','linewidth',1)
plot(bL'+[.5;.5],xlim','k-','linewidth',1)
selPW=((targets.newCluster-targets.newCluster')==0);
subplot(2,2,4)
imagesc(-log10(pWT).*selPW,[0 10])
saveas(gcf,'CrMapUre2.svg')
colormap(gca,brewermap(128,'BuPu'))

 

