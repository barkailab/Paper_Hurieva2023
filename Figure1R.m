%% Load samples
load('20231005_medians_BH_RM.mat','medianSumPromNewAll')
     
GP=load('group_imp.mat');
subtelomereGenes = GP.groups{1, 23}{1, 2}{1, 61};
fldNames = fieldnames(medianSumPromNewAll); %get field names from the struct
for i = 1:length(fldNames)
    medianSumProm.(fldNames{i}) = medianSumPromNewAll.(fldNames{i});
end
strainNamesNew = fldNames;
clearvars fldNames i


%% Strains order for plots
currStrainsTrunc = {'GLN3_250N_WT'; 'GLN3_200N_WT'; 'GLN3_100N_WT'; 'GLN3_050N_WT';
                                     'GLN3_orig';
                                     'GLN3_050C_WT'; 'GLN3_100C_WT'; 'GLN3_150C_WT'; 'GLN3_200C_WT'; 'GLN3_250C_WT'; 'GLN3_300C_WT'; 'GLN3_350C_WT';
                                      'GLN3_050N350C_WT'; 'GLN3_100N350C_WT'; 'GLN3_200N350C_WT'; 'GLN3_250N350C_WT';'GLN3_nonDBD_WT'};
                                  

%% Selected sumProms and normProfiles
for i = 1:length(currStrainsTrunc)
       medianSumPromShort(:,i) = medianSumProm.(currStrainsTrunc{i}); 
end
clearvars i



%% Targets and promoter information

load('gln3sites.mat')
load('promoterLengthsORF')
load('promoterIDXvecF')





%% Figure 1C - Promoter preferences

load('nBS.mat') % col1: #hap2 BS in promoter, col2: #gln3 BS in promoter
load('promType.mat') % annotation for each promoter, 1,2: good promoter, >2: bad promoter (e.g. telomeric region)

figure('Color', [1 1 1], 'Renderer', 'painters')
[~,idx]=ismember({'GLN3_orig','GLN3_250N350C_WT'},currStrainsTrunc)
scatter(log2(medianSumPromShort(promType<3, idx(1))+700),log2(medianSumPromShort(promType<3, idx(2))+700), 300, nBS(promType<3,2), 'filled','MarkerEdgeColor',[1 1 1].*.75)
caxis([-0.5 8.5])
colormap(gca, brighten(brewermap(9,'OrRd'),0.1))
axis tight
xlim(xlim+[-.2 .2])
ylim(ylim+[-.2 .2])
ylabel("Gln3 DBD-only")
xlabel("Gln3 FL")
ylabel(colorbar(),'# invitro motifs')
title(sprintf('%.2f', corr(log2(medianSumPromShort(promType<3,idx(1))+700),log2(medianSumPromShort(promType<3,idx(2))+700))))
saveas(gcf, 'fig1B_comp.svg');

[~,idx]=ismember({'GLN3_orig','GLN3_nonDBD_WT'},currStrainsTrunc)
scatter(log2(medianSumPromShort(promType<3, idx(1))+700),log2(medianSumPromShort(promType<3, idx(2))+700), 300, nBS(promType<3,2), 'filled','MarkerEdgeColor',[1 1 1].*.75)
caxis([-0.5 8.5])
colormap(gca, brighten(brewermap(9,'OrRd'),0.1))
axis tight
xlim(xlim+[-.2 .2])
ylim(ylim+[-.2 .2])
ylabel("Gln3 DBD-only")
xlabel("Gln3 FL")
ylabel(colorbar(),'# invitro motifs')
title(sprintf('%.2f', corr(log2(medianSumPromShort(promType<3,idx(1))+700),log2(medianSumPromShort(promType<3,idx(2))+700))))
saveas(gcf, 'fig1B_comp.svg');


%% Figure 1D - correlation to WT individual repeats 
clearvars -except promType GP 
repeatData=load('repeatData.mat')


currStrainsTrunc = {'GLN3_250N_WT'; 'GLN3_200N_WT'; 'GLN3_100N_WT'; 'GLN3_050N_WT';
'GLN3_orig';
'GLN3_050C_WT'; 'GLN3_100C_WT'; 'GLN3_150C_WT'; 'GLN3_200C_WT'; 'GLN3_250C_WT'; 'GLN3_300C_WT'; 'GLN3_350C_WT';
'GLN3_050N350C_WT'; 'GLN3_100N350C_WT'; 'GLN3_200N350C_WT'; 'GLN3_250N350C_WT'};

[~,repeatData.samples.strainOrder]=ismember(repeatData.samples.strain,currStrainsTrunc)
repeats=repeatData.samples(repeatData.samples.strainOrder>0,:)
sumProm=repeatData.sumPromRep(:,repeatData.samples.strainOrder>0)
[repeats,idx]=sortrows(repeats,"strainOrder")
sumProm=sumProm(:,idx);

sumPromOrig=mean(sumProm(promType<3,ismember(repeats.strain,'GLN3_orig')),2);
corrWt= corr(log2(700+sumProm(promType<3,:)),log2(700+sumPromOrig))
meanCr=accumarray(repeats.strainOrder,corrWt,[],@mean)
stdCr=accumarray(repeats.strainOrder,corrWt,[],@std)
figure('Color', [1 1 1], 'Renderer', 'painters')
errorbar(1:length(currStrainsTrunc),meanCr,stdCr,'linewidth',2,'CapSize',0)
hold on
scatter(1:length(currStrainsTrunc),meanCr,'.')
ylabel("Correlation to Gln3 FL")
xlabel("Gln3 truncation")

currStrainsTruncUre2=regexprep(currStrainsTrunc,{'_WT','_orig'},{'_ure2','_norm_ure2'})

[~,repeatData.samples.strainOrder]=ismember(repeatData.samples.strain,currStrainsTruncUre2)
repeats=repeatData.samples(repeatData.samples.strainOrder>0,:)
sumProm=repeatData.sumPromRep(:,repeatData.samples.strainOrder>0)
[repeats,idx]=sortrows(repeats,"strainOrder")
sumProm=sumProm(:,idx);

%sumPromOrig=mean(sumProm(promType<3,ismember(repeats.strain,'GLN3_norm_ure2')),2);
corrWt= corr(log2(700+sumProm(promType<3,:)),log2(700+sumPromOrig))
hold on
meanCr=accumarray(repeats.strainOrder,corrWt,[],@mean)
stdCr=accumarray(repeats.strainOrder,corrWt,[],@std)
errorbar(1:length(currStrainsTruncUre2),meanCr,stdCr,'linewidth',2,'CapSize',0)
scatter(1:length(currStrainsTruncUre2),meanCr,'.')
ylabel("Correlation to Gln3 FL")
xlabel("Gln3 truncation")
xlim([0 16]+.5)
saveas(gcf,'crwithError.svg')

%% Figure S1A
clear all
load('20231005_medians_BH_RM.mat','medianSumPromNewAll')
load('promType.mat')
load('nBS.mat') % col1: #hap2 BS in promoter, col2: #gln3 BS in promoter
subplot(1,1,1)
scatter(log2(medianSumPromNewAll.GLN3_orig(promType<3)+700),log2(medianSumPromNewAll.GLN3_norm_ure2(promType<3)+700), 300, nBS(promType<3,2), 'filled','MarkerEdgeColor',[1 1 1].*.75)
caxis([-0.5 8.5])
axis tight
xlim(xlim+[-.2 .2])
ylim(ylim+[-.2 .2])



%% Figure S1B
clear all
load('repeatData.mat')
load('promType.mat') % annotation for each promoter, 1,2: good promoter, >2: bad promoter (e.g. telomeric region)
load('targetFresh.mat')

finalStrains={'GLN3_orig','GLN3_norm_ure2','Gln3_Full_TEF_Mnase_HO_loc_nodel_native_loc_TF','GLN3_250N350C_WT','GLN3_250N350C_ure2','Gln3DBD_Gcn4AD_TEF_Mnase_HO_loc_del_native_loc_TF'}
[~,samples.uID]=ismember(regexp(samples.name,'.*(?=_.\.mat)','match','once'),finalStrains)
files=samples(samples.uID>0,:)
sumProm=sumPromRep(:,samples.uID>0);


mFac=12./accumarray(files.uID,1)
sampleOrder=cell2mat(arrayfun(@(x)sort(repmat(find(files.uID==x),mFac(x),1))',1:max(files.uID),'UniformOutput',0))
bL=find(diff(sampleOrder)~=0)
figure('Color', [1 1 1], 'Renderer', 'painters')
imagesc(corr(log2(700+sumProm(targets.geneId,sampleOrder)),'rows','pairwise'))
colorbar()
hold on
plot(xlim',bL+[.5;.5],'r-','linewidth',1)
plot(bL+[.5;.5],ylim','r-','linewidth',1)
title('targets-log')
%colormap(gca,brighten(brewermap(128,'blues'),0.3))
caxis([0 1])


%% Figure 1B - requires profile data and Figure S1C
% Individual promoters 

all_norm = '/home/labs/barkailab/bohdana/001_sequencing_data/004_medians/normProfiles/'; 
for i = 1:numel(currStrainsTrunc)
    load([all_norm currStrainsTrunc{i} '.mat'])
    normProfileAll(:,i) = normProfile;
    clearvars normProfile nucReads oldName sumProm
end
clearvars i

DBD = load('/home/labs/barkailab/bohdana/001_sequencing_data/004_medians/normProfiles/GLN3_250N350C_WT.mat', 'normProfile');
FL = load('/home/labs/barkailab/bohdana/001_sequencing_data/004_medians/normProfiles/GLN3_orig.mat', 'normProfile');
nonDBD = load('/home/labs/barkailab/bohdana/001_sequencing_data/004_medians/normProfiles/GLN3_nonDBD_WT.mat', 'normProfile');

metaPro=metaProfilePromLenDivya(normProfileAll,'promEnd', 'position', 'afterTss', 100, 'promLen', promoterLengthsORF);
metaBS=metaProfilePromLenDivya(full(sparse(round(gln3sites.pos(gln3sites.score>9)),1,1,height(normProfileAll),1)),'promEnd','position','afterTss',100,'promLen',promoterLengthsORF);

figure('Color',[1 1 1], 'Renderer', 'painters')

plotCols=[0.2532    0.6006    0.7758;0.9929    0.7564    0.5529];
intGene = {'MEP2'};
[~,geneId]=ismember(intGene,GP.gene_infoR64.nameNew);


g=0;
for strain={ 'GLN3_orig',  'GLN3_250N350C_WT','GLN3_nonDBD_WT'}
    g=g+1;
    subplot(3,1,g)
    hold off
    plot([-max(promoterLengthsORF):100], metaPro(geneId,:,ismember(currStrainsTrunc,strain)))
    hold on
    scatter(find(metaBS(geneId,:)>0)-max(promoterLengthsORF), zeros(sum(metaBS(geneId,:)>0),1), [], plotCols(2,:),'filled')
    hold on
    plot((GP.gene_infoR64.felixTss(geneId, 2)-GP.gene_infoR64.position(geneId, 2)).*[1;1].*GP.gene_infoR64.dir(geneId),ylim,'k--')
    axis tight
    ylim([0 1500])
    %saveas(gcf, 'fig1_normProfDBD.svg');
    sgtitle(intGene)
end
sum(promType<3)
clearvars metaPro metaBS

%  Mean motifs requires profile data
bsSur=300;
idx = find(gln3sites.prom==1 & gln3sites.score>9);

figure('Color', [1 1 1], 'Renderer', 'painters')
imMat = reshape(DBD.normProfile(gln3sites.pos(idx)+[-bsSur:bsSur], :),[],bsSur*2+1);
plot(-bsSur:bsSur,movmean(mean(imMat,'omitnan'),21))
ylim([0 10])
saveas(gcf, 'fig1_motifDBD.svg');

clearvars imMat



%% Figure 1F - requires profile data and Figure S1D
% motif binding individual repeats
clear all
samples=struct2table([dir('/home/labs/barkailab/bohdana/001_sequencing_data/003_checProfiles/GLN3_*.mat')])
samples.strain=regexp(samples.name,'.*(?=_[a-zA-Z0-9]\.mat$)','match','once')
goodRpts=readtable('/home/labs/barkailab/bohdana/001_sequencing_data/001_Tables/CombChecAll.xlsx')
currStrainsTrunc = {'GLN3_250N_WT'; 'GLN3_200N_WT'; 'GLN3_100N_WT'; 'GLN3_050N_WT';
    'GLN3_orig';
    'GLN3_050C_WT'; 'GLN3_100C_WT'; 'GLN3_150C_WT'; 'GLN3_200C_WT'; 'GLN3_250C_WT'; 'GLN3_300C_WT'; 'GLN3_350C_WT';
    'GLN3_050N350C_WT'; 'GLN3_100N350C_WT'; 'GLN3_200N350C_WT'; 'GLN3_250N350C_WT'};
keepStrains=ismember(samples.strain,currStrainsTrunc)&ismember(samples.name,goodRpts.sampleNames);
samples=samples(keepStrains,:)
[~,samples.strainOrder]=ismember(samples.strain,currStrainsTrunc);
samples=sortrows(samples,'strainOrder','ascend')
for i=1:height(samples)
    temp=load([samples.folder{i},'/',samples.name{i}]);
    normProfileAll(:,i)=temp.normProfile;
end
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/gln3sites.mat')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/promoterLengthsORF.mat','promoterLengthsORF')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/promoterIDXvecFimp.mat','promoterIDXvecF')
clearvars -except currStrainsTrunc gln3sites normProfileAll promoterIDXvecF promoterLengthsORF samples
gln3sites.type=promoterIDXvecF(round(gln3sites.pos));

bsSize=25;
selBs = gln3sites.score>9 & gln3sites.type==1;
bsProfiles = reshape(normProfileAll(gln3sites.pos(selBs)+[-bsSize:bsSize],:), [], 2*bsSize+1, width(normProfileAll));
bsSig=squeeze(mean(bsProfiles,2));

imMat = bsSig ./ mean(normProfileAll(promoterIDXvecF==1, :));
nAll=12;
imMatFull=nan(height(imMat),nAll.*numel(currStrainsTrunc))
for i=1:numel(currStrainsTrunc)
    idx=sort(repmat(find(samples.strainOrder==i),nAll./sum(samples.strainOrder==i),1),'ascend')';
    imMatFull(:,i*12-11:i*12)=imMat(:,idx);
end
hold off
wtStrain=find(ismember(currStrainsTrunc,'GLN3_orig'))
[~,idx]=sort(mean(imMat(:,samples.strainOrder==5),2),'descend')
imagesc(log2(imMatFull(idx,:)))
caxis([1 7])
colormap(gca,[1 1 1;brighten(brewermap(128,'OrRd'),-0.3)])
hold on
plot(reshape([.5:width(imMatFull)-.5; 1.5:width(imMatFull)+.5], [], 1), reshape(repmat(sum(imMatFull>2),2,1), [], 1), '-', 'Color', [1 1 1].*0.5)
hold on
plot([0:nAll:width(imMatFull)]+[.5;.5],ylim',':','Color',[1 1 1].*.75)
ylabel('promoter motif occurences')
xlabel('Gln3 truncations')
ylabel(colorbar(), 'log2 read enr.')
yticks([])
saveas(gcf,'motifInd.svg')
nBS=sum(imMat>2);



meanImMat=cell2mat(arrayfun(@(x)mean(imMat(:,samples.strainOrder==x),2),1:max(samples.strainOrder),'UniformOutput',false))
imNan=meanImMat;
imNan(imNan<2)=nan;
yVec=arrayfun(@(x)mean(sum(imMat(:,samples.strainOrder==x)>=2),2),1:max(samples.strainOrder));
eVec=arrayfun(@(x)std(sum(imMat(:,samples.strainOrder==x)>=2),[],2),1:max(samples.strainOrder));


subplot(1,1,1)
subplot(1,3,1)
scatter(1:16,yVec,200,diff(quantile(log2(imNan),[0.5,0.99])),'filled')
hold on
errorbar(1:16, yVec,eVec,'LineStyle','none','CapSize',0)
ylabel('# bound motifs')
title('top1% vs. median')
colorbar()

subplot(1,3,2)
scatter(1:16,yVec,200,diff(quantile(log2(imNan),[0.5,0.99])),'filled')
hold on
errorbar(1:16, yVec,eVec,'LineStyle','none','CapSize',0)

ylabel('# bound motifs')
title('std log space')
colorbar()

subplot(1,3,3)
scatter(1:16,yVec,200,diff(quantile(log2(imNan),[0.5,0.99])),'filled')
hold on
errorbar(1:16, yVec,eVec,'LineStyle','none','CapSize',0)
ylabel('# bound motifs')
title('std lin space')
colorbar()


subplot(1,1,1)
for i=1:width(imMat)
    subplot(2,1,1)
    Violin(max(log2(imMat(:,i)),0),i,'ShowData',false)
    
    subplot(2,1,2)
    Violin(max(log2(imMat(imMat(:,i)>2,i)),0),i,'ShowData',false)
end
subplot(1,1,1)
for i=1:width(imMat)
subplot(4,4,i)
histogram(max(log2(imMat(:,i)),-1),[-1.1:.2:10.1])
end


%%
GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');
load('seqOrd.mat')
