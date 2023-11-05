%% Load samples
load(['/home/labs/barkailab/bohdana/001_sequencing_data/004_medians/20230704_medians_BH_RM'],'medianSumPromNewAll')
     
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
                                      'GLN3_050N350C_WT'; 'GLN3_100N350C_WT'; 'GLN3_200N350C_WT'; 'GLN3_250N350C_WT';};
                                  

%% Selected sumProms and normProfiles
for i = 1:length(currStrainsTrunc)
       medianSumPromShort(:,i) = medianSumProm.(currStrainsTrunc{i}); 
end
clearvars i


all_norm = '/home/labs/barkailab/bohdana/001_sequencing_data/004_medians/20230629_normProfiles/'; 
for i = 1:numel(currStrainsTrunc)
    load([all_norm currStrainsTrunc{i} '.mat'])
    normProfileAll(:,i) = normProfile;
    clearvars normProfile nucReads oldName sumProm
end
clearvars i


%% Targets and promoter information
addpath('/home/labs/barkailab/bohdana/001_sequencing_data/000_Scripts/functions')

addpath('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rossi2021')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/gln3sites.mat')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/promoterLengthsORF.mat','promoterLengthsORF')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/promoterIDXvecFimp.mat','promoterIDXvecF')


%% Individual promoters (Figure 1A, left)
DBD = load('/home/labs/barkailab/bohdana/001_sequencing_data/004_medians/20230629_normProfiles/GLN3_250N350C_WT.mat', 'normProfile');
FL = load('/home/labs/barkailab/bohdana/001_sequencing_data/004_medians/20230629_normProfiles/GLN3_orig.mat', 'normProfile');
nonDBD = load('/home/labs/barkailab/bohdana/001_sequencing_data/004_medians/20230629_normProfiles/GLN3_nonDBD_WT.mat', 'normProfile');

metaPro=metaProfilePromLenDivya(DBD.normProfile, 'promEnd', 'position', 'afterTss', 100, 'promLen', promoterLengthsORF);
metaBS=metaProfilePromLenDivya(full(sparse(round(gln3sites.pos(gln3sites.score>9)),1,1,height(DBD.normProfile),1)),'promEnd','position','afterTss',100,'promLen',promoterLengthsORF);

intGene = {'AGP1'};
[~,geneId]=ismember(intGene,GP.gene_infoR64.nameNew);
plotCols=[0.2532    0.6006    0.7758;0.9929    0.7564    0.5529];

figure('Color',[1 1 1], 'Renderer', 'painters')
plot([-max(promoterLengthsORF):100], metaPro(geneId,:))
hold on
scatter(find(metaBS(geneId,:)>0)-max(promoterLengthsORF), zeros(sum(metaBS(geneId,:)>0),1), [], plotCols(2,:),'filled')
hold on
plot((GP.gene_infoR64.felixTss(geneId, 2)-GP.gene_infoR64.position(geneId, 2)).*[1;1].*GP.gene_infoR64.dir(geneId),ylim,'k--')
ylim([0 8000]) 
saveas(gcf, 'fig1_normProfDBD.svg');

clearvars metaPro metaBS


%% Mean motifs (Figure 1A, right)
bsSur=300;
idx = find(gln3sites.prom==1 & gln3sites.score>9);

figure('Color', [1 1 1], 'Renderer', 'painters')
imMat = reshape(DBD.normProfile(gln3sites.pos(idx)+[-bsSur:bsSur], :),[],bsSur*2+1);
plot(-bsSur:bsSur,movmean(mean(imMat,'omitnan'),21))
ylim([0 10])
saveas(gcf, 'fig1_motifDBD.svg');

clearvars imMat


%% Promoter preferences: DBD vs FL (Figure 1B)
metaPro=metaProfilePromLenDivya(promoterIDXvecF,'promEnd','position','afterTss',0,'promLen',promoterLengthsORF);
promType=mode(metaPro,2);

onesBS=metaProfilePromLenDivya(full(sparse(round(gln3sites.pos(gln3sites.score>9)),1,1,GP.chrIdx(18),1)),'promEnd','position','afterTss',100,'promLen',promoterLengthsORF);
nBS = sum(onesBS, 2, 'omitnan');

figure('Color', [1 1 1], 'Renderer', 'painters')
scatter(log2(medianSumPromShort(promType<3, 5)+700),log2(medianSumPromShort(promType<3, 16)+700), 300, nBS(promType<3), 'filled','MarkerEdgeColor',[1 1 1].*.75)
caxis([-0.5 8.5])
colormap(gca, brighten(brewermap(9,'OrRd'),0.1))
axis tight
xlim(xlim+[-.2 .2])
ylim(ylim+[-.2 .2])
ylabel("Gln3 DBD-only")
xlabel("Gln3 FL")
ylabel(colorbar(),'# invitro motifs')
title(sprintf('%.2f', corr(log2(medianSumPromShort(promType<3,1)+700),log2(medianSumPromShort(promType<3,2)+700))))
saveas(gcf, 'fig1B_comp.svg');


%% Correlation to FL (Figure 1D)
sumPromlog2 = log2(medianSumPromShort+700);
corr_toWT = corr(sumPromlog2, sumPromlog2(:,5), 'rows', 'pairwise');

figure('Color', [1 1 1], 'Renderer', 'painters')
plot(1:length(currStrainsTrunc), corr_toWT, 'linewidth',2)
ylabel("Correlation to Gln3 FL")
xlabel("Gln3 truncation")
saveas(gcf, 'fig1_corrWT.svg');


%% Motif binding (Figure 1E)
gln3sites.type=promoterIDXvecF(round(gln3sites.pos));

bsSize=25;
selBs = gln3sites.score>9 & gln3sites.type==1;
bsProfiles = reshape(normProfileAll(gln3sites.pos(selBs)+[-bsSize:bsSize],:), [], 2*bsSize+1, width(normProfileAll));
bsSig=squeeze(mean(bsProfiles,2));

imMat = bsSig ./ mean(normProfileAll(promoterIDXvecF==1, :));
[~,idx]=sort(sum(imMat,2,'omitnan'),'descend');

figure('Color', [1 1 1], 'Renderer', 'painters')
imagesc(log2(imMat(idx,:)))
caxis([1 7])
colormap(gca,[1 1 1;brighten(brewermap(128,'OrRd'),-0.3)])
hold on
plot(reshape([.5:width(imMat)-.5; 1.5:width(imMat)+.5], [], 1), reshape(repmat(sum(log2(imMat(idx,:))>1),2,1), [], 1), '-', 'Color', [1 1 1].*0.5)
hold on
plot([1:1:width(imMat)]+[.5;.5],ylim',':','Color',[1 1 1].*.75)
ylabel('promoter motif occurences')
xlabel('Gln3 truncations')
ylabel(colorbar(), 'log2 read enr.')
yticks([])
xlim([.5, width(normProfileAll)+.5])
saveas(gcf, 'fig1_motifBS.svg');

