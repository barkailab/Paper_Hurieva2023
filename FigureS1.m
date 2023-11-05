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
currStrainsTruncUre2 = {'GLN3_250N_ure2'; 'GLN3_200N_ure2'; 'GLN3_100N_ure2'; 'GLN3_050N_ure2';
                                     'GLN3_norm_ure2';
                                     'GLN3_050C_ure2'; 'GLN3_100C_ure2'; 'GLN3_150C_ure2'; 'GLN3_200C_ure2'; 'GLN3_250C_ure2'; 'GLN3_300C_ure2'; 'GLN3_350C_ure2';
                                      'GLN3_050N350C_ure2'; 'GLN3_100N350C_ure2'; 'GLN3_200N350C_ure2';};
                                  
currStrainsUre2 =  { 'GLN3_orig'; 'GLN3_norm_ure2';};
                                  

%% Selected sumProms and normProfiles
currStrainsShort = currStrainsTruncUre2;

for i = 1:length(currStrainsShort)
       medianSumPromShort(:,i) = medianSumProm.(currStrainsShort{i}); 
end
clearvars i


all_norm = '/home/labs/barkailab/bohdana/001_sequencing_data/004_medians/20230629_normProfiles/'; 
for i = 1:numel(currStrainsShort)
    load([all_norm currStrainsShort{i} '.mat'])
    normProfileAll(:,i) = normProfile;
    clearvars normProfile nucReads oldName sumProm
end
clearvars i


%% Targets and promoter information
addpath('/home/labs/barkailab/bohdana/001_sequencing_data/000_Scripts/functions')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/targetFresh.mat')

addpath('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rossi2021')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/gln3sites.mat')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/promoterLengthsORF.mat','promoterLengthsORF')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/promoterIDXvecFimp.mat','promoterIDXvecF')


%% Promoter preferences (Figure S1A)
metaPro=metaProfilePromLenDivya(promoterIDXvecF,'promEnd','position','afterTss',0,'promLen',promoterLengthsORF);
promType=mode(metaPro,2);

onesBS=metaProfilePromLenDivya(full(sparse(round(gln3sites.pos(gln3sites.score>9)),1,1,GP.chrIdx(18),1)),'promEnd','position','afterTss',100,'promLen',promoterLengthsORF);
nBS = sum(onesBS, 2, 'omitnan');

figure('Color', [1 1 1], 'Renderer', 'painters')
scatter(log2(medianSumPromShort(promType<3, 1)+700),log2(medianSumPromShort(promType<3, 2)+700), 300, nBS(promType<3), 'filled','MarkerEdgeColor',[1 1 1].*.75)
caxis([-0.5 8.5])
colormap(gca, brighten(brewermap(9,'OrRd'),0.1))
axis tight
xlim(xlim+[-.2 .2])
ylim(ylim+[-.2 .2])
ylabel("Gln3 dUre2")
xlabel("Gln3 wt")
ylabel(colorbar(),'# invitro motifs')
title(sprintf('%.2f', corr(log2(medianSumPromShort(promType<3,1)+700),log2(medianSumPromShort(promType<3,2)+700))))
saveas(gcf, 'fig1_compUre2.svg');


%% Correlation to FL (Figure S1B)
sumPromlog2 = log2(medianSumPromShort+700);
corr_toWT = corr(sumPromlog2, sumPromlog2(:,5), 'rows', 'pairwise');

figure('Color', [1 1 1], 'Renderer', 'painters')
plot(1:length(currStrainsTruncUre2), corr_toWT, 'linewidth',2)
ylabel("Correlation to Gln3 FL")
xlabel("Gln3 truncation")
%saveas(gcf, 'fig1_corrWT_ure2.svg');


%% Motif binding (Figure S1C)
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
%saveas(gcf, 'fig1_motifBS_ure2.svg');
