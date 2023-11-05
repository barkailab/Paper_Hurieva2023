%% Load samples
load(['/home/labs/barkailab/bohdana/001_sequencing_data/004_medians/20231005_medians_BH_RM'], 'medianSumPromNewAll')
     
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
currStrainsPlot = {'GLN3_orig'; 'GLN3_050C_WT'};                               
currStrainsPlot2 = {'GLN3_150_250C_WT'; 'GLN3_d503_505'};         


 %% Combine all medians in one table
currStrains = currStrainsPlot;

for i = 1:length(currStrains)
       medianSumPromShort(:,i) = medianSumProm.(currStrains{i}); 
end
clearvars i


%% Scatter plot of promoters (Figure S3A, C)
figure('Color',[1 1 1])
for i=1:4
    scatter(log2(medianSumPromShort(targets.geneId(ismember(targets.cluster,i)),1)+700),log2(medianSumPromShort(targets.geneId(ismember(targets.cluster,i)),2)+700),100,cluster_colors(targets.cluster(ismember(targets.cluster,i)),:),'filled','MarkerEdgeColor',[1 1 1].*.35)
    hold on
end
xlabel("Gln3 FL")
ylabel("GLN3 1-680")
plot([10 19], [10 19], '--', 'color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off')
axis tight
set(gcf,'renderer','painters')
saveas(gcf, 'figs3A.svg');

