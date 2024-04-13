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
currStrains481 = {'GLN3_orig';'GLN3_dFWY481_580';'GLN3_dS481_580';'GLN3_dM481_580';'GLN3_dPG481_580';'GLN3_dLIV481_580';'GLN3_dKR481_580';
                                    'GLN3_edPG481_580';'GLN3_edLIV481_580';'GLN3_edNQT481_580';'GLN3_edKR481_580';
                                    'GLN3_150_250C_WT'; 'GLN3_350C_WT'; 'GLN3_d101_280_d381_730'};
                                
currStrainsPoint = {'GLN3_orig'; 'GLN3_d482_484';'GLN3_d503_505';'GLN3_d537_541';'GLN3_d569_571';'GLN3_150_250C_WT'};

currStrains101_280 = {'GLN3_350C_WT'; 'GLN3_dKR101_280_d381_730'; 'GLN3_dS101_280_d381_730';  'GLN3_dPG101_280_d381_730'; 'GLN3_dM101_280_d381_730';
                                         'GLN3_dFWY101_280_d381_730'; 'GLN3_dLIV101_280_d381_730'; 'GLN3_dNQT101_280_d381_730';
                                         'GLN3_edPG101_280_d381_730'; 'GLN3_edLIV101_280_d381_730'; 'GLN3_edNQT101_280_d381_730'; 'GLN3_edFYW101_280_d381_730';
                                         'GLN3_d101_280_d381_730'};
                  
currStrainsChimera = {'GLN3_orig'; 'GLN3_gla_wt'; 'GLN3_pCglab494_601'; 'GLN3_pCglab602_690'; 'GLN3_pCglab602_823'; 'GLN3_pCglab494_823';
                                    'GLN3_150_250C_WT'};
                                
currStrainsCompare = {'GLN3_orig'; 'GLN3_350C_WT'; 
                                         'GLN3_edLIV481_580'; 'GLN3_dLIV481_580'; 'GLN3_edPG481_580'; 'GLN3_dPG481_580';
                                          'GLN3_edKR481_580'; 'GLN3_dKR481_580';
                                           'GLN3_edFYW101_280_d381_730'; 'GLN3_dFWY101_280_d381_730'; 'GLN3_edLIV101_280_d381_730'; 'GLN3_dLIV101_280_d381_730';
                                            'GLN3_edNQT101_280_d381_730'; 'GLN3_dNQT101_280_d381_730'; 'GLN3_edPG101_280_d381_730'; 'GLN3_dPG101_280_d381_730';
                                            'GLN3_dFWY481_580';'GLN3_dS481_580';'GLN3_dM481_580';
                                             'GLN3_dKR101_280_d381_730'; 'GLN3_dS101_280_d381_730';  'GLN3_dM101_280_d381_730';};
                                
                                        
%% Combine all medians in one table
currStrains = currStrains101_280;

for i = 1:length(currStrains)
       medianSumPromShort(:,i) = medianSumProm.(currStrains{i}); 
end
clearvars i


%% Labels
currStrainsLabels101 = strrep(currStrains, 'GLN3_250N350C_WT', 'DBD');
currStrainsLabels101 = strrep(currStrainsLabels101, 'GLN3_', '');
currStrainsLabels101 = strrep(currStrainsLabels101, '101_280', ' ');
currStrainsLabels101 = strrep(currStrainsLabels101, 'd381_730', ' ');
currStrainsLabels101 = strrep(currStrainsLabels101, 'ed', '');
currStrainsLabels101 = strrep(currStrainsLabels101, 'd', '');
currStrainsLabels101 = strrep(currStrainsLabels101, '_', ' ');
currStrainsLabels101{14} = 'd101 280';


%% Colors and targets 
addpath('/home/labs/barkailab/bohdana/001_sequencing_data/000_Scripts/functions')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/targetFresh.mat')

color_scheme = brewermap(11, 'Spectral');
cluster_colors=color_scheme([8:11],:);
compare_colors=color_scheme([8 11],:);

all_colorsCt = [3 1 1 1 1 1 1 2 2 2 2 3]; %c-terminal
all_colors101 = [3 1 1 1 1 1 1 1 2 2 2 2 3 3]; %n-terminal
all_colorsChim = [2 1 1 1 1 2]; %chimera

color_scheme = brewermap(20, 'Spectral');
mut_colors = color_scheme([14 17 19],:);


%% Disorder plot (Figure 4A, top)
load(['/home/labs/barkailab/divyakr/matlab/Disorder tendency/disTendency.mat'])
index_hap2 = find(ismember(GP.gene_infoR64.name, "GLN3"));
gln3 = disTendency(:, index_gln3);
gln3 = gln3(1:730);

figure('Color',[1 1 1],'Renderer','painters')
plot(1:730, gln3)
hold on
plot([0 265], [0.5 0.5], '--', 'color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off')
xlim([0 265])
saveas(gcf, 'gln3_dis.svg')


%% AA composition (Figure 4A, middle)
aaType(aa2int('KRDEFWYLIVNQTSHPGMAC'))=[1,1,2,2,3,3,3,4,4,4,5,5,5,6,7,8,8,9,10,10];
aaType = aaType';
aaColor=[1,2,3,4,5,5,5,6,6,6];
typeColor=[
    128, 170, 255
    255, 128, 128
    198, 236, 217
   237, 214, 34
   255, 204, 255
    180 180 180]/255;
xLabel = {'KR' 'DE' 'FWY' 'LIV' 'NQT' 'S' 'H' 'PG' 'M' 'AC'}';

Nterm = 'NHMDDDIAMFDSLATTQPIDIAASNQQNGEIAQLWDFNVDQFNMTPSNSSGSATISAPNSFTSDIPQYNHGSLGNSVSKSSLFPYNSSTSNSNINQPSINNNSNTNAQSHHSFNIYKLQNNNSSSSAMNITNNNNSNNSNIQHPFLKKSDSIGLSSSNTTNSVRKNSLIKPMSSTSLANF';
Cterm = 'RRSSTSSNTSSSSKSSSRSVVPILPKPSPNSANSQQFNMNMNLMNTTNNVSAGNSVASSPRIISSANFNSNSPLQQNLLSNSFQRQGMNIPRRKMSRNAS';
protSeq = Cterm;
aaDist=accumarray(aa2int(protSeq)',1)/numel(protSeq)*100;
aaDist = aaDist';
groupDist=accumarray(aaType,aaDist);

figure('Color',[1 1 1],'Renderer','painters')
for i = 1:10
    bar(i, groupDist(i), 'FaceColor',typeColor(aaColor(i), :))
    hold on
end
set(gca,'Xtick',1:10,'Xticklabel',xLabel,'FontSize',9)
xlim([0 10.5])
ylabel('% from total sequence')
title('C terminal 481-580')
saveas(gcf, 'aa_comp_Cterm.fig')

clearvars protSeq aaDist groupDist i


%% Scatter plot of promoters (Figure 4A, bottom)
figure('Color',[1 1 1])
for i=1:4
    scatter(log2(medianSumPromShort(targets.geneId(ismember(targets.cluster,i)),1)+700),log2(medianSumPromShort(targets.geneId(ismember(targets.cluster,i)),12)+700),100,cluster_colors(targets.cluster(ismember(targets.cluster,i)),:),'filled','MarkerEdgeColor',[1 1 1].*.35)
    hold on
end
xlabel("Gln3 FL")
ylabel("GLN3 dR2")
plot([10 19], [10 19], '--', 'color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off')
axis tight
set(gcf,'renderer','painters')
saveas(gcf, 'fig4_scatter.svg');

    
%% Scatter plot of correlation to WT/DBD (Figure 4B, D, E, G)
medianSumPromTarg = medianSumPromShort(targets.geneId, :);

corr_toWT= corr(medianSumPromTarg,medianSumPromTarg(:,1), 'rows', 'pairwise');
corr_toNterm = corr(medianSumPromTarg,medianSumPromTarg(:,14), 'rows', 'pairwise');

figure('Color',[1 1 1])
all_colors = all_colors101;
scatter(corr_toWT(all_colors==1), corr_toNterm(all_colors==1), 180, mut_colors(all_colors(all_colors==1), :), 'filled', 'MarkerEdgeColor',[1 1 1].*.35,'DisplayName','deletion')
hold on
scatter(corr_toWT(all_colors==2), corr_toNterm(all_colors==2), 180, mut_colors(all_colors(all_colors==2), :), 'filled','MarkerEdgeColor',[1 1 1].*.35, 'DisplayName','spread')
hold on
scatter(corr_toWT(all_colors==3), corr_toNterm(all_colors==3), 180, mut_colors(all_colors(all_colors==3), :), 'filled','MarkerEdgeColor',[1 1 1].*.35, 'DisplayName','extremes')
text(corr_toWT, corr_toNterm, currStrainsLabels101)
xlabel("correlation to Gln3 dCt")
ylabel("correlation to Gln3 dR3 dCt")
legend('Location', 'northeast')
saveas(gcf, '4b.svg');


%% Comparison plot: d/ed (Figure 4C)
medianSumPromTarg = medianSumPromShort(targets.geneId, :);

currStrainsLabelsCompare = {'LIV'; 'PG'; 'KR'; 'FWY'; 'LIV'; 'NQT'; 'PG'; 'FWY'; 'S'; 'M'; 'KR'; 'S'; 'M'};
all_colorsComp = [2 2 2 3 3 3 3 2 2 2 3 3 3]; %compare
all_colors=all_colorsComp;

corr_delWT481 = corr(medianSumPromTarg(:, [4 6 8]), medianSumPromTarg(:,1), 'rows', 'pairwise');
corr_delWT101 = corr(medianSumPromTarg(:, [10 12 14 16]), medianSumPromTarg(:,2), 'rows', 'pairwise');
corr_delWT = [corr_delWT481; corr_delWT101];

corr_edDel = diag(corr(medianSumPromTarg(:, [3 5 7 9 11 13 15]), medianSumPromTarg(:, [4 6 8 10 12 14 16]), 'rows', 'pairwise'));

corr_edWT481 = corr(medianSumPromTarg(:, [3 5 7]), medianSumPromTarg(:,1), 'rows', 'pairwise');
corr_edWT101 = corr(medianSumPromTarg(:, [9 11 13 15]), medianSumPromTarg(:,2), 'rows', 'pairwise');
corr_edWT = [corr_edWT481; corr_edWT101];


ratio = (1-corr_edWT)./(1- corr_delWT);
ratioAll = [ratio; -0.1; -0.1; -0.1; -0.1; -0.1; -0.1];

corr_delWT481_2 = corr(medianSumPromTarg(:, [17 18 19]), medianSumPromTarg(:,1), 'rows', 'pairwise');
corr_delWT101_2 = corr(medianSumPromTarg(:, [20 21 22]), medianSumPromTarg(:,2), 'rows', 'pairwise');

corr_delAll = [corr_delWT; corr_delWT481_2; corr_delWT101_2];

figure('Color',[1 1 1])
scatter(1- corr_delAll(all_colors==2), ratioAll(all_colors==2), 180, cluster_colors(all_colors(all_colors==2), :), 'filled', 'MarkerEdgeColor',[1 1 1].*.35, 'DisplayName','R2: 481-580')
hold on
scatter(1- corr_delAll(all_colors==3), ratioAll(all_colors==3), 180, cluster_colors(all_colors(all_colors==3), :), 'filled', 'MarkerEdgeColor',[1 1 1].*.35, 'DisplayName','R3: 101-280')
hold on
text(1-corr_delAll, ratioAll, currStrainsLabelsCompare)
plot([0 1], [0 0], '--', 'color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off')
hold on
plot([0 1], [1 1], '--', 'color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off')
xlabel("1 - deletion corr to WT")
ylabel("relative effect of the spread")
legend('Location', 'northeast')
saveas(gcf, '4c.svg');


%% ST/DE in phospho region (Figure 4G, top)
Scer = 'YSSSFMAASLQQLHEQQQVDVNSNTNTNSNRQNWNSSNSVSTNSRSSNFVSQKPNFDIFNTPVDSPSVSRPSSRKSHTSLLS';
Cglab = 'YSSSFINTNPQQTQDIGNFNDDQSIQGNGTGVNSQNIRRINSNYDSPQPNFDLFRLDSNKDSPENIPDVLRSDSRLSQKSQVSHTSLLS';
ST = 'ST';
DE = 'DE';

imagesc(ismember(aa2int(Cglab),aa2int(ST)))
saveas(gcf, 'composST_glab.svg');
