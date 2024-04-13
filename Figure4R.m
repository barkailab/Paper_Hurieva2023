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
c=0
for term={'Nterm','Cterm'}
    c=c+1;
    subplot(1,2,c)
    protSeq = eval(term{1});
    aaDist=accumarray(aa2int(protSeq)',1)/numel(protSeq)*100;
    aaDist = aaDist';
    groupDist=accumarray(aaType,aaDist);
    for i = 1:10
        bar(i, groupDist(i), 'FaceColor',typeColor(aaColor(i), :))
        hold on
    end
    set(gca,'Xtick',1:10,'Xticklabel',xLabel,'FontSize',9)
    xlim([0 10.5])
    ylabel('% from total sequence')
    title(term{1})
end

saveas(gcf, 'aa_comp_Cterm.fig')


%% Load samples
clear all

load('targetFresh.mat')
GP=load('group_imp.mat');
load('20231005_medians_BH_RM.mat','medianSumPromNewAll')     
%% Colors and targets 

cluster_colors=[0.913725490196078,0.596078431372549,0.478431372549020;0.400000000000000,0.760784313725490,0.647058823529412;0.196078431372549,0.533333333333333,0.741176470588235;0.721568627450980,0.670588235294118,0.827450980392157];
color_scheme = brewermap(11, 'Spectral');
compare_colors=color_scheme([8 11],:);


color_scheme = brewermap(20, 'Spectral');
mut_colors = color_scheme([14 17 19],:);

%% Scatter plot of promoters (Figure 4A, bottom)
cmpStrains={'GLN3_350C_WT','GLN3_d101_280_d381_730';'GLN3_orig','GLN3_150_250C_WT'}
figure('Color',[1 1 1])
subplot(1,1,1)
for j=1:2
    subplot(1,2,j)
    for i=1:4
        scatter(log2(medianSumPromNewAll.(cmpStrains{j,1})(targets.geneId(ismember(targets.cluster,i)))+700),log2(medianSumPromNewAll.(cmpStrains{j,2})(targets.geneId(ismember(targets.cluster,i)))+700),100,cluster_colors(targets.cluster(ismember(targets.cluster,i)),:),'filled','MarkerEdgeColor',[1 1 1].*.35)
        hold on
    end
    xlabel(cmpStrains{j,1})
    ylabel(cmpStrains{j,2})
    plot([10 19], [10 19], '--', 'color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off')
    axis tight
end
set(gcf,'renderer','painters')
saveas(gcf, 'fig4_scatter.svg');

%% Scatter plot of correlation to with and without IDR (Figure 4B, D, E, G)

currStrains481 = {'GLN3_orig';'GLN3_dFWY481_580';'GLN3_dS481_580';'GLN3_dM481_580';'GLN3_dPG481_580';'GLN3_dLIV481_580';'GLN3_dKR481_580';
                                    'GLN3_edPG481_580';'GLN3_edLIV481_580';'GLN3_edNQT481_580';'GLN3_edKR481_580';
                                    'GLN3_150_250C_WT'};
                                
currStrainsPoint = {'GLN3_orig'; 'GLN3_d482_484';'GLN3_d503_505';'GLN3_IPKLP';'GLN3_d537_541';'GLN3_d569_571';'GLN3_150_250C_WT'};

currStrains101_280 = {'GLN3_350C_WT'; 'GLN3_dKR101_280_d381_730'; 'GLN3_dS101_280_d381_730';  'GLN3_dPG101_280_d381_730'; 'GLN3_dM101_280_d381_730';
                                         'GLN3_dFWY101_280_d381_730'; 'GLN3_dLIV101_280_d381_730'; 'GLN3_dNQT101_280_d381_730';
                                         'GLN3_edPG101_280_d381_730'; 'GLN3_edLIV101_280_d381_730'; 'GLN3_edNQT101_280_d381_730'; 'GLN3_edFYW101_280_d381_730';
                                         'GLN3_d101_280_d381_730'};
                  
currStrainsChimera = {'GLN3_orig'; 'GLN3_gla_wt'; 'GLN3_pCglab494_601'; 'GLN3_pCglab602_690'; 'GLN3_pCglab602_823'; 'GLN3_pCglab494_823';
                                    'GLN3_150_250C_WT'};
all_colorsCt = [3 1 1 1 1 1 1 2 2 2 2 3]; %c-terminal
all_colors101 = [3 1 1 1 1 1 1 1 2 2 2 2 3 ]; %n-terminal
all_colorsChim = [2 1 1 1 1 2]; %chimera
all_colorsPoint = [2 1 1 1 1 1 2]; %point

allStrains={currStrains101_280,currStrains481,currStrainsPoint,currStrainsChimera};
allColors={all_colors101,all_colorsCt,all_colorsPoint,all_colorsChim};
for g=1:numel(allStrains)
    currStrains = allStrains{g};
    clear medianSumPromShort
    for i = 1:length(currStrains)
        medianSumPromShort(:,i) = medianSumPromNewAll.(currStrains{i});
    end
    clearvars i

    currStrainsLabels101 = strrep(currStrains, 'GLN3_250N350C_WT', 'DBD');
    currStrainsLabels101 = strrep(currStrainsLabels101, 'GLN3_', '');
    currStrainsLabels101 = strrep(currStrainsLabels101, '101_280', ' ');
    currStrainsLabels101 = strrep(currStrainsLabels101, 'd381_730', ' ');
    currStrainsLabels101 = strrep(currStrainsLabels101, 'ed', '');
    currStrainsLabels101 = strrep(currStrainsLabels101, 'd', '');
    currStrainsLabels101 = strrep(currStrainsLabels101, '_', ' ');



    medianSumPromTarg = medianSumPromShort(targets.geneId, :);

    corr_toWT= corr(medianSumPromTarg,medianSumPromTarg(:,1), 'rows', 'pairwise');
    corr_toNterm = corr(medianSumPromTarg,medianSumPromTarg(:,end), 'rows', 'pairwise');
    subplot(2,2,g)
    hold off
    currColors = allColors{g};
    scatter(corr_toWT(currColors==1), corr_toNterm(currColors==1), 180, mut_colors(currColors(currColors==1), :), 'filled', 'MarkerEdgeColor',[1 1 1].*.35,'DisplayName','deletion')
    hold on
    scatter(corr_toWT(currColors==2), corr_toNterm(currColors==2), 180, mut_colors(currColors(currColors==2), :), 'filled','MarkerEdgeColor',[1 1 1].*.35, 'DisplayName','spread')
    hold on
    scatter(corr_toWT(currColors==3), corr_toNterm(currColors==3), 180, mut_colors(currColors(currColors==3), :), 'filled','MarkerEdgeColor',[1 1 1].*.35, 'DisplayName','extremes')
    text(corr_toWT, corr_toNterm, currStrainsLabels101)
    xlabel("correlation to with IDR")
    ylabel("correlation to without IDR")
    legend('Location', 'northeast')
    axis tight
end
saveas(gcf, '4b.svg');


%% Comparison plot: d/ed (Figure 4C)
clearvars -except medianSumPromNewAll targets cluster_colors
currStrainsCompare = {'GLN3_orig'; 'GLN3_350C_WT'; 
                                         'GLN3_edLIV481_580'; 'GLN3_dLIV481_580'; 'GLN3_edPG481_580'; 'GLN3_dPG481_580';
                                          'GLN3_edKR481_580'; 'GLN3_dKR481_580';
                                           'GLN3_edFYW101_280_d381_730'; 'GLN3_dFWY101_280_d381_730'; 'GLN3_edLIV101_280_d381_730'; 'GLN3_dLIV101_280_d381_730';
                                           'GLN3_edNQT101_280_d381_730'; 'GLN3_dNQT101_280_d381_730'; 'GLN3_edPG101_280_d381_730'; 'GLN3_dPG101_280_d381_730';
                                           'GLN3_dFWY481_580';'GLN3_dS481_580';'GLN3_dM481_580';
                                           'GLN3_dKR101_280_d381_730'; 'GLN3_dS101_280_d381_730';  'GLN3_dM101_280_d381_730';};

for i = 1:length(currStrainsCompare)
    medianSumPromShort(:,i) = medianSumPromNewAll.(currStrainsCompare{i});
end

medianSumPromTarg = medianSumPromShort(targets.geneId, :);

currStrainsLabelsCompare = {'LIV'; 'PG'; 'KR'; 'FWY'; 'LIV'; 'NQT'; 'PG'; 'FWY'; 'S'; 'M'; 'KR'; 'S'; 'M'};
all_colorsComp = [2 2 2 3 3 3 3 2 2 2 3 3 3]; %compare
currColors=all_colorsComp;

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
scatter(1- corr_delAll(currColors==2), ratioAll(currColors==2), 180, cluster_colors(currColors(currColors==2), :), 'filled', 'MarkerEdgeColor',[1 1 1].*.35, 'DisplayName','R2: 481-580')
hold on
scatter(1- corr_delAll(currColors==3), ratioAll(currColors==3), 180, cluster_colors(currColors(currColors==3), :), 'filled', 'MarkerEdgeColor',[1 1 1].*.35, 'DisplayName','R3: 101-280')
hold on
text(1-corr_delAll, ratioAll, currStrainsLabelsCompare)
plot([0 1], [0 0], '--', 'color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off')
hold on
plot([0 1], [1 1], '--', 'color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off')
xlabel("1 - deletion corr to WT")
ylabel("relative effect of the spread")
legend('Location', 'northeast')
saveas(gcf, '4c.svg');
%% Figure 4F,G, S3B
clear all
load('targetFresh.mat')
GP=load('group_imp.mat');
repeatData=load('repeatData.mat')

cluster_colors=[0.913725490196078,0.596078431372549,0.478431372549020;0.400000000000000,0.760784313725490,0.647058823529412;0.196078431372549,0.533333333333333,0.741176470588235;0.721568627450980,0.670588235294118,0.827450980392157];
color_scheme = brewermap(20, 'Spectral');
mut_colors = color_scheme([14 17 19],:);


                                
currStrainsPoint = {'GLN3_orig'; 'GLN3_d482_484';'GLN3_d503_505';'GLN3_IPKLP';'GLN3_d537_541';'GLN3_d569_571';'GLN3_150_250C_WT'};

currStrainsChimera = {'GLN3_orig'; 'GLN3_pCglab494_601'; 'GLN3_pCglab602_690'; 'GLN3_pCglab602_823'; 'GLN3_pCglab494_823';
                                    'GLN3_150_250C_WT'};

currStrains481 = {'GLN3_orig'; 'GLN3_dM481_580'; 'GLN3_dS481_580'; 'GLN3_dFWY481_580';'GLN3_dPG481_580';'GLN3_edPG481_580';
                                    'GLN3_edNQT481_580'; 'GLN3_dLIV481_580';  'GLN3_edLIV481_580'; 'GLN3_dKR481_580'; 'GLN3_edKR481_580';
                                    'GLN3_150_250C_WT'};
                 
currStrains101_280 = {'GLN3_350C_WT'; 'GLN3_dKR101_280_d381_730'; 'GLN3_dS101_280_d381_730';  'GLN3_dM101_280_d381_730'; 'GLN3_dPG101_280_d381_730';
                                       'GLN3_edPG101_280_d381_730'; 'GLN3_dFWY101_280_d381_730'; 'GLN3_edFYW101_280_d381_730';
                                       'GLN3_dLIV101_280_d381_730'; 'GLN3_edLIV101_280_d381_730';  'GLN3_dNQT101_280_d381_730'; 'GLN3_edNQT101_280_d381_730'; 
                                       'GLN3_d55_87_d381_730';'GLN3_d101_280_d381_730'};
                 


all_colorsPoint = [2 1 1 1 1 1  2]; %motif
all_colorsChim = [2 1 1 1 1 2]; %chimera
all_colorsCt = [3 1 1 1 1 2 2 1 2 1 2 3]; %c-terminal
all_colors101 = [3 1 1 1 1 2 1 2 1 2 1 2 3 3]; %n-terminal

allStrains={currStrains101_280,currStrains481,currStrainsPoint,currStrainsChimera};
allColors={all_colors101,all_colorsCt,all_colorsPoint,all_colorsChim};
intClusters=[{3},{2},{2},{2}]
c=0

for g=1:4
    [~,repeatData.samples.strainOrder]=ismember(repeatData.samples.strain,allStrains{g})
    repeats=repeatData.samples(repeatData.samples.strainOrder>0,:)
    sumPromShort=repeatData.sumPromRep(:,repeatData.samples.strainOrder>0)
    [repeats,idx]=sortrows(repeats,"strainOrder")
    sumPromShort=sumPromShort(:,idx);
    clearvars idx

    [currStrainsUnique, ~, idxUn] = unique(repeats.strain, 'stable');

    currStrainsLabels101 = strrep(currStrainsUnique, 'GLN3_250N350C_WT', 'DBD');
    currStrainsLabels101 = strrep(currStrainsLabels101, 'GLN3_', '');
    currStrainsLabels101 = strrep(currStrainsLabels101, '101_280', ' ');
    currStrainsLabels101 = strrep(currStrainsLabels101, 'd381_730', ' ');
    currStrainsLabels101 = strrep(currStrainsLabels101, 'ed', '');
    currStrainsLabels101 = strrep(currStrainsLabels101, 'd', '');
    currStrainsLabels101 = strrep(currStrainsLabels101, '_', ' ');
    currColors = allColors{g};
    for cluster = intClusters{g}
        c=c+1;
        sumBindingCluster = sum(sumPromShort(targets.geneId(targets.cluster==cluster), :), 'omitnan')./sum(sumPromShort(targets.geneId, :), 'omitnan');
        stdStrain = accumarray(idxUn, sumBindingCluster, [], @(x)std(x, 'omitnan'));
        meanStrain = accumarray(idxUn, sumBindingCluster, [], @(x)mean(x, 'omitnan'));
        subplot(2,2,c)
        hold off
        barh(1:numel(meanStrain), meanStrain.*(currColors==1)', 'Barwidth',0.5, 'FaceColor', mut_colors(1, :),'Barwidth',0.8,'Linestyle','none')
        hold on
        errorbar(meanStrain(currColors==1), find(currColors==1), stdStrain(currColors==1), '.','horizontal', 'Color', 'k')
        barh(1:numel(meanStrain), meanStrain.*(currColors==2)', 'Barwidth',0.5, 'FaceColor', mut_colors(2, :),'Barwidth',0.8,'Linestyle','none')
        errorbar(meanStrain(currColors==2), find(currColors==2), stdStrain(currColors==2), '.','horizontal', 'Color', 'k')

        barh(1:numel(meanStrain), meanStrain.*(currColors==3)', 'Barwidth',0.5, 'FaceColor', mut_colors(3, :),'Barwidth',0.8,'Linestyle','none')
        errorbar(meanStrain(currColors==3), find(currColors==3), stdStrain(currColors==3), '.','horizontal', 'Color', 'k')

        %xlabel("fraction of binding to cluster "+  cluster)
        yticklabels(strrep(currStrainsLabels101,'_',' '))
        set(gca,'YDir','reverse'); % Flips the Y Axis
        ylim([0 numel(meanStrain)]+.5)
        clearvars sumBindingCluster meanBindingCluster stdStrain meanStrain
    end
end
saveas(gcf, 'Fig4BarsR.svg')


%% ST/DE in phospho region (Figure 4G, top)
Scer = 'YSSSFMAASLQQLHEQQQVDVNSNTNTNSNRQNWNSSNSVSTNSRSSNFVSQKPNFDIFNTPVDSPSVSRPSSRKSHTSLLS';
Cglab = 'YSSSFINTNPQQTQDIGNFNDDQSIQGNGTGVNSQNIRRINSNYDSPQPNFDLFRLDSNKDSPENIPDVLRSDSRLSQKSQVSHTSLLS';
ST = 'ST';
DE = 'DE';

imagesc(ismember(aa2int(Cglab),aa2int(ST)))
saveas(gcf, 'composST_glab.svg');

%% Figure S3A,C
clear all
load('20231005_medians_BH_RM.mat','medianSumPromNewAll')
load('targetFresh.mat')
cluster_colors=[0.913725490196078,0.596078431372549,0.478431372549020;0.400000000000000,0.760784313725490,0.647058823529412;0.196078431372549,0.533333333333333,0.741176470588235;0.721568627450980,0.670588235294118,0.827450980392157];

cmpStrains={'GLN3_orig','GLN3_050C_WT';'GLN3_150_250C_WT','GLN3_d503_505'}
figure('Color',[1 1 1])
subplot(1,1,1)
for j=1:2
    subplot(1,2,j)
    for i=1:4
        scatter(log2(medianSumPromNewAll.(cmpStrains{j,1})(targets.geneId(ismember(targets.cluster,i)))+700),log2(medianSumPromNewAll.(cmpStrains{j,2})(targets.geneId(ismember(targets.cluster,i)))+700),100,cluster_colors(targets.cluster(ismember(targets.cluster,i)),:),'filled','MarkerEdgeColor',[1 1 1].*.35)
        hold on
    end
    xlabel(cmpStrains{j,1})
    ylabel(cmpStrains{j,2})
    plot([10 19], [10 19], '--', 'color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off')
    axis tight
end