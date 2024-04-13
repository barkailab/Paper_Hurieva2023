%% (i) Load individual repeats
pathToBH = '/home/labs/barkailab/bohdana/001_sequencing_data/003_checProfiles';
strainsBH = dir(pathToBH);
strainsBH = strainsBH(3:end); 
CombineChecBH = readtable('/home/labs/barkailab/bohdana/001_sequencing_data/001_Tables/CombChecAll.xlsx');
newToCombineBH = CombineChecBH.sampleNames;

newToCombineBH= strrep(newToCombineBH, '.mat', '');
for i = 1:length(newToCombineBH)
    rm = regexp(char(newToCombineBH{i}) , '(_\d*$|_[ABCDabcd]$|_rep\d$|_\dF$)', 'split');
    nameN(i) = rm(1);
    clearvars rm
end
strainNamesNew = unique(nameN);
clearvars nameN i

for i =1:length(strainNamesNew) %among all unique names
    clearvars sampleNames sumPromStrain %clear from previous iteration
   a =  ['^' strainNamesNew{i} '(_\d*$|_[ABCDabcd]$|_rep\d$|_\dF$)'];
   sampleNames = newToCombineBH(~cellfun(@isempty, regexp(newToCombineBH,a,'match'))); %select all files with this strain
   clearvars a
    for j = 1:length(sampleNames)
        load([pathToBH '/' sampleNames{j} '.mat'])
        sumPromRepeatAll.(sampleNames{j}) = sumPromNew;
        %normProfileRepeat.(sampleNames{j}) = normProfile;
        clearvars normProfile meta merScore merScoreNew merScore5 pwm pwmNew sumProm sumPromNew
    end  
    NumberOfSamples.(strainNamesNew{i}) = length(sampleNames);
end
clearvars i j sampleNames

% filter out subtelomere genes
GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');
subtelomereGenes = GP.groups{1, 23}{1, 2}{1, 61};
fldNames = fieldnames(sumPromRepeatAll); %get field names from the struct
for i = 1:length(fldNames)
    sumPromRepeat.(fldNames{i}) = sumPromRepeatAll.(fldNames{i});
    sumPromRepeat.(fldNames{i}) (subtelomereGenes) = nan;
end
strainNamesNewRep = fldNames;
clearvars fldNames i sumPromRepeatAll nucReads totalReads


%% Strains order for plots                                
currStrainsPoint = {'GLN3_orig'; 'GLN3_d482_484';'GLN3_d503_505';'GLN3_d537_541';'GLN3_d569_571';'GLN3_150_250C_WT'};

currStrainsChimera = {'GLN3_orig'; 'GLN3_pCglab494_601'; 'GLN3_pCglab602_690'; 'GLN3_pCglab602_823'; 'GLN3_pCglab494_823';
                                    'GLN3_150_250C_WT'};

currStrains481 = {'GLN3_orig'; 'GLN3_dM481_580'; 'GLN3_dS481_580'; 'GLN3_dFWY481_580';'GLN3_dPG481_580';'GLN3_edPG481_580';
                                    'GLN3_edNQT481_580'; 'GLN3_dLIV481_580';  'GLN3_edLIV481_580'; 'GLN3_dKR481_580'; 'GLN3_edKR481_580';
                                    'GLN3_150_250C_WT'};
                 
currStrains101_280 = {'GLN3_350C_WT'; 'GLN3_dKR101_280_d381_730'; 'GLN3_dS101_280_d381_730';  'GLN3_dM101_280_d381_730'; 'GLN3_dPG101_280_d381_730';
                                       'GLN3_edPG101_280_d381_730'; 'GLN3_dFWY101_280_d381_730'; 'GLN3_edFYW101_280_d381_730';
                                       'GLN3_dLIV101_280_d381_730'; 'GLN3_edLIV101_280_d381_730';  'GLN3_dNQT101_280_d381_730'; 'GLN3_edNQT101_280_d381_730'; 
                                       'GLN3_d55_87_d381_730';'GLN3_d101_280_d381_730'};
                 

%% (i) Get good repeats
count = 1;
for i =1:length(currStrains101_2) %among all unique names
    clearvars sampleNames
   a =  ['^' currStrains101_2{i} '(_\d*$|_[ABCDabcd]$|_rep\d$|_\dF$)'];
   sampleNames = strainNamesNewRep(~cellfun(@isempty, regexp(strainNamesNewRep,a,'match'))); %select all files with this strain
   clearvars a
    for j = 1:length(sampleNames)
        currStrainsRepeat{count} = sampleNames{j};
        count = count+1;
    end  
end
clearvars i j count sampleNames ans

strainNameRep = regexp(currStrainsRepeat,'.*(?=_.$)','match','once');
[currStrainsUnique, ~, idxUn] = unique(strainNameRep, 'stable');


%% Combine indiv repeats
for i = 1:length(currStrainsRepeat)
       sumPromShort(:,i) = sumPromRepeat.(currStrainsRepeat{i}); 
end
clearvars i


%% Labels
currStrainsLabels101 = strrep(currStrainsUnique, 'GLN3_250N350C_WT', 'DBD');
currStrainsLabels101 = strrep(currStrainsLabels101, 'GLN3_', '');
currStrainsLabels101 = strrep(currStrainsLabels101, '101_280', ' ');
currStrainsLabels101 = strrep(currStrainsLabels101, 'd381_730', ' ');
currStrainsLabels101 = strrep(currStrainsLabels101, 'ed', '');
currStrainsLabels101 = strrep(currStrainsLabels101, 'd', '');
currStrainsLabels101 = strrep(currStrainsLabels101, '_', ' ');
currStrainsLabels101{14} = 'd101 280';


%% Colors and targets 
addpath('/home/labs/barkailab/bohdana/001_sequencing_data/000_Scripts')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/targetFresh.mat')

color_scheme = brewermap(11, 'Spectral');
cluster_colors=color_scheme([8:11],:);

all_colorsPoint = [2 1 1 1 1 2]; %motif
all_colorsChim = [2 1 1 1 1 2]; %chimera
all_colorsCt = [3 1 1 1 1 2 2 1 2 1 2 3]; %c-terminal
all_colors101 = [3 1 1 1 1 2 1 2 1 2 1 2 3 3]; %n-terminal

color_scheme = brewermap(20, 'Spectral');
mut_colors = color_scheme([14 17 19],:);


%% Bargraph with error bars
figure('Color',[1 1 1])
yAxis = categorical(currStrainsUnique);
yAxis = reordercats(yAxis, currStrainsUnique);
all_colors = all_colorsPoint;
for i = 1:4
    cluster = i;
    sumBindingCluster = sum(sumPromShort(targets.geneId(targets.cluster==cluster), :), 'omitnan')./sum(sumPromShort(targets.geneId, :), 'omitnan');
    stdStrain = accumarray(idxUn, sumBindingCluster, [], @(x)std(x, 'omitnan'));
    meanStrain = accumarray(idxUn, sumBindingCluster, [], @(x)mean(x, 'omitnan'));
    subplot(2,2,i)
    barh(yAxis(all_colors==1), meanStrain(all_colors==1), 'Barwidth',0.5, 'FaceColor', mut_colors(1, :))
    hold on
    errorbar(meanStrain(all_colors==1), yAxis(all_colors==1), stdStrain(all_colors==1), '.','horizontal', 'Color', 'k')
    hold on
    barh(yAxis(all_colors==2), meanStrain(all_colors==2), 'Barwidth',0.5, 'FaceColor', mut_colors(2, :))
    hold on
    errorbar(meanStrain(all_colors==2), yAxis(all_colors==2), stdStrain(all_colors==2), '.','horizontal', 'Color', 'k')
    hold on
%     barh(yAxis(all_colors==3), meanStrain(all_colors==3), 'Barwidth',0.5, 'FaceColor', mut_colors(3, :))
%     hold on
%     errorbar(meanStrain(all_colors==3), yAxis(all_colors==3), stdStrain(all_colors==3), '.','horizontal', 'Color', 'k')
%     hold on
    xlabel("fraction of binding to cluster "+  cluster)
    yticklabels(currStrainsPoint)
    set(gca,'YDir','reverse'); % Flips the Y Axis
    clearvars sumBindingCluster meanBindingCluster stdStrain meanStrain
end
saveas(gcf, 'Fig4.fig')


clearvars -except sumPromRepeat strainNamesNewRep NumberOfSamples
