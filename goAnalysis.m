function out=goAnalysis(intGenes,varargin)
    ip=inputParser();
    ip.addParameter('BG',[]);
    ip.addParameter('pTh',0.05);
    ip.addParameter('adjust',true);
    ip.addParameter('type',[]);
    ip.parse(varargin{:})
    if numel(ip.Results.type)==0
        load('gotermsYeastAdd.mat','goMap','GOtable')
    elseif ip.Results.type==1
       load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/groupForGO.mat','goMap','GOtable')
    elseif  ip.Results.type==2
       load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/keggForGo.mat','goMap','GOtable')
    elseif  ip.Results.type==3
       load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/group2ForGO.mat','goMap','GOtable')
       goMap=logical(goMap);
    end
    if numel(ip.Results.BG)==0
        BG=find(full(any(goMap,2)));
    else
        BG=ip.Results.BG;
    end
    GOtable.n=full(sum(goMap(intGenes,:)))';
    GOtable.pVals=1-arrayfun(@(b)hygecdf(GOtable.n(b)-1,numel(BG),sum(goMap(BG,b)),numel(intGenes)),1:size(goMap,2))';
    GOtable.genes=arrayfun(@(b)intGenes(goMap(intGenes,b)),1:size(goMap,2),'UniformOutput',false)';
    GOtable.log2Enr=(log2(full(sum(goMap(intGenes,:)))+.1)-log2(numel(intGenes)+0.1)-log2(full(sum(goMap(BG,:)))+0.1)+log2(numel(BG)+0.1))';
    out=sortrows(GOtable(GOtable.n>1&GOtable.pVals<ip.Results.pTh,:),{'n'},'descend');
end

%% create GO table
function createGO()
go2=table()
c=0;
goMap=sparse(6701,sum(arrayfun(@(x)numel(GP.groups{x}{1}),1:26)));
for i=1:numel(GP.groups)
    for j=1:numel(GP.groups{i}{2})
        c=c+1;
        goLine(c,:)=table(i,GP.groups{i}{1}(j),j,'VariableNames',{'Var4','Var5','Var6'});
        goMap(GP.groups{i}{2}{j}(GP.groups{i}{2}{j}>0),c)=true;
    end
end
GOtable=goLine
end