% set the paths of datafile and glmnet package accordingly
fpath = '.\';% fpath = 'set the file path';
dpath = strcat(fpath,'example_data.\');
glmnetpath = strcat(fpath,'glmnet_matlab\');

% import input data
lfc = dlmread(strcat(dpath,'lfc_mouse-pancreas-beta_13010genesX87samples.txt'),'\t');

tobject = readtable(strcat(dpath,'table_of_samples.txt'),'delimiter','\t');
grplist = tobject.group;
tp = tobject.time;
GList = readtable(strcat(dpath,'list_of_genes.txt'),'ReadVariableNames',false);
GList = GList.Var1;

% import TF-gene and protein-protein interactions data
tftg = readtable(strcat(dpath,'edges-TFTG_mouse_pancreas_fromCellNet.txt'),'delimiter','\t');
tftg = table2cell(tftg);
ppi = readtable(strcat(dpath,'edges_ppi_fromSTRING_short.txt'),'delimiter','\t');
ppi = table2cell(ppi);


%%
% calculate slopes of log2FC for each time point
slope = generateSlope(lfc,tp,grplist);

% generate protein-gene network (PGN)
ppi_thre = 500;
tftg_thre = 0;
ptf_thre = 0;

[pgn,ppn2] = generatePGN(GList,tftg,ppi,ppi_thre,tftg_thre,ptf_thre);

%%
addpath(glmnetpath)

kfold = 10;
par = 1;
[Pscore,A] = salmon(lfc,slope,pgn,grplist,kfold,par);
