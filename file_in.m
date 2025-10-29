clear all
%% 
% 结构导入

opts = spreadsheetImportOptions("NumVariables", 4);
opts.Sheet = "河道结构";
opts.DataRange = "A1:D95";opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4"];
opts.VariableTypes = ["double", "double", "double", "double"];structure = readtable("I:\基于SWAT的泾河BMPs效果评估与河道管理\rt\新建文件夹\参数样本库.xlsx", opts, "UseExcel", false);
structure = table2array(structure);
clear  opts
%% 
% 输入数据导入

opts = spreadsheetImportOptions("NumVariables", 16);
opts.Sheet = "Sheet2";
opts.DataRange = "A2:P420769";
opts.VariableNames = ["SUB", "PETmm", "WYLD", "SYLD", "NO2", "CHLA", "CBOD", "DISOX", "ORGN", "ORGP", "NO3", "NH4", "MINP", "TMVAP","PCP", "RA"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
inputs = readtable("I:\基于SWAT的泾河BMPs效果评估与河道管理\rt\新建文件夹\输入数据库.xlsx", opts, "UseExcel", false);
inputs = table2array(inputs);
clear  opts
%% 
% 参数导入

opts = spreadsheetImportOptions("NumVariables", 1);
opts.Sheet = "bsn";
opts.DataRange = "B1:B24";
opts.VariableNames = "VarName2";
opts.VariableTypes = "double";
paramater_bsn = readtable("I:\基于SWAT的泾河BMPs效果评估与河道管理\rt\新建文件夹\参数样本库.xlsx", opts, "UseExcel", false);
paramater_bsn = table2array(paramater_bsn);
clear  opts


opts = spreadsheetImportOptions("NumVariables", 25);
opts.Sheet = "rte";
opts.DataRange = "B2:Z49";
opts.VariableNames = ["chd", "chn2", "chs2", "chw2", "chl2", "chk2", "ch_cov2", "sub_ha", "ch_revap", "alpha_bank", "ch_onco", "ch_opco", "bc1", "bc2", "latitude", "rs1", "rk1", "rk3", "rk4", "bc3", "bc4", "rs2", "rs3", "rs4", "rs5"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
paratemer_rte = readtable("I:\基于SWAT的泾河BMPs效果评估与河道管理\rt\新建文件夹\参数样本库.xlsx", opts, "UseExcel", false);
paratemer_rte = table2array(paratemer_rte);
clear  opts

opts = spreadsheetImportOptions("NumVariables", 12);
opts.Sheet = "rte";
opts.DataRange = "AB2:AM49";
opts.VariableNames = ["cherodmo1", "cherodmo2", "cherodmo3", "cherodmo4", "cherodmo5", "cherodmo6", "cherodmo7", "cherodmo8", "cherodmo9", "cherodmo10", "cherodmo11", "cherodmo12"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
ch_erodmo = readtable("I:\基于SWAT的泾河BMPs效果评估与河道管理\rt\新建文件夹\参数样本库.xlsx", opts, "UseExcel", false);
ch_erodmo = table2array(ch_erodmo);
clear  opts


rchstor=zeros(1,48);
bankst=zeros(1,48);
sedst=zeros(1,48);
ch_orgp=zeros(1,48);
rch_cbod=zeros(1,48);
alage=zeros(1,48);
organicn=zeros(1,48);
nitriten=zeros(1,48);
ammonian=zeros(1,48);
nitraten=zeros(1,48);
organicp=zeros(1,48);
dislovp=zeros(1,48);
rch_dox=zeros(1,48);
depch=zeros(1,48);