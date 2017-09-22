%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WEIGHT: Streamline Volume Density - VD (7).
% PURPOSE: Compute NORMALIZED local Clustering Coefficients for each network
% and compare between clean and unclean.


%%%%%%%%% Info %%%%%%%%%%%
%%% NOS_Variables - Output from Global Analyses
%%% E_local_rand_cl/un -  82x1000x275 Matrix.  DIM1: Nodes, DIM2: the local
%%%                       efficiency of 'jth' random network for the
%%%                       subject specified by DIM3.
%%% NOS_Node_Out_A - Output from local NON-normalized analyses
%%% LocalEfficiency_cl/un - 82x275 Matrix. Each row stores the efficiency
                            coefficient of the 'jth' subject (DIM2).
%%% NOTE: the clustering versions of these have the more thorough code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%--- Add paths
addpath /UndergraduateThesis
addpath /UndergraduateThesis/RUN
addpath /UndergraduateThesis/RUN/'Brain Connectivity Toolbox'/
addpath /UndergraduateThesis/RUN/2_Node_Level_Differences/
addpath /UndergraduateThesis/RUN/2_Node_Level_Differences/Output_Global
addpath /UndergraduateThesis/RUN/2_Node_Level_Differences/Output_Local
addpath /UndergraduateThesis/RUN/2_Node_Level_Differences/Normalized_Local_Script/

%--- Introduce connectivity Matrices
Outliers_Unclean_Out
clc
clearvars -except 'Clean3' 'Unclean4' 'connectivity' 'regionDescriptions' 'ROIs' 'subjects' 'weightDescriptions'
conn_clean= connectivity(:,:,:, Clean3);
conn_unclean= connectivity(:,:,:, Unclean4);
m = size(conn_clean(:,:,7,:),4);
n = size(conn_unclean(:,:,7,:),4);

%--- Load data generated in GLOBAL analyses
load ('VD_Variables',  'E_local_rand_cl', 'E_local_rand_un')

%--- Loading data generated in LOCAL non-normalized anayses
load ('VD_Node_Out_A',  'LocalEfficiency_un', 'LocalEfficiency_cl')



%%%%%%%%%%%  Normalized Local Efficiency  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--- CAT: Matrix with the local EFFICIENCIES for the random
%networks of each subject.
CAT_cl= E_local_rand_cl;
CAT_un= E_local_rand_un;

%--- Collapse across the third DIM
CAT2_cl = mean(CAT_cl,2);
CAT2_un = mean(CAT_un,2);

%--- Take each vector from the 3rd DIM and store each subject next to
% each other in a matrix.
LERand_CL = zeros(82, 275);
for i= 1:size(CAT2_cl,3)
  LERand_CL(:,i)= CAT2_cl(:,1,i);%Local Clust Coeff Random Networks - For Clean
end

LERand_UN = zeros(82, 275);
for i= 1:size(CAT2_un,3)
 LERand_UN(:,i)= CAT2_un(:,1,i);%Local Clust Coeff Random Networks - For Clean
end

%---Normalized local efficiencies
Loc_Effi_Normalized_cl = LocalEfficiency_cl./LERand_CL;
Loc_Effi_Normalized_un = LocalEfficiency_un./LERand_UN;

%making NaNs = 0
Loc_Effi_Normalized_cl(isnan(Loc_Effi_Normalized_cl))=0;
Loc_Effi_Normalized_un(isnan(Loc_Effi_Normalized_un))=0;


%------------ ------------ ttests------------ ------------

%%%%%%%%%%% WITH bonferroni %%%%%%%%%%%
LEffi_NORM_Comparison = zeros(82, 9);
for i = 1: 82
    %bonferroni correction - P= 0.05/82
[H,P,CI,STATS]= ttest(Loc_Effi_Normalized_cl(i,:), Loc_Effi_Normalized_un(i,:), 'alpha', 0.00060975609);
LEffi_NORM_Comparison(i,1) =  H;
LEffi_NORM_Comparison(i,2) =  P;
LEffi_NORM_Comparison(i,3) = STATS.tstat;
LEffi_NORM_Comparison(i,4) = mean(Loc_Effi_Normalized_cl(i,:));
LEffi_NORM_Comparison(i,5) = std(Loc_Effi_Normalized_cl(i,:));
LEffi_NORM_Comparison(i,6) = mean(Loc_Effi_Normalized_un(i,:));
LEffi_NORM_Comparison(i,7) = std(Loc_Effi_Normalized_un(i,:));
LEffi_NORM_Comparison(i,8) = STATS.sd;
LEffi_NORM_Comparison(i,9) = i; %gives you the number of the node
end

%%%%%%%%%%% NO BONFERRONI  %%%%%%%%%%%%%%
LEffi_NORM_Comparison_bon = zeros(82, 9);
for i = 1: 82
    %bonferroni correction - P= 0.05/82
[H,P,CI,STATS]= ttest(Loc_Effi_Normalized_cl(i,:), Loc_Effi_Normalized_un(i,:));
LEffi_NORM_Comparison_bon(i,1) =  H;
LEffi_NORM_Comparison_bon(i,2) =  P;
LEffi_NORM_Comparison_bon(i,3) = STATS.tstat;
LEffi_NORM_Comparison_bon(i,4) = mean(Loc_Effi_Normalized_cl(i,:));
LEffi_NORM_Comparison_bon(i,5) = std(Loc_Effi_Normalized_cl(i,:));
LEffi_NORM_Comparison_bon(i,6) = mean(Loc_Effi_Normalized_un(i,:));
LEffi_NORM_Comparison_bon(i,7) = std(Loc_Effi_Normalized_un(i,:));
LEffi_NORM_Comparison_bon(i,8) = STATS.sd;
LEffi_NORM_Comparison_bon(i,9) = i; %gives you the number of the node
end


save('VD_EFFI_Norm_Local_Output')
